"""Gap-independent protein entry: same protein core, oriented upstream / downstream neighbour genes compared.

Discovery starts from every annotated CDS of every assembly in a subset manifest (columns sample, fasta, gff3, faa),
not from a gap list, reads or an empty site. Evidence is assembly_only; nothing here is an insertion call.

genes    GFF3 CDS + RNA units per contig (s4_annotate.parse_gff), S4 coding check of every CDS against its FAA protein,
         gene order per contig and circularity (region Is_circular=true). Proteins named by full sequence sha256.
cluster  MMseqs2 at the frozen S4 settings (family 0.9/0.8, homology group 0.5/0.8) exactly as stage1_catalog runs
         them, then every member aligned to its final representative; failing members are reassigned to a passing
         representative or become singletons, with the failing values kept (stage1_catalog.reassign_failed).
call     Core occurrences:
           single_gene  one complete CDS, S4 rule: not partial / pseudo / origin-crossing, translation matches FAA
                        (match or match_with_transl_except) and no coding flag (no_start, no_stop, internal_stop, ...);
                        core key = protein family
           two_gene     two consecutive complete CDS on one contig; core key = both families in order with relative
                        orientation, canonical over the two strands; spacing kept per occurrence, not in the key
         Core span = first start to last end (two-gene includes the spacing). Two scopes are compared separately:
           all             every occurrence (wide protein search, no length limit)
           target_lt2500   only occurrences with span < 2500 bp (the 1-2 complete coding unit target subset)
         All core keys with >= 2 occurrences in a scope are kept, including several copies in one assembly. Context = the K genes
         on each side in core orientation (wrapping only on circular contigs), labelled by homology group (CDS with
         protein) or RNA type + product; CDS without a protein carry no label and never count as shared.
         Pair relation (prototype diagnostic, not orthology):
           similar_context              a shared neighbour label on both sides
           different_context_candidate  no shared label on either side, every neighbour of both occurrences labelled
           one_side_shared_unresolved   shared on one side; the other side fully labelled and not shared
           insufficient_context_labels  a side without a shared label has an unlabelled neighbour (unknown != different)
           edge_or_missing_context      fewer than K neighbours on a side of either occurrence
         Context classes: occurrences joined by similar_context links (provisional, not S6 site groups).
"""
import argparse
import csv
import gzip
import hashlib
import itertools
import json
import sys
from argparse import Namespace
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'fna_pool'))
from s4_annotate import MEMBER_MIN, MMSEQS_SETTINGS, coding_check, parse_gff, read_faa  # noqa: E402
from stage1_catalog import cmd_cluster, reassign_failed  # noqa: E402

VERSION = 'core_context_protein_v0.2'
K = 3
CODING_OK = {'match', 'match_with_transl_except'}
RELATIONS = ['similar_context', 'different_context_candidate', 'one_side_shared_unresolved', 'insufficient_context_labels',
             'edge_or_missing_context']
TARGET_MAX = 2500  # exclusive, bp of core span
MAX_PAIRS = 20000
GENE_COLS = ['gene_uid', 'sample', 'contig', 'index', 'start0', 'end0', 'strand', 'type', 'unit_id', 'protein_id', 'protein_sha256',
             'product', 'partial', 'pseudo', 'crosses_origin', 'coding_status', 'coding_flags', 'complete_cds']


def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def read_tsv(path):
    op = gzip.open if str(path).endswith('.gz') else open
    with op(path, 'rt') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def write_tsv(path, rows, cols):
    op = gzip.open if str(path).endswith('.gz') else open
    with op(path, 'wt') as f:
        w = csv.DictWriter(f, cols, delimiter='\t', lineterminator='\n', extrasaction='ignore')
        w.writeheader(); w.writerows(rows)


def read_fna(path):
    out, name, buf = {}, None, []
    for line in open(path):
        if line.startswith('>'):
            if name:
                out[name] = ''.join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name:
        out[name] = ''.join(buf)
    return out


def is_complete(u, status, flags, has_protein):
    """S4 rule for a complete CDS; any coding flag (e.g. no_start) or an origin-crossing unit disqualifies it."""
    return (u['type'] == 'CDS' and has_protein and not u['partial'] and not u['pseudo'] and not u['crosses_origin']
            and status in CODING_OK and not flags)


def cmd_genes(a):
    a.work.mkdir(parents=True, exist_ok=False)
    rows, proteins, inputs = [], {}, []
    for m in read_tsv(a.manifest):
        files = {k: m[k] for k in ('fasta', 'gff3', 'faa')}
        for k, p in files.items():
            want = m.get({'fasta': 'fna_sha256', 'gff3': 'gff3_sha256', 'faa': 'faa_sha256'}[k])
            got = digest(p)
            if want and want != got:
                raise ValueError('input hash differs from manifest: ' + p)
            inputs.append(dict(sample=m['sample'], kind=k, path=p, sha256=got))
        seqs, faa = read_fna(m['fasta']), read_faa(m['faa'])
        circular = {l.split('\t')[0] for l in open(m['gff3']) if '\tregion\t' in l and 'Is_circular=true' in l}
        units = parse_gff(m['gff3'], {c: len(s) for c, s in seqs.items()})
        for contig in sorted(units):
            us = sorted(units[contig], key=lambda u: (u['segments'][0][0], u['segments'][-1][1], u['id']))
            for i, u in enumerate(us):
                status, flags = ('rna', []) if u['type'] != 'CDS' else coding_check(u, seqs[contig], faa)
                p = faa.get(u['protein']) if u['type'] == 'CDS' else None
                h = hashlib.sha256(p.encode()).hexdigest() if p else ''
                if p:
                    proteins[h] = p
                rows.append(dict(gene_uid='%s|%s|%d' % (m['sample'], contig, i), sample=m['sample'], contig=contig, index=i,
                    start0=u['segments'][0][0], end0=u['segments'][-1][1], strand=u['strand'], type=u['type'], unit_id=u['id'],
                    protein_id=u['protein'], protein_sha256=h, product=u['product'], partial=u['partial'], pseudo=u['pseudo'],
                    crosses_origin=u['crosses_origin'], coding_status=status, coding_flags=','.join(flags),
                    complete_cds=is_complete(u, status, flags, bool(h)),
                    circular_contig=contig in circular, contig_length=len(seqs[contig])))
    write_tsv(a.work / 'genes.tsv.gz', rows, GENE_COLS + ['circular_contig', 'contig_length'])
    with open(a.work / 'proteins.faa', 'w') as f:
        for h in sorted(proteins):
            f.write('>%s\n%s\n' % (h, proteins[h]))
    (a.work / 'genes_run.json').write_text(json.dumps(dict(version=VERSION, manifest=str(a.manifest),
        manifest_sha256=digest(a.manifest), inputs=inputs, genes=len(rows), proteins=len(proteins),
        coding_status=Counter(r['coding_status'] for r in rows)), indent=1) + '\n')


def cmd_cluster_entry(a):
    cmd_cluster(Namespace(work=a.work, mmseqs=a.mmseqs, threads=a.threads))


def assignments(a):
    """{level: {protein sha256: group id}} after the member-to-final-representative check, plus repair records."""
    w, assign, fixes, checks = a.work / 'mmseqs', {}, [], Counter()
    prefix = {'family': 'PF_', 'homology_group': 'PH_'}
    for level in MMSEQS_SETTINGS:
        reps = {}
        for line in open(w / (level + '_cluster.tsv')):
            rep, mem = line.split()
            reps[mem] = rep
        aln = {}
        for line in open(w / (level + '_members.m8')):
            q, t, fid, qc, tc, _ = line.rstrip('\n').split('\t')
            if (t, q) not in aln or float(fid) > aln[(t, q)][0]:
                aln[(t, q)] = (float(fid), float(tc), float(qc))  # query = representative, target = member
        mid, mc = MEMBER_MIN[level]
        check = {}
        for m, r in reps.items():
            fid, mcov, rcov = (1.0, 1.0, 1.0) if m == r else aln.get((m, r), (0.0, 0.0, 0.0))
            check[(level, m)] = (r, fid, mcov, rcov, fid >= mid and mcov >= mc and rcov >= mc)
            checks[(level, check[(level, m)][4])] += 1
        for m, rec in reassign_failed(Namespace(work=a.work, mmseqs=a.mmseqs, threads=a.threads), level, reps, check).items():
            fixes.append(rec)
            reps[m] = rec['new_rep']
        assign[level] = {m: prefix[level] + r[:12] for m, r in reps.items()}
        if len({r for r in reps.values()}) != len(set(assign[level].values())):
            raise ValueError('representative prefix collision at ' + level)
    return assign, fixes, {'%s_%s' % (lv, 'pass' if ok else 'fail_before_repair'): n for (lv, ok), n in checks.items()}


def label(g, hg):
    if g['type'] != 'CDS':
        return '%s:%s' % (g['type'], g['product'])
    return hg.get(g['protein_sha256']) or None


def context(genes, lo, hi, forward, labels):
    """K neighbours on each side of genes[lo..hi] in core orientation: lists of (label, same direction as core, gene_uid)."""
    n, circ = len(genes), genes[0]['circular_contig'] == 'True'

    def side(step, start):
        out, j = [], start
        while len(out) < K:
            j += step
            if not 0 <= j < n:
                if not circ:
                    break
                j %= n
            if lo <= j <= hi:
                break
            g = genes[j]
            out.append((labels[g['gene_uid']], (g['strand'] == '+') == forward, g['gene_uid']))
        return out
    left, right = side(-1, lo), side(+1, hi)
    return (left, right) if forward else (right, left)


def pair_key(f1, s1, f2, s2):
    """Canonical two-gene key over both strand readings, and whether the + strand reading is the canonical one."""
    flip = {'+': '-', '-': '+'}
    fwd, rev = '%s%s,%s%s' % (f1, s1, f2, s2), '%s%s,%s%s' % (f2, flip[s2], f1, flip[s1])
    return min(fwd, rev), fwd <= rev


def pair_span(g, h):
    """Span of two consecutive genes from the first start to the last end, so the spacing between them counts."""
    return max(int(g['end0']), int(h['end0'])) - min(int(g['start0']), int(h['start0']))


def relation(a, b):
    if min(len(a['up']), len(a['down']), len(b['up']), len(b['down'])) < K:
        return 'edge_or_missing_context'
    shared = [bool({x[0] for x in a[s] if x[0]} & {x[0] for x in b[s] if x[0]}) for s in ('up', 'down')]
    if all(shared):
        return 'similar_context'
    labelled = [all(x[0] for x in a[s] + b[s]) for s in ('up', 'down')]
    if any(not sh and not lb for sh, lb in zip(shared, labelled)):
        return 'insufficient_context_labels'
    if not any(shared):
        return 'different_context_candidate'
    return 'one_side_shared_unresolved'


def compare(groups, scope, keep):
    """Pair relations and provisional context classes of every core key with >= 2 kept occurrences in this scope."""
    pairs, cores = [], []
    for (structure, key), xs in sorted(groups.items()):
        xs = [o for o in xs if keep(o)]
        if len(xs) < 2:
            continue
        parent = {o['occurrence_id']: o['occurrence_id'] for o in xs}

        def root(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]; x = parent[x]
            return x
        total = len(xs) * (len(xs) - 1) // 2
        st, diff_links = Counter(), []
        for p, q in itertools.islice(itertools.combinations(xs, 2), MAX_PAIRS):
            r = relation(p, q)
            st[r] += 1
            pairs.append(dict(scope=scope, structure=structure, core_key=key, occurrence_a=p['occurrence_id'],
                              occurrence_b=q['occurrence_id'], same_sample=p['sample'] == q['sample'], relation=r))
            if r == 'similar_context':
                parent[root(p['occurrence_id'])] = root(q['occurrence_id'])
            elif r == 'different_context_candidate':
                diff_links.append((p['occurrence_id'], q['occurrence_id']))
        for o in xs:
            o['context_class_' + scope] = root(o['occurrence_id'])
        per_sample = Counter(o['sample'] for o in xs)
        spans = [o['span'] for o in xs]
        cores.append(dict(scope=scope, structure=structure, core_key=key, occurrences=len(xs), samples=len(per_sample),
            max_copies_in_one_sample=max(per_sample.values()), span_range='%d..%d' % (min(spans), max(spans)),
            context_classes_provisional=len({root(o['occurrence_id']) for o in xs}),
            context_classes_linked_by_different_candidate=len({root(x) for l in diff_links for x in l}),
            pairs_total=total, pairs_checked=sum(st.values()), pair_review_complete=total <= MAX_PAIRS,
            **{r: st[r] for r in RELATIONS}, products=Counter(o['products'] for o in xs).most_common(1)[0][0][:200],
            spacing_range='' if structure == 'single_gene' else '%d..%d' % (min(o['spacing'] for o in xs), max(o['spacing'] for o in xs))))
    return pairs, cores


def cmd_call(a):
    a.out.mkdir(parents=True, exist_ok=True)
    genes = read_tsv(a.work / 'genes.tsv.gz')
    assign, fixes, checks = assignments(a)
    fam, hg = assign['family'], assign['homology_group']
    labels = {g['gene_uid']: label(g, hg) for g in genes}
    by_contig = defaultdict(list)
    for g in genes:
        by_contig[(g['sample'], g['contig'])].append(g)
    occ = []
    for (sample, contig), gs in sorted(by_contig.items()):
        gs.sort(key=lambda g: int(g['index']))
        for i, g in enumerate(gs):
            if g['complete_cds'] != 'True':
                continue
            f = fam[g['protein_sha256']]
            up, down = context(gs, i, i, g['strand'] == '+', labels)
            occ.append(dict(structure='single_gene', core_key=f, sample=sample, contig=contig, genes=g['gene_uid'],
                start0=g['start0'], end0=g['end0'], span=int(g['end0']) - int(g['start0']), strand=g['strand'], spacing='', homology_groups=hg[g['protein_sha256']],
                products=g['product'], up=up, down=down))
            if i + 1 < len(gs) and gs[i + 1]['complete_cds'] == 'True':
                h = gs[i + 1]
                key, forward = pair_key(f, g['strand'], fam[h['protein_sha256']], h['strand'])
                up2, down2 = context(gs, i, i + 1, forward, labels)
                occ.append(dict(structure='two_gene', core_key=key, sample=sample, contig=contig,
                    genes=g['gene_uid'] + ';' + h['gene_uid'], start0=g['start0'], end0=max(int(g['end0']), int(h['end0'])),
                    span=pair_span(g, h),
                    strand='+' if forward else '-', spacing=int(h['start0']) - int(g['end0']),
                    homology_groups=hg[g['protein_sha256']] + ';' + hg[h['protein_sha256']],
                    products=g['product'] + ' | ' + h['product'], up=up2, down=down2))
    for o in occ:
        o['in_target_length'] = o['span'] < TARGET_MAX
    groups = defaultdict(list)
    for o in occ:
        groups[(o['structure'], o['core_key'])].append(o)
    groups = {k: v for k, v in groups.items() if len(v) >= 2}
    for (structure, key), xs in groups.items():
        for n, o in enumerate(xs):
            o['occurrence_id'] = '%s|%s|%d' % (structure[0], key, n)
    kept = [o for _, xs in sorted(groups.items()) for o in xs]
    pairs, cores = [], []
    for scope, keep in (('all', lambda o: True), ('target_lt2500', lambda o: o['in_target_length'])):
        p, c = compare(groups, scope, keep)
        pairs += p; cores += c
    for o in kept:
        for s in ('up', 'down'):
            o[s + '_labels'] = ','.join('%s%s' % (x[0] or 'nolabel', '' if x[1] else '(rev)') for x in o[s])
            o[s + '_genes'] = ','.join(x[2] for x in o[s])
    occ_cols = ['occurrence_id', 'structure', 'core_key', 'sample', 'contig', 'genes', 'start0', 'end0', 'span', 'in_target_length',
                'strand', 'spacing', 'homology_groups', 'products', 'up_labels', 'down_labels', 'up_genes', 'down_genes',
                'context_class_all', 'context_class_target_lt2500']
    write_tsv(a.out / 'protein_occurrences.tsv.gz', kept, occ_cols)
    write_tsv(a.out / 'protein_pairs.tsv.gz', pairs, ['scope', 'structure', 'core_key', 'occurrence_a', 'occurrence_b', 'same_sample', 'relation'])
    core_cols = ['scope', 'structure', 'core_key', 'occurrences', 'samples', 'max_copies_in_one_sample', 'span_range',
                 'context_classes_provisional',
                 'context_classes_linked_by_different_candidate', 'pairs_total', 'pairs_checked', 'pair_review_complete'] + RELATIONS + \
                ['spacing_range', 'products']
    write_tsv(a.out / 'protein_cores.tsv', cores, core_cols)
    write_tsv(a.out / 'member_check_repairs.tsv', fixes, list(fixes[0]) if fixes else ['level', 'member'])
    stat = {}
    for scope in ('all', 'target_lt2500'):
        for structure in ('single_gene', 'two_gene'):
            cs = [c for c in cores if c['structure'] == structure and c['scope'] == scope]
            ps = [p for p in pairs if p['structure'] == structure and p['scope'] == scope]
            os_ = [o for o in occ if o['structure'] == structure and (scope == 'all' or o['in_target_length'])]
            stat['%s:%s' % (scope, structure)] = dict(occurrences_considered=len(os_),
                core_keys_with_2plus=len(cs), occurrences_kept=sum(c['occurrences'] for c in cs),
                core_keys_multicopy=sum(c['max_copies_in_one_sample'] >= 2 for c in cs),
                core_keys_with_different_context_candidate=sum(c['different_context_candidate'] > 0 for c in cs),
                core_keys_2plus_context_classes=sum(c['context_classes_provisional'] >= 2 for c in cs),
                pairs=len(ps), relations={r: sum(p['relation'] == r for p in ps) for r in RELATIONS},
                different_candidate_pairs_same_sample=sum(p['relation'] == RELATIONS[1] and p['same_sample'] for p in ps),
                truncated_pair_reviews=sum(not c['pair_review_complete'] for c in cs))
    run = json.loads((a.work / 'mmseqs' / 'cluster_run.json').read_text())
    gr = json.loads((a.work / 'genes_run.json').read_text())
    summary = dict(version=VERSION, discovery_entry='all_faa_gff3_proteins', gap_inputs=[], reads_inputs=[],
        evidence_source='assembly_only', K_neighbours=K, genes=gr['genes'], proteins_distinct=gr['proteins'],
        coding_status=gr['coding_status'], complete_cds=sum(g['complete_cds'] == 'True' for g in genes),
        cds_excluded_by_coding_flags=sum(g['type'] == 'CDS' and g['coding_status'] in CODING_OK and g['coding_flags'] != ''
                                         for g in genes),
        cds_excluded_crossing_origin=sum(g['type'] == 'CDS' and g['crosses_origin'] == 'True' for g in genes),
        target_max_span_exclusive=TARGET_MAX,
        member_check=checks, member_repairs=len(fixes), mmseqs_version=run['mmseqs_version'], mmseqs_settings=run['settings'],
        families=len(set(fam.values())), homology_groups=len(set(hg.values())), structures=stat,
        annotation_complete=False, site_groups_computed=False, production_ready=False,
        limitations=['Neighbour homology-group overlap is a context diagnostic, not an orthology or site test.',
                     'Paralogs (e.g. tuf, rRNA-adjacent genes) and rearrangements also give different contexts.',
                     'Single-gene and two-gene cores are not complete insertion boundaries.',
                     'Context classes are provisional and are not S6 site groups.'])
    (a.out / 'protein_summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest='mode', required=True)
    for name in ('genes', 'cluster', 'call'):
        s = sub.add_parser(name)
        s.add_argument('--work', type=Path, required=True)
        if name == 'genes':
            s.add_argument('--manifest', type=Path, required=True)
        if name in ('cluster', 'call'):
            s.add_argument('--mmseqs', default='/global/scratch/users/kh36969/longread_callout/env/s4/bin/mmseqs')
            s.add_argument('--threads', type=int, default=8)
        if name == 'call':
            s.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    {'genes': cmd_genes, 'cluster': cmd_cluster_entry, 'call': cmd_call}[a.mode](a)


if __name__ == '__main__':
    main()
