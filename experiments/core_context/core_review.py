"""Review of the protein-entry U-core-D candidates: sample split, sequence similarity, boundary beyond the CDS, nesting.

Cores reviewed: every target_lt2500 core key (single- and two-gene) with >= 1 different_context_candidate pair in one
subset's committed protein_entry/; only its target-scope occurrences and pairs. Evidence stays assembly_only; nothing
here is an insertion call or a site group.

1 sample split   every different_context_candidate pair is within one assembly (multi-copy) or across assemblies; both
                 are kept and counted separately (within-assembly copies show the pattern but not strain variation).
2 similarity     per pair: protein identity / coverage per core gene in core order (exact sha256, then blastp) and DNA
                 identity / coverage of the oriented coding span (blastn). A protein family (MMseqs 0.9/0.8) is not an
                 identical sequence; the minimum over the two genes is reported for two-gene cores.
3 boundary       per pair: the oriented coding span +/- W bp of both copies aligned (blastn, plus strand). The HSP with the
                 largest overlap of both coding spans is the seed; colinear HSPs within CHAIN_GAP bp are chained to it.
                 The chained HSPs are kept: a side's extension (shared sequence beyond the coding span) is measured
                 only if that core end lies inside one aligned piece in both copies, or the outermost piece stops
                 exactly at it (a measured 0). An end not reached, or in / at an unaligned gap between chained pieces,
                 is unknown, never 0; gap ends are counted separately. Per occurrence and side:
                 resolved only if every different-context partner gives a known extension, none reaches the window /
                 contig end, and they spread <= SPREAD bp; partners not covering the end are counted. Estimated span
                 = coding span + median extension of each side, only when both sides are resolved. The <2500 bp target
                 is re-checked on it. A resolved side is a candidate conserved end, not a mapped insertion boundary. End hints, not calls: identity of the outer 20 bp
                 with the reverse complement of the other end (inverted repeat), longest exact direct repeat (2-12 bp)
                 just outside both ends (target-site duplication; 2-3 bp matches are expected by chance).
                 Cross-assembly pairs also carry the genome ANI of the two assemblies (skani edges, frozen in S2) and
                 core DNA identity minus ANI: a copy that moved recently is more similar than its genomes, a vertically
                 inherited gene in a rearranged neighbourhood is about as similar as its genomes (diagnostic, no cutoff).
5 copy number    complete-CDS target-scope copies of the core per subset assembly (0 where absent) and a pattern:
                 same_count_in_all_assemblies / same_count_where_present / count_varies. Reported only: a constant count
                 does not prove ordinary paralogy and a varying count does not prove movement; nothing is excluded.
4 nesting        a neighbour gene is mobile-like if its product matches a mobile keyword (review_protein.hint:
                 transposase / IS, integrase / recombinase / resolvase, phage) or its protein family is itself a reviewed
                 single-gene core. Host context = the first K neighbours per side after skipping mobile-like genes
                 (scanning at most MAX_SCAN genes). Pairs are re-related on host context with the same rule;
                 near-different but host-similar pairs look like nested or adjacent mobile DNA, not different host sites.
"""
import argparse
import json
import statistics
import subprocess
import tempfile
from argparse import Namespace
from collections import Counter, defaultdict
from pathlib import Path

import protein_context as pc
from review_protein import colour, esc, hint

VERSION = 'core_context_review_v0.4'
W = 3000
CHAIN_GAP = 50
EDGE = 5
SPREAD = 20
MAX_SCAN = 10
MOBILE_HINTS = {'transposase_or_IS', 'integrase_recombinase_resolvase', 'phage_like'}
COMP = str.maketrans('ACGTacgtNn', 'TGCAtgcaNn')
OUTFMT = '6 qseqid sseqid pident length qstart qend sstart send qlen slen bitscore'


def revcomp(s):
    return s.translate(COMP)[::-1]


def cslice(seq, lo, hi):
    """seq[lo:hi] where lo / hi may run past either end of a circular sequence."""
    n, out, i = len(seq), [], lo
    while i < hi:
        j = i % n
        take = min(hi - i, n - j)
        out.append(seq[j:j + take]); i += take
    return ''.join(out)


def window(seq, s, e, forward, w, circular):
    """Coding span [s, e) +/- w bp in core orientation, wrapping only on circular contigs that are long enough.
    Returns the sequence and the core's [start, end) in it."""
    if circular and len(seq) > e - s + 2 * w:
        lo, hi = s - w, e + w
    else:
        lo, hi = max(0, s - w), min(len(seq), e + w)
    sub = cslice(seq, lo, hi)
    return (sub, s - lo, e - lo) if forward else (revcomp(sub), hi - e, hi - s)


def chain(hsps, a0, a1, b0, b1, gap=CHAIN_GAP):
    """Shared block around the core from plus/plus HSPs (1-based inclusive qs, qe, ss, se) of A vs B, cores at [a0, a1)
    and [b0, b1) (0-based). Seed = HSP with the largest overlap of both cores; colinear HSPs within gap bp (small
    overlaps allowed) are added on each side. Returns the chained HSPs (qs0, qe, ss0, se), 0-based half-open and
    ordered along A, or None; the block is their outer extent, the gaps between them are not aligned."""
    def ov(x0, x1, y0, y1):
        return max(0, min(x1, y1) - max(x0, y0))
    hs = [(qs - 1, qe, ss - 1, se) for qs, qe, ss, se in hsps if qe > qs and se > ss]
    seeds = [(min(ov(h[0], h[1], a0, a1), ov(h[2], h[3], b0, b1)), h) for h in hs]
    seeds = [x for x in seeds if x[0] > 0]
    if not seeds:
        return None
    pieces = [max(seeds)[1]]
    while True:
        c = pieces[0]
        left = [h for h in hs if h[0] < c[0] and h[2] < c[2] and -gap <= c[0] - h[1] <= gap and -gap <= c[2] - h[3] <= gap]
        if not left:
            break
        pieces.insert(0, min(left, key=lambda h: (h[0], h[2])))
    while True:
        c = pieces[-1]
        right = [h for h in hs if h[1] > c[1] and h[3] > c[3] and -gap <= h[0] - c[1] <= gap and -gap <= h[2] - c[3] <= gap]
        if not right:
            break
        pieces.append(max(right, key=lambda h: (h[1], h[3])))
    return pieces


def end_cover(pieces, pos_a, pos_b, upstream):
    """How a chain covers one core end (A position pos_a, B position pos_b; the first core base for upstream, one past
    the last for downstream): 'covered' (inside one aligned piece in both copies), 'covered_at_chain_end' (the
    outermost piece stops exactly at the end: a measured 0), 'end_in_chain_gap' (the end is not inside an aligned piece,
    or sits where an inner piece starts / stops, so it borders an unaligned gap) or 'not_reached'."""
    q0, s0 = (0, 2) if upstream else (1, 3)
    if (pieces[0][0] > pos_a or pieces[0][2] > pos_b) if upstream else (pieces[-1][1] < pos_a or pieces[-1][3] < pos_b):
        return 'not_reached'
    for k, h in enumerate(pieces):
        if h[0] < pos_a < h[1] and h[2] < pos_b < h[3]:
            return 'covered'
        if h[q0] == pos_a and h[s0] == pos_b and k == (0 if upstream else len(pieces) - 1):
            return 'covered_at_chain_end'
    return 'end_in_chain_gap'


def extensions(pieces, a0, a1, alen, b0, b1, blen, edge=EDGE):
    """Shared bp beyond the core on each side for A and B, how each core end is covered (end_cover), whether a side's
    shared stretch reaches a sequence end, and how many unaligned chain gaps lie beyond the end. Extensions are None
    (unknown) unless the end is covered in both copies; an end in a chain gap or not reached is never a 0."""
    qs, qe, ss, se = pieces[0][0], pieces[-1][1], pieces[0][2], pieces[-1][3]
    cu = end_cover(pieces, a0, b0, True)
    cd = end_cover(pieces, a1, b1, False)
    up, down = cu.startswith('covered'), cd.startswith('covered')
    return dict(up_cover=cu, down_cover=cd,
                a_up=a0 - qs if up else None, b_up=b0 - ss if up else None,
                a_down=qe - a1 if down else None, b_down=se - b1 if down else None,
                up_gaps_beyond_end=sum(h[1] <= a0 for h in pieces[:-1]) if up else '',
                down_gaps_beyond_end=sum(h[0] >= a1 for h in pieces[1:]) if down else '',
                up_reaches_end=up and (qs <= edge or ss <= edge), down_reaches_end=down and (qe >= alen - edge or se >= blen - edge))


def side_boundary(values, reached, spread=SPREAD):
    """Per occurrence side from one extension per different-context partner (None = end not covered):
    (median of known extensions or None, status, partners with unknown extension)."""
    known = [v for v in values if v is not None]
    unknown = len(values) - len(known)
    if not known:
        return None, 'end_not_covered', unknown
    med = int(statistics.median(known))
    if reached:
        return med, 'reaches_window_or_contig_end', unknown
    if unknown:
        return med, 'end_not_covered_by_some_partners', unknown
    if max(known) - min(known) > spread:
        return med, 'inconsistent_between_partners', unknown
    return med, 'resolved', unknown


def end_hints(seq, e0, e1, n=20, tsd_max=12):
    """Inverted-repeat identity of the outer n bp and longest exact direct repeat just outside [e0, e1) (hints only)."""
    left, right = seq[e0:e0 + n], revcomp(seq[e1 - n:e1])
    tir = sum(x == y for x, y in zip(left, right)) / n if len(left) == n == len(right) else None
    tsd = 0
    for k in range(tsd_max, 1, -1):
        if e0 - k >= 0 and e1 + k <= len(seq) and seq[e0 - k:e0].upper() == seq[e1:e1 + k].upper():
            tsd = k
            break
    return tir, tsd


def scan(genes, lo, hi, forward, labels, mobile, max_scan=MAX_SCAN, keep=None):
    """Neighbours per side in core orientation until keep (default K) non-mobile genes or max_scan genes are seen:
    lists of (label, same direction as core, gene_uid, mobile-like). Wraps only on circular contigs."""
    n, circ = len(genes), genes[0]['circular_contig'] == 'True'

    def side(step, start):
        out, j = [], start
        while sum(not x[3] for x in out) < (keep or pc.K) and len(out) < max_scan:
            j += step
            if not 0 <= j < n:
                if not circ:
                    break
                j %= n
            if lo <= j <= hi:
                break
            g = genes[j]
            out.append((labels[g['gene_uid']], (g['strand'] == '+') == forward, g['gene_uid'], mobile(g)))
        return out
    left, right = side(-1, lo), side(+1, hi)
    return (left, right) if forward else (right, left)


def host(scanned):
    return [x[:3] for x in scanned if not x[3]][:pc.K]


def blast(prog, records, task):
    """All-vs-all within one core: best HSP per (query, subject) for proteins, every HSP for DNA."""
    with tempfile.NamedTemporaryFile('w', suffix='.fa') as f:
        for k, s in records.items():
            f.write('>%s\n%s\n' % (k, s))
        f.flush()
        cmd = [prog, '-query', f.name, '-subject', f.name, '-outfmt', OUTFMT, '-evalue', '1e-5']
        if prog == 'blastn':
            cmd += ['-task', task, '-strand', 'plus', '-dust', 'no']
        out = subprocess.run(cmd, check=True, capture_output=True, text=True).stdout
    hits = defaultdict(list)
    for line in out.splitlines():
        q, s, pid, ln, qs, qe, ss, se, ql, sl, bits = line.split('\t')
        if q != s:
            hits[(q, s)].append(dict(pident=float(pid), length=int(ln), qs=int(qs), qe=int(qe), ss=int(ss), se=int(se),
                                     qlen=int(ql), slen=int(sl), bits=float(bits)))
    return hits


def best(hits):
    if not hits:
        return None
    h = max(hits, key=lambda h: h['bits'])
    return h['pident'], h['length'] / h['qlen'], h['length'] / h['slen']


def copy_pattern(copies):
    """Descriptive copy-number pattern over the subset assemblies (0 = absent); not a paralogy or mobility call."""
    present = {c for c in copies if c}
    if len(present) > 1:
        return 'count_varies'
    return 'same_count_in_all_assemblies' if all(copies) else 'same_count_where_present'


def load_ani(path, fasta):
    """{frozenset(sample pair): ANI} for the manifest samples from skani edges keyed by FASTA path."""
    smp = {v: k for k, v in fasta.items()}
    out = {}
    for r in pc.read_tsv(path):
        x, y = smp.get(r['Ref_file']), smp.get(r['Query_file'])
        if x and y and x != y:
            out[frozenset((x, y))] = float(r['ANI'])
    return out


def ani_bin(d):
    """Descriptive bins of core DNA identity minus genome ANI (percentage points); not an admission cutoff."""
    if d == '':
        return 'no_cross_pair_with_ani'
    return 'le_1' if d <= 1 else ('1_to_5' if d <= 5 else 'gt_5')


def core_genes(o):
    gs = o['genes'].split(';')
    return gs if o['strand'] == '+' else gs[::-1]


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--subset-dir', type=Path, required=True)
    p.add_argument('--work', type=Path, required=True, help='protein work dir (genes.tsv.gz, proteins.faa, mmseqs/)')
    p.add_argument('--manifest', type=Path, required=True)
    p.add_argument('--mmseqs', default='/global/scratch/users/kh36969/longread_callout/env/s4/bin/mmseqs')
    p.add_argument('--threads', type=int, default=8)
    p.add_argument('--skani', type=Path, default=Path('/global/scratch/users/kh36969/longread_callout/fna_pool/'
                   'scaleup_stage1_20261026/comamonadaceae/s2/skani_edges.tsv'))
    a = p.parse_args()
    pe = a.subset_dir / 'protein_entry'
    out = pe / 'core_review'
    out.mkdir(exist_ok=True)

    genes = pc.read_tsv(a.work / 'genes.tsv.gz')
    gene = {g['gene_uid']: g for g in genes}
    assign, _, _ = pc.assignments(Namespace(work=a.work, mmseqs=a.mmseqs, threads=a.threads))
    fam, hg = assign['family'], assign['homology_group']
    labels = {g['gene_uid']: pc.label(g, hg) for g in genes}
    by_contig = defaultdict(list)
    for g in genes:
        by_contig[(g['sample'], g['contig'])].append(g)
    for gs in by_contig.values():
        gs.sort(key=lambda g: int(g['index']))
    prot = {}
    for line in open(a.work / 'proteins.faa'):
        if line.startswith('>'):
            name = line[1:].strip()
        else:
            prot[name] = line.strip()

    cores = [c for c in pc.read_tsv(pe / 'protein_cores.tsv')
             if c['scope'] == 'target_lt2500' and int(c['different_context_candidate']) > 0]
    keys = {(c['structure'], c['core_key']) for c in cores}
    single_keys = {k for s, k in keys if s == 'single_gene'}
    occs = {o['occurrence_id']: o for o in pc.read_tsv(pe / 'protein_occurrences.tsv.gz')
            if o['in_target_length'] == 'True' and (o['structure'], o['core_key']) in keys}
    pairs = [r for r in pc.read_tsv(pe / 'protein_pairs.tsv.gz')
             if r['scope'] == 'target_lt2500' and (r['structure'], r['core_key']) in keys]

    def mobile(g):
        if hint(g['product']) in MOBILE_HINTS:
            return 'keyword'
        if g['type'] == 'CDS' and g['protein_sha256'] and fam.get(g['protein_sha256']) in single_keys:
            return 'reviewed_core_family'
        return ''

    for o in occs.values():
        idx = [int(u.rsplit('|', 1)[1]) for u in o['genes'].split(';')]
        gs = by_contig[(o['sample'], o['contig'])]
        su, sd = scan(gs, min(idx), max(idx), o['strand'] == '+', labels, mobile)
        o['scan'] = dict(up=su, down=sd)
        o['host'] = dict(up=host(su), down=host(sd))
        o['near'] = dict(up=[x[:3] for x in su[:pc.K]], down=[x[:3] for x in sd[:pc.K]])
        o['near_mobile'] = Counter(x[3] for s in ('up', 'down') for x in o['scan'][s][:pc.K] if x[3])

    manifest = pc.read_tsv(a.manifest)
    fasta = {r['sample']: r['fasta'] for r in manifest}
    samples = [r['sample'] for r in manifest]
    ani = load_ani(a.skani, fasta)
    seqs = {}
    for smp in sorted({o['sample'] for o in occs.values()}):
        seqs[smp] = pc.read_fna(fasta[smp])
    for o in occs.values():
        c = by_contig[(o['sample'], o['contig'])][0]
        s, e, fwd = int(o['start0']), int(o['end0']), o['strand'] == '+'
        o['ext'], o['c0'], o['c1'] = window(seqs[o['sample']][o['contig']], s, e, fwd, W, c['circular_contig'] == 'True')
        o['core_dna'] = o['ext'][o['c0']:o['c1']]

    pair_rows = []
    diff_by_core = defaultdict(list)
    for r in pairs:
        if r['relation'] == 'different_context_candidate':
            diff_by_core[(r['structure'], r['core_key'])].append(r)
    for key, ps in sorted(diff_by_core.items()):
        need = {x for r in ps for x in (r['occurrence_a'], r['occurrence_b'])}
        ids = {oid: str(i) for i, oid in enumerate(sorted(need))}
        npos = 1 if key[0] == 'single_gene' else 2
        php = [blast('blastp', {ids[oid] + '_%d' % i: prot[gene[core_genes(occs[oid])[i]]['protein_sha256']] for oid in need}, None)
               for i in range(npos)]
        dcore = blast('blastn', {ids[oid]: occs[oid]['core_dna'] for oid in need}, 'blastn')
        dext = blast('blastn', {ids[oid]: occs[oid]['ext'] for oid in need}, 'blastn')
        for r in ps:
            A, B = occs[r['occurrence_a']], occs[r['occurrence_b']]
            ia, ib = ids[A['occurrence_id']], ids[B['occurrence_id']]
            same_prot = all(gene[x]['protein_sha256'] == gene[y]['protein_sha256'] for x, y in zip(core_genes(A), core_genes(B)))
            pb = [best(php[i][(ia + '_%d' % i, ib + '_%d' % i)]) for i in range(npos)]
            db = best(dcore[(ia, ib)])
            blk = chain([(h['qs'], h['qe'], h['ss'], h['se']) for h in dext[(ia, ib)]], A['c0'], A['c1'], B['c0'], B['c1'])
            ext = extensions(blk, A['c0'], A['c1'], len(A['ext']), B['c0'], B['c1'], len(B['ext'])) if blk else \
                dict(up_cover='no_core_alignment', down_cover='no_core_alignment', a_up=None, b_up=None, a_down=None,
                     b_down=None, up_gaps_beyond_end='', down_gaps_beyond_end='', up_reaches_end=False, down_reaches_end=False)
            hrel = pc.relation(A['host'], B['host'])
            row = dict(structure=key[0], core_key=key[1], occurrence_a=A['occurrence_id'], occurrence_b=B['occurrence_id'],
                sample_a=A['sample'], sample_b=B['sample'], same_sample=r['same_sample'],
                identical_protein=same_prot,
                protein_identity_min=min(x[0] for x in pb) if all(pb) else '',
                protein_coverage_min=round(min(min(x[1], x[2]) for x in pb), 3) if all(pb) else '',
                identical_core_dna=A['core_dna'].upper() == B['core_dna'].upper(),
                core_dna_identity=db[0] if db else '', core_dna_coverage=round(min(db[1], db[2]), 3) if db else '',
                pair_ani=ani.get(frozenset((A['sample'], B['sample'])), ''),
                near_relation=r['relation'], host_relation=hrel,
                near_mobile_neighbours=sum(A['near_mobile'].values()) + sum(B['near_mobile'].values()), **ext)
            row['core_dna_minus_ani'] = round(row['core_dna_identity'] - row['pair_ani'], 2) \
                if row['pair_ani'] != '' and row['core_dna_identity'] != '' else ''
            pair_rows.append(row)
            for o, pre in ((A, 'a'), (B, 'b')):
                o.setdefault('partners', []).append((ext[pre + '_up'], ext['up_reaches_end'], ext[pre + '_down'],
                                                     ext['down_reaches_end'], ext['up_cover'], ext['down_cover']))

    occ_rows = []
    for oid, o in sorted(occs.items()):
        ps = o.get('partners', [])
        up = side_boundary([x[0] for x in ps], any(x[1] for x in ps))
        down = side_boundary([x[2] for x in ps], any(x[3] for x in ps))
        span = int(o['span'])
        resolved = up[1] == down[1] == 'resolved'
        est = span + up[0] + down[0] if resolved else None
        tir, tsd = end_hints(o['ext'], o['c0'] - up[0], o['c1'] + down[0]) if resolved else (None, None)
        occ_rows.append(dict(occurrence_id=oid, structure=o['structure'], core_key=o['core_key'], sample=o['sample'],
            contig=o['contig'], start0=o['start0'], end0=o['end0'], strand=o['strand'], coding_span=span,
            different_partners=len(ps), up_partners_end_not_covered=up[2], down_partners_end_not_covered=down[2],
            up_partners_end_in_chain_gap=sum(x[4] == 'end_in_chain_gap' for x in ps),
            down_partners_end_in_chain_gap=sum(x[5] == 'end_in_chain_gap' for x in ps),
            ext_up=up[0] if up[0] is not None else '', up_status=up[1],
            ext_down=down[0] if down[0] is not None else '', down_status=down[1],
            estimated_span=est if est is not None else '', boundary_resolved=resolved,
            estimated_lt2500='' if est is None else est < pc.TARGET_MAX,
            tir_identity_20bp='' if tir is None else round(tir, 2), direct_repeat_outside_bp='' if tsd is None else tsd,
            near_mobile_keyword=o['near_mobile']['keyword'], near_mobile_core_family=o['near_mobile']['reviewed_core_family'],
            near_up=','.join('%s%s' % (x[0] or 'nolabel', '' if x[1] else '(rev)') for x in o['near']['up']),
            near_down=','.join('%s%s' % (x[0] or 'nolabel', '' if x[1] else '(rev)') for x in o['near']['down']),
            host_up=','.join(x[2] for x in o['host']['up']), host_down=','.join(x[2] for x in o['host']['down']),
            host_up_products=' / '.join(gene[x[2]]['product'] for x in o['host']['up']),
            host_down_products=' / '.join(gene[x[2]]['product'] for x in o['host']['down'])))
    orow = {r['occurrence_id']: r for r in occ_rows}

    core_rows = []
    for c in sorted(cores, key=lambda c: (c['structure'], c['core_key'])):
        key = (c['structure'], c['core_key'])
        ps = [r for r in pair_rows if (r['structure'], r['core_key']) == key]
        os_ = [r for r in occ_rows if (r['structure'], r['core_key']) == key and r['different_partners']]
        within = [r for r in ps if r['same_sample'] == 'True']
        cross = [r for r in ps if r['same_sample'] != 'True']
        hd = lambda xs: sum(r['host_relation'] == 'different_context_candidate' for r in xs)
        nested = sum(r['host_relation'] == 'similar_context' for r in ps)
        if hd(cross):
            status = 'cross_assembly_host_different'
        elif hd(within):
            status = 'within_assembly_host_different_only'
        elif nested:
            status = 'nested_or_adjacent_mobile_only'
        else:
            status = 'host_context_unresolved'
        pid = [r['protein_identity_min'] for r in ps if r['protein_identity_min'] != '']
        did = [r['core_dna_identity'] for r in ps if r['core_dna_identity'] != '']
        est = [r['estimated_span'] for r in os_ if r['boundary_resolved']]
        dani = [r['core_dna_minus_ani'] for r in cross if r['core_dna_minus_ani'] != '']
        per = Counter(o['sample'] for o in occs.values() if (o['structure'], o['core_key']) == key)
        copies = [per[x] for x in samples]
        core_rows.append(dict(structure=key[0], core_key=key[1], products=c['products'], occurrences=c['occurrences'],
            copies_per_assembly=','.join(map(str, copies)), copy_pattern=copy_pattern(copies),
            samples=c['samples'], occurrences_in_different_pairs=len(os_),
            samples_in_different_pairs=len({r['sample'] for r in os_}),
            different_pairs_within_assembly=len(within), different_pairs_cross_assembly=len(cross),
            host_different_within=hd(within), host_different_cross=hd(cross), host_similar_nested=nested,
            pairs_identical_protein=sum(r['identical_protein'] for r in ps),
            protein_identity_range='%.1f..%.1f' % (min(pid), max(pid)) if pid else '',
            pairs_identical_core_dna=sum(r['identical_core_dna'] for r in ps),
            core_dna_identity_range='%.1f..%.1f' % (min(did), max(did)) if did else '',
            cross_median_core_dna_minus_ani=round(statistics.median(dani), 2) if dani else '',
            coding_span_range=c['span_range'], occurrences_boundary_resolved=len(est),
            estimated_span_range='%d..%d' % (min(est), max(est)) if est else '',
            resolved_estimated_lt2500=sum(x < pc.TARGET_MAX for x in est),
            occurrences_near_mobile=sum(bool(r['near_mobile_keyword'] or r['near_mobile_core_family']) for r in os_),
            status=status, product_hint=hint(c['products'])))

    pc.write_tsv(out / 'review_pairs.tsv.gz', pair_rows, list(pair_rows[0]) if pair_rows else [])
    occ_cols = list(occ_rows[0]) if occ_rows else []
    pc.write_tsv(out / 'review_occurrences.tsv.gz', [r for r in occ_rows if r['different_partners']], occ_cols)
    pc.write_tsv(out / 'review_cores.tsv', core_rows, list(core_rows[0]) if core_rows else [])

    summary = dict(version=VERSION, subset=a.subset_dir.name, evidence_source='assembly_only', scope='target_lt2500',
        parameters=dict(window_bp=W, chain_gap_bp=CHAIN_GAP, end_margin_bp=EDGE, spread_bp=SPREAD, max_scan=MAX_SCAN,
                        K=pc.K, mobile_hints=sorted(MOBILE_HINTS), blastn_task='blastn', evalue=1e-5),
        blast_version=subprocess.run(['blastn', '-version'], capture_output=True, text=True).stdout.split('\n')[0],
        site_groups_computed=False, production_ready=False, structures={})
    for st in ('single_gene', 'two_gene'):
        cs = [c for c in core_rows if c['structure'] == st]
        ps = [r for r in pair_rows if r['structure'] == st]
        os_ = [r for r in occ_rows if r['structure'] == st and r['different_partners']]
        res = [r for r in os_ if r['boundary_resolved']]
        summary['structures'][st] = dict(
            cores=len(cs), cores_by_status=dict(Counter(c['status'] for c in cs)),
            cores_with_cross_assembly_different_pair=sum(c['different_pairs_cross_assembly'] > 0 for c in cs),
            cores_with_within_assembly_different_pair=sum(c['different_pairs_within_assembly'] > 0 for c in cs),
            cores_within_assembly_only=sum(c['different_pairs_cross_assembly'] == 0 for c in cs),
            different_pairs=dict(within_assembly=sum(r['same_sample'] == 'True' for r in ps),
                                 cross_assembly=sum(r['same_sample'] != 'True' for r in ps)),
            host_relation_of_different_pairs={'%s:%s' % ('within' if r == 'True' else 'cross', h): n for (r, h), n in
                                              sorted(Counter((r['same_sample'], r['host_relation']) for r in ps).items())},
            pairs_identical_protein=sum(r['identical_protein'] for r in ps),
            pairs_protein_identity_lt95=sum(r['protein_identity_min'] != '' and r['protein_identity_min'] < 95 for r in ps),
            pairs_identical_core_dna=sum(r['identical_core_dna'] for r in ps),
            pairs_without_core_dna_hit=sum(r['core_dna_identity'] == '' for r in ps),
            cross_pairs_without_ani=sum(r['same_sample'] != 'True' and r['pair_ani'] == '' for r in ps),
            cores_by_cross_median_core_dna_minus_ani=dict(Counter(ani_bin(c['cross_median_core_dna_minus_ani']) for c in cs)),
            occurrences=len(os_), occurrences_boundary_resolved=len(res),
            occurrences_with_partner_end_in_chain_gap=sum(bool(r['up_partners_end_in_chain_gap'] or
                                                               r['down_partners_end_in_chain_gap']) for r in os_),
            pair_end_cover={side: dict(Counter(r[side + '_cover'] for r in ps)) for side in ('up', 'down')},
            cores_by_copy_pattern=dict(Counter(c['copy_pattern'] for c in cs)),
            occurrences_with_partner_not_covering_an_end=sum(bool(r['up_partners_end_not_covered'] or
                                                                   r['down_partners_end_not_covered']) for r in os_),
            occurrences_side_status={s: dict(Counter(r[s + '_status'] for r in os_)) for s in ('up', 'down')},
            resolved_estimated_lt2500=sum(r['estimated_lt2500'] is True for r in res),
            resolved_estimated_ge2500=sum(r['estimated_lt2500'] is False for r in res),
            resolved_extension_total_range='%d..%d' % (min(r['ext_up'] + r['ext_down'] for r in res),
                                                       max(r['ext_up'] + r['ext_down'] for r in res)) if res else '',
            occurrences_near_mobile_neighbour=sum(bool(r['near_mobile_keyword'] or r['near_mobile_core_family']) for r in os_))
    (out / 'review_summary.json').write_text(json.dumps(summary, indent=2) + '\n')

    cases = pick_cases(core_rows)
    rows_svg = []
    for tag, c in cases:
        key = (c['structure'], c['core_key'])
        os_ = sorted((occs[r['occurrence_id']] for r in occ_rows if (r['structure'], r['core_key']) == key and r['different_partners']),
                     key=lambda o: (o['context_class_target_lt2500'], o['sample'], o['occurrence_id']))
        xs, seen = [], set()
        for o in os_:
            if o['context_class_target_lt2500'] not in seen and len(xs) < 4:
                seen.add(o['context_class_target_lt2500']); xs.append(o)
        head = '%s %s: %s | %s occ / %s samples, copies %s; different pairs within %s, cross %s; host-different within %s, cross %s; ' \
               'protein id %s, core DNA id %s' % (tag, c['core_key'], c['products'][:50], c['occurrences'], c['samples'], c['copies_per_assembly'],
               c['different_pairs_within_assembly'], c['different_pairs_cross_assembly'], c['host_different_within'],
               c['host_different_cross'], c['protein_identity_range'], c['core_dna_identity_range'])
        rows_svg.append((head, [(o, orow[o['occurrence_id']]) for o in xs]))
    (out / 'cases_review.svg').write_text(draw(rows_svg, gene, '%s: reviewed U-core-D cases (target scope)' % a.subset_dir.name))
    case_rows = [dict(case=h.split(' ')[0], core_key=r['core_key'], occurrence_id=r['occurrence_id'], sample=r['sample'],
                      contig=r['contig'], start0=r['start0'], end0=r['end0'], strand=r['strand'], coding_span=r['coding_span'],
                      up_partners_end_not_covered=r['up_partners_end_not_covered'],
                      down_partners_end_not_covered=r['down_partners_end_not_covered'],
                      ext_up=r['ext_up'], up_status=r['up_status'], ext_down=r['ext_down'], down_status=r['down_status'],
                      estimated_span=r['estimated_span'], tir_identity_20bp=r['tir_identity_20bp'],
                      direct_repeat_outside_bp=r['direct_repeat_outside_bp'],
                      near_mobile=r['near_mobile_keyword'] + r['near_mobile_core_family'],
                      host_up_products_near_to_far=r['host_up_products'], host_down_products_near_to_far=r['host_down_products'])
                 for h, xs in rows_svg for _, r in xs]
    pc.write_tsv(out / 'case_occurrences.tsv', case_rows, list(case_rows[0]) if case_rows else ['case'])
    summary['cases'] = [h for h, _ in rows_svg]
    summary['case_rule'] = pick_cases.__doc__
    (out / 'review_summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))


def pick_cases(core_rows):
    """Fixed rule: A1/A2 = single-gene cores with the most cross-assembly host-different pairs; B = single-gene core
    with the most within-assembly host-different pairs and none across; C = core with the most near-different /
    host-similar (nested) pairs; D = two-gene core with the most host-different pairs. Ties: samples, then key; a core
    already picked is skipped for the next rank."""
    def top(xs, f):
        xs = [c for c in xs if f(c) > 0]
        return sorted(xs, key=lambda c: (-f(c), -int(c['samples']), c['core_key']))
    single = [c for c in core_rows if c['structure'] == 'single_gene']
    out = []

    def add(tag, xs, n=1):
        xs = [c for c in xs if all(c is not o[1] for o in out)][:n]
        out.extend((tag if n == 1 else '%s%d' % (tag, i + 1), c) for i, c in enumerate(xs))
    add('A', top(single, lambda c: c['host_different_cross']), 2)
    add('B', top([c for c in single if not c['host_different_cross']], lambda c: c['host_different_within']))
    add('C', top(core_rows, lambda c: c['host_similar_nested']))
    add('D', top([c for c in core_rows if c['structure'] == 'two_gene'], lambda c: c['host_different_cross'] + c['host_different_within']))
    return out


def draw(cases, gene, title):
    W_, y, box, gap = 1500, 74, 74, 4
    parts = []
    for heading, rows in cases:
        parts.append('<text x="14" y="%d" font-size="11" font-weight="bold">%s</text>' % (y, esc(heading)))
        y += 10
        for o, r in rows:
            bnd = 'conserved-end candidate %s bp (coding %s, +%s/+%s, %s/%s)' % (r['estimated_span'] or '?', r['coding_span'],
                  r['ext_up'] if r['ext_up'] != '' else '?', r['ext_down'] if r['ext_down'] != '' else '?', r['up_status'][:12], r['down_status'][:12])
            parts.append('<text x="14" y="%d" font-size="10">%s %s:%s-%s (%s)</text>' % (y + 14, o['sample'], o['contig'],
                         int(o['start0']) + 1, o['end0'], o['strand']))
            parts.append('<text x="14" y="%d" font-size="9" fill="#555">%s</text>' % (y + 25, esc(bnd)))
            up, down = o['scan']['up'][:6][::-1], o['scan']['down'][:6]
            x = 330 + (6 - len(up)) * (box + gap)
            for item in up + ['CORE'] + down:
                core = item == 'CORE'
                lab, same, uid, mob = ('core', True, '', '') if core else item
                fill = '#edc06f' if core else colour(lab or 'nolabel')
                pts = ('%d,%d %d,%d %d,%d %d,%d %d,%d' % (x, y + 6, x + box - 9, y + 6, x + box, y + 15, x + box - 9, y + 24, x, y + 24)
                       if same else '%d,%d %d,%d %d,%d %d,%d %d,%d' % (x + 9, y + 6, x + box, y + 6, x + box, y + 24, x + 9, y + 24, x, y + 15))
                stroke = 'stroke="#c0392b" stroke-width="2" stroke-dasharray="3,2"' if mob else 'stroke="#666" stroke-width="0.5"'
                tip = 'core' if core else '%s | %s%s' % (lab or 'nolabel', gene[uid]['product'], ' | mobile-like: ' + mob if mob else '')
                parts.append('<polygon points="%s" fill="%s" %s><title>%s</title></polygon>' % (pts, fill, stroke, esc(tip)))
                txt = 'core' if core else gene[uid]['product'][:13]
                parts.append('<text x="%d" y="%d" font-size="7">%s</text>' % (x + 10, y + 18, esc(txt)))
                x += box + gap
            y += 32
        y += 14
    head = ['<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d">' % (W_, y + 10),
            '<rect width="100%" height="100%" fill="white"/>',
            '<text x="14" y="24" font-size="16">%s</text>' % esc(title),
            '<text x="14" y="42" font-size="10">Assembly only. Up to 6 scanned neighbour genes per side in core orientation, '
            'coloured by homology group; red dashed outline = mobile-like neighbour skipped for host context. Conserved-end candidate '
            '= coding span + flank shared with all different-context copies (not an insertion boundary).</text>',
            '<text x="14" y="56" font-size="10">A = cross-assembly host-different, B = within-assembly only, C = nested / '
            'adjacent mobile under the current rule (near different, host similar), D = two-gene core. Hover a gene for its product.</text>']
    return '\n'.join(head + parts + ['</svg>']) + '\n'


if __name__ == '__main__':
    main()
