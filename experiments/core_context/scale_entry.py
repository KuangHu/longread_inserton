"""Assembly-only scale-up of the protein entry over the verified Comamonadaceae pool, in batches with one protein-family
mapping. Evidence stays assembly_only; production_ready=false; nothing here is an insertion call or a site group.

manifest     eligible pool assemblies (verified FNA + GFF3 + FAA, not atypical, every assembly level) from the frozen
             S2 catalog; file sha256 from each verify.json. Batch 1 = Complete Genome; batch 2 = Chromosome + Scaffold;
             batch 3 = Contig. Reads play no part in selection.
genes        per assembly (protein_context.genes_one: GFF3 units, S4 coding check, input hashes); resumable per assembly.
pool_proteins  every distinct protein (sha256) of every manifest assembly's FAA, read-only (no discovery).
families     one mapping for all batches. --pool (used here): all pool proteins clustered once at the frozen S4
             settings with the member check and repair (stage1_catalog.cmd_cluster + protein_context.assignments),
             before any batch is compared, so a family never depends on batch order. The incremental mode (--batch:
             batch 1 clustered, later batches searched against existing representatives with the same member
             thresholds, the rest clustered) is kept for additions outside the pool but not used: on the smoke run
             (5 subset assemblies split 3 + 2) it changed 10 of 242 occurrences with a different-context partner
             relative to one-pass clustering.
occurrences  single- and two-gene cores with K=3 contexts per assembly (protein_context.build_occurrences, unchanged),
             written to core-hash shards; occurrence ids are batch-independent (structure|sample|contig|gene index).
compare      per shard, over every batch given: exact, resumable pair accounting (scale_compare), for both scopes
             (all, target_lt2500), single- and two-gene, same-assembly and cross-assembly pairs.
summarize    completeness (pairs total / accounted / not compared per core), candidate counts, the regression set
             (frozen ucd_candidates_v0.1 located by coordinates), resource use.
"""
import argparse
import gzip
import hashlib
import json
import os
import subprocess
import sys
import time
from argparse import Namespace
from collections import Counter, defaultdict
from pathlib import Path

import protein_context as pc
import scale_compare as sc

sys.path.insert(0, str(pc.HERE.parent / 'fna_pool'))
from s4_annotate import MEMBER_MIN  # noqa: E402
from stage1_catalog import cmd_cluster  # noqa: E402

VERSION = 'scale_entry_v0.1'
SHARDS = 64
LEVEL_BATCH = {'Complete Genome': 1, 'Chromosome': 2, 'Scaffold': 2, 'Contig': 3}
SCOPES = (('all', lambda o: True), ('target_lt2500', lambda o: o['in_target_length'] == 'True'))
OCC_COLS = ['occurrence_id', 'structure', 'core_key', 'sample', 'contig', 'genes', 'start0', 'end0', 'span',
            'in_target_length', 'strand', 'spacing', 'homology_groups', 'products', 'up', 'down']
SEARCH = ['--alignment-mode', '3', '-e', '10', '-s', '7.5', '--max-seqs', '1000',
          '--format-output', 'query,target,fident,qcov,tcov,alnlen']


def shard_of(structure, key, n=SHARDS):
    return int(hashlib.sha1(('%s|%s' % (structure, key)).encode()).hexdigest()[:8], 16) % n


def batch_samples(a):
    return [m for m in pc.read_tsv(a.manifest) if int(m['batch']) == a.batch]


def cmd_manifest(a):
    rows, skipped = [], Counter()
    for r in pc.read_tsv(a.catalog):
        d = a.genomes / r['assembly']
        if r['atypical']:
            skipped['atypical'] += 1; continue
        if not (d / 'genomic.gff').exists() or not (d / 'protein.faa').exists():
            skipped['no_gff3_faa'] += 1; continue
        v = json.loads((d / 'verify.json').read_text())
        if v['status'] != 'verified' or not all(v['checks'].values()):
            skipped['verify_failed'] += 1; continue
        sha = {f['file']: f['sha256'] for f in v['files']}
        fna = [f for f in sha if f.endswith('_genomic.fna')]
        rows.append(dict(batch=LEVEL_BATCH[r['assembly_level']], sample=r['sample_id'], assembly=r['assembly'],
            assembly_level=r['assembly_level'], genus=r['genus'], organism=r['organism'], species_cluster=r['species_cluster'],
            clone_group=r['clone_group'], fasta=str(d / fna[0]), fna_sha256=sha[fna[0]], gff3=str(d / 'genomic.gff'),
            gff3_sha256=sha['genomic.gff'], faa=str(d / 'protein.faa'), faa_sha256=sha['protein.faa'],
            total_bp=r['total_bp'], records=r['records']))
    rows.sort(key=lambda m: (int(m['batch']), m['assembly']))
    pc.write_tsv(a.out, rows, list(rows[0]))
    print(json.dumps(dict(eligible=len(rows), skipped=skipped, batches=Counter(m['batch'] for m in rows),
                          levels=Counter(m['assembly_level'] for m in rows)), indent=1))


def cmd_genes(a):
    out = a.work / 'genes'
    out.mkdir(parents=True, exist_ok=True)
    ms = batch_samples(a)[a.task::a.ntasks]
    for m in ms:
        done = out / (m['sample'] + '.json')
        if done.exists():
            continue
        t = time.time()
        rows, prot, inputs = pc.genes_one(m)
        pc.write_tsv(out / (m['sample'] + '.tsv.gz'), rows, pc.GENE_COLS + ['circular_contig', 'contig_length'])
        with gzip.open(out / (m['sample'] + '.faa.gz'), 'wt') as f:
            for h in sorted(prot):
                f.write('>%s\n%s\n' % (h, prot[h]))
        done.write_text(json.dumps(dict(sample=m['sample'], assembly=m['assembly'], batch=a.batch, inputs=inputs,
                                        genes=len(rows), proteins=len(prot), seconds=round(time.time() - t, 1),
                                        coding_status=Counter(r['coding_status'] for r in rows))) + '\n')


def read_faa_gz(path):
    out, name = {}, None
    with gzip.open(path, 'rt') as f:
        for line in f:
            if line.startswith('>'):
                name = line[1:].strip()
            else:
                out[name] = line.strip()
    return out


def load_catalog(work):
    p = work / 'families' / 'catalog.tsv.gz'
    return {r['protein_sha256']: r for r in pc.read_tsv(p)} if p.exists() else {}


def cluster_new(work_dir, proteins, a):
    """Cluster `proteins` among themselves at the frozen settings; {level: {sha: (group id, representative)}}."""
    work_dir.mkdir(parents=True, exist_ok=False)
    with open(work_dir / 'proteins.faa', 'w') as f:
        for h in sorted(proteins):
            f.write('>%s\n%s\n' % (h, proteins[h]))
    ns = Namespace(work=work_dir, mmseqs=a.mmseqs, threads=a.threads)
    cmd_cluster(ns)
    assign, fixes, checks = pc.assignments(ns)
    pc.write_tsv(work_dir / 'member_check_repairs.tsv', fixes, list(fixes[0]) if fixes else ['level', 'member'])
    reps = {}
    for level, m in assign.items():
        rep_of = {}
        for line in open(work_dir / 'mmseqs' / (level + '_cluster.tsv')):
            r, x = line.split()
            rep_of[x] = r
        for f in fixes:
            if f['level'] == level:
                rep_of[f['member']] = f['new_rep']
        reps[level] = {h: (m[h], rep_of[h]) for h in m}
    return reps, checks, len(fixes)


def search_existing(work_dir, new, reps_faa, level, a):
    """Best existing representative passing the level's member thresholds for each new protein."""
    mid, mc = MEMBER_MIN[level]
    q, t, hits = work_dir / (level + '_new.faa'), work_dir / (level + '_reps.faa'), work_dir / (level + '_hits.m8')
    for path, recs in ((q, new), (t, reps_faa)):
        with open(path, 'w') as f:
            for h in sorted(recs):
                f.write('>%s\n%s\n' % (h, recs[h]))
    subprocess.run([a.mmseqs, 'easy-search', str(q), str(t), str(hits), str(work_dir / ('tmp_' + level)),
                    '--threads', str(a.threads)] + SEARCH, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    best = {}
    for line in open(hits):
        x, r, fid, qc, tc, _ = line.rstrip('\n').split('\t')
        fid, qc, tc = float(fid), float(qc), float(tc)
        if fid >= mid and qc >= mc and tc >= mc:
            key = (-fid, -min(qc, tc), r)
            if x not in best or key < best[x][0]:
                best[x] = (key, r, fid, qc, tc)
    return best


def cmd_pool_proteins(a):
    seen, n, per = {}, 0, []
    for m in pc.read_tsv(a.manifest):
        k = 0
        for name, seq in pc.read_faa(m['faa']).items():
            h = hashlib.sha256(seq.encode()).hexdigest()
            seen.setdefault(h, seq); k += 1
        n += k
        per.append(k)
    a.work.mkdir(parents=True, exist_ok=True)
    with gzip.open(a.work / 'pool_proteins.faa.gz', 'wt') as f:
        for h in sorted(seen):
            f.write('>%s\n%s\n' % (h, seen[h]))
    (a.work / 'pool_proteins.json').write_text(json.dumps(dict(assemblies=len(per), faa_records=n, distinct=len(seen),
                                                               manifest_sha256=pc.digest(a.manifest)), indent=1) + '\n')


def cmd_families(a):
    fdir = a.work / 'families'
    fdir.mkdir(parents=True, exist_ok=True)
    cat = load_catalog(a.work)
    if a.pool:
        if cat:
            raise ValueError('a family catalog exists already')
        proteins = read_faa_gz(a.work / 'pool_proteins.faa.gz')
        a.batch = 0
    else:
        proteins = {}
        for m in batch_samples(a):
            proteins.update(read_faa_gz(a.work / 'genes' / (m['sample'] + '.faa.gz')))
    new = {h: s for h, s in proteins.items() if h not in cat}
    info = dict(batch=a.batch, proteins_in_batch=len(proteins), already_mapped=len(proteins) - len(new), new=len(new))
    bdir = fdir / ('pool' if a.pool else 'b%d' % a.batch)
    if not cat:
        reps, checks, nfix = cluster_new(bdir, new, a)
        rows = [dict(protein_sha256=h, family=reps['family'][h][0], family_rep=reps['family'][h][1],
                     homology_group=reps['homology_group'][h][0], homology_group_rep=reps['homology_group'][h][1],
                     batch=a.batch, assigned='cluster', family_identity='', homology_group_identity='') for h in sorted(new)]
        info.update(member_check=checks, member_repairs=nfix)
    else:
        bdir.mkdir(parents=True, exist_ok=False)
        seq_of = dict(proteins)
        old_reps = {lv: {r[lv + '_rep'] for r in cat.values()} for lv in ('family', 'homology_group')}
        need = set().union(*old_reps.values()) - set(seq_of)
        for b in range(1, a.batch):
            for m in [m for m in pc.read_tsv(a.manifest) if int(m['batch']) == b]:
                if not need:
                    break
                got = read_faa_gz(a.work / 'genes' / (m['sample'] + '.faa.gz'))
                for h in need & set(got):
                    seq_of[h] = got[h]
                need -= set(got)
        if need:
            raise ValueError('representative sequences missing: %d' % len(need))
        rows, per = [], {}
        for lv in ('family', 'homology_group'):
            per[lv] = search_existing(bdir, new, {h: seq_of[h] for h in old_reps[lv]}, lv, a)
        left = {h: s for h, s in new.items() if h not in per['family'] or h not in per['homology_group']}
        reps = cluster_new(bdir / 'new', left, a)[0] if left else {'family': {}, 'homology_group': {}}
        id_of = {lv: {r[lv + '_rep']: r[lv] for r in cat.values()} for lv in ('family', 'homology_group')}
        used = {lv: set(id_of[lv].values()) for lv in id_of}
        for h in sorted(new):
            row = dict(protein_sha256=h, batch=a.batch)
            how = []
            for lv in ('family', 'homology_group'):
                if h in per[lv]:
                    _, r, fid, qc, tc = per[lv][h]
                    row.update({lv: id_of[lv][r], lv + '_rep': r, lv + '_identity': fid}); how.append(lv + ':search')
                else:
                    gid, r = reps[lv][h]
                    if gid in used[lv] and r not in id_of[lv]:
                        raise ValueError('new %s id collides with an existing one: %s' % (lv, gid))
                    row.update({lv: gid, lv + '_rep': r, lv + '_identity': ''}); how.append(lv + ':cluster')
            row['assigned'] = ','.join(how)
            rows.append(row)
        info.update(assigned=Counter(r['assigned'] for r in rows))
    cols = ['protein_sha256', 'family', 'family_rep', 'homology_group', 'homology_group_rep', 'batch', 'assigned',
            'family_identity', 'homology_group_identity']
    allrows = list(cat.values()) + rows
    pc.write_tsv(fdir / 'catalog.tmp.tsv.gz', allrows, cols)
    os.replace(fdir / 'catalog.tmp.tsv.gz', fdir / 'catalog.tsv.gz')
    info.update(families=len({r['family'] for r in allrows}), homology_groups=len({r['homology_group'] for r in allrows}))
    (fdir / ('pool.json' if a.pool else 'b%d.json' % a.batch)).write_text(json.dumps(info, indent=1) + '\n')
    print(json.dumps(info, indent=1))


def cmd_occurrences(a):
    cat = load_catalog(a.work)
    fam = {h: r['family'] for h, r in cat.items()}
    hg = {h: r['homology_group'] for h, r in cat.items()}
    out = a.work / 'occ' / ('b%d' % a.batch)
    out.mkdir(parents=True, exist_ok=False)
    fh = [gzip.open(out / ('shard_%03d.tsv.gz' % k), 'wt') for k in range(SHARDS)]
    for f in fh:
        f.write('\t'.join(OCC_COLS) + '\n')
    stat = Counter()
    for m in batch_samples(a):
        genes = pc.read_tsv(a.work / 'genes' / (m['sample'] + '.tsv.gz'))
        labels = {g['gene_uid']: pc.label(g, hg) for g in genes}
        by_contig = defaultdict(list)
        for g in genes:
            by_contig[(g['sample'], g['contig'])].append(g)
        for o in pc.build_occurrences(by_contig, fam, hg, labels):
            first = o['genes'].split(';')[0].rsplit('|', 1)[1]
            o['occurrence_id'] = '%s|%s|%s|%s' % (o['structure'][0], o['sample'], o['contig'], first)
            o['up'], o['down'] = json.dumps(o['up']), json.dumps(o['down'])
            fh[shard_of(o['structure'], o['core_key'])].write('\t'.join(str(o[c]) for c in OCC_COLS) + '\n')
            stat[o['structure']] += 1
    for f in fh:
        f.close()
    (out / 'occurrences.json').write_text(json.dumps(dict(batch=a.batch, assemblies=len(batch_samples(a)),
                                                          occurrences=stat), indent=1) + '\n')


def read_shard(work, batches, k):
    occs = []
    for b in batches:
        occs += pc.read_tsv(work / 'occ' / ('b%d' % b) / ('shard_%03d.tsv.gz' % k))
    return occs


def labels_of(side):
    return [x[0] for x in json.loads(side)]


def cmd_compare(a):
    batches = [int(x) for x in a.batches.split(',')]
    cdir = a.work / 'compare' / ('b' + '-'.join(map(str, batches))) / ('shard_%03d' % a.shard)
    cdir.mkdir(parents=True, exist_ok=True)
    prog = cdir / 'progress.jsonl'
    done = defaultdict(list)
    if prog.exists():
        for line in open(prog):
            try:
                r = json.loads(line)
            except ValueError:
                continue      # a record cut by a kill is not a finished block
            if r.get('v') == sc.RECORD_VERSION:      # older record formats are recomputed
                done[(r['scope'], r['structure'], r['core_key'], r['digest'])].append(r)
    occs = read_shard(a.work, batches, a.shard)
    groups = defaultdict(list)
    for o in occs:
        groups[(o['structure'], o['core_key'])].append(o)
    budget = a.max_sig_pairs
    used = 0
    core_rows, sig_rows = [], []
    with open(prog, 'a') as pf:
        for (structure, key), xs in sorted(groups.items()):
            for scope, keep in SCOPES:
                ys = sorted((o for o in xs if keep(o)), key=lambda o: o['occurrence_id'])
                if len(ys) < 2:
                    continue
                sigs, mult, index, dg = sc.collapse([(o['occurrence_id'], o['sample'], labels_of(o['up']), labels_of(o['down']))
                                                     for o in ys])
                recs = {(r['i0'], r['i1']): r for r in done[(scope, structure, key, dg)]}
                for i0, i1 in sc.blocks(mult):
                    if (i0, i1) in recs:
                        continue
                    work_pairs = sum(len(sigs) - i for i in range(i0, i1))
                    if budget is not None and used + work_pairs > budget and used > 0:
                        continue
                    rec = sc.run_block(sigs, mult, i0, i1)
                    rec.update(scope=scope, structure=structure, core_key=key, digest=dg)
                    pf.write(json.dumps(rec) + '\n'); pf.flush(); os.fsync(pf.fileno())
                    recs[(i0, i1)] = rec
                    used += work_pairs
                summary, srows = sc.finalize(sigs, mult, list(recs.values()))
                per = Counter(o['sample'] for o in ys)
                core_rows.append(dict(scope=scope, structure=structure, core_key=key, **summary,
                                      max_copies_in_one_assembly=max(per.values()),
                                      products=Counter(o['products'] for o in ys).most_common(1)[0][0][:200],
                                      span_range='%s..%s' % (min(int(o['span']) for o in ys), max(int(o['span']) for o in ys))))
                named = sc.name_partners([(o['occurrence_id'], o['sample']) for o in ys], index, srows)
                for o in ys:
                    s = srows[index[o['occurrence_id']]]
                    sig_rows.append(dict(scope=scope, structure=structure, core_key=key, occurrence_id=o['occurrence_id'],
                                         sample=o['sample'], sig=s['sig'], context_class=s['context_class'],
                                         **named[o['occurrence_id']]))
    pc.write_tsv(cdir / 'cores.tsv.gz', core_rows, list(core_rows[0]) if core_rows else ['scope'])
    pc.write_tsv(cdir / 'occurrence_sigs.tsv.gz', sig_rows, list(sig_rows[0]) if sig_rows else ['scope'])
    (cdir / 'run.json').write_text(json.dumps(dict(shard=a.shard, batches=batches, cores=len(core_rows),
        complete=sum(r['complete'] for r in core_rows), sig_pairs_this_run=used, budget=budget), indent=1) + '\n')


def cmd_summarize(a):
    batches = [int(x) for x in a.batches.split(',')]
    tag = 'b' + '-'.join(map(str, batches))
    cores = []
    missing = []
    for k in range(SHARDS):
        p = a.work / 'compare' / tag / ('shard_%03d' % k) / 'cores.tsv.gz'
        if p.exists():
            cores += pc.read_tsv(p)
        else:
            missing.append(k)
    out = dict(version=VERSION, batches=batches, evidence_source='assembly_only', production_ready=False,
               shards_missing=missing, completeness={}, candidates={})
    for scope in ('all', 'target_lt2500'):
        for st in ('single_gene', 'two_gene'):
            cs = [c for c in cores if c['scope'] == scope and c['structure'] == st]
            key = '%s:%s' % (scope, st)
            out['completeness'][key] = dict(cores=len(cs), complete=sum(c['complete'] == 'True' for c in cs),
                pairs_total=sum(int(c['pairs_total']) for c in cs), pairs_accounted=sum(int(c['pairs_accounted']) for c in cs),
                pairs_not_compared=sum(int(c['pairs_not_compared']) for c in cs),
                cores_above_old_cap=sum(int(c['pairs_total']) > pc.MAX_PAIRS for c in cs),
                largest_core_occurrences=max((int(c['occurrences']) for c in cs), default=0),
                largest_core_signatures=max((int(c['signatures']) for c in cs), default=0))
            d = [c for c in cs if int(c['different_context_candidate']) > 0]
            cross = lambda c: int(c['different_context_candidate']) > int(c['different_context_candidate_same_assembly'])
            out['candidates'][key] = dict(cores_2plus=len(cs), cores_with_different_pair=len(d),
                with_cross_assembly_different=sum(cross(c) for c in d),
                with_same_assembly_different=sum(int(c['different_context_candidate_same_assembly']) > 0 for c in d),
                cross_assembly_only=sum(cross(c) and int(c['different_context_candidate_same_assembly']) == 0 for c in d),
                single_copy_per_assembly=sum(c['max_copies_in_one_assembly'] == '1' for c in d),
                different_pairs=sum(int(c['different_context_candidate']) for c in d),
                different_pairs_same_assembly=sum(int(c['different_context_candidate_same_assembly']) for c in d))
    tdir = a.repo_out
    tdir.mkdir(parents=True, exist_ok=True)
    cand = [c for c in cores if c['scope'] == 'target_lt2500' and int(c['different_context_candidate']) > 0]
    pc.write_tsv(tdir / ('target_cores_with_different_pair_%s.tsv.gz' % tag), cand, list(cand[0]) if cand else ['scope'])
    if a.regression:
        out['regression'] = regression(a, batches, tag)
    (tdir / ('summary_%s.json' % tag)).write_text(json.dumps(out, indent=1) + '\n')
    print(json.dumps(out, indent=1))


def regression(a, batches, tag):
    """Locate every frozen candidate pair by coordinates and recheck it in the scaled run."""
    cands = pc.read_tsv(a.regression)

    def loc_key(c, loc):
        contig, rest = loc.rsplit(':', 1)
        s, e = rest.split('(')[0].split('-')
        return (c['sample'], contig, str(int(s) - 1), e)
    want = {loc_key(c, c[x]) for c in cands for x in ('loc_a', 'loc_b')}
    occ_loc, shards = {}, {}
    for k in range(SHARDS):
        for o in read_shard(a.work, batches, k):
            key = (o['sample'], o['contig'], o['start0'], o['end0'])
            if o['structure'] == 'single_gene' and key in want:
                occ_loc[key] = o
                shards[o['occurrence_id']] = k
    sig = {}
    rows = []
    for c in cands:
        res = dict(candidate=c['candidate'], batch=c['batch'], subset=c['subset'], old_core_key=c['core_key'], sample=c['sample'])
        found = [occ_loc.get(loc_key(c, c[x])) for x in ('loc_a', 'loc_b')]
        if not all(found):
            res.update(status='occurrence_not_found'); rows.append(res); continue
        A, B = found
        k = shards[A['occurrence_id']]
        if k not in sig:
            sig[k] = {(r['scope'], r['occurrence_id']): r for r in
                      pc.read_tsv(a.work / 'compare' / tag / ('shard_%03d' % k) / 'occurrence_sigs.tsv.gz')}
        rel = sc.relation_sig((sc.side_sig(labels_of(A['up'])), sc.side_sig(labels_of(A['down']))),
                              (sc.side_sig(labels_of(B['up'])), sc.side_sig(labels_of(B['down']))))
        sa, sb = sig[k].get(('target_lt2500', A['occurrence_id'])), sig[k].get(('target_lt2500', B['occurrence_id']))
        res.update(new_core_key_a=A['core_key'], new_core_key_b=B['core_key'], same_core=A['core_key'] == B['core_key'],
                   relation_k3=rel, occurrence_a=A['occurrence_id'], occurrence_b=B['occurrence_id'],
                   a_different_partners_cross_same='%s/%s' % (sa['diff_partners_cross'], sa['diff_partners_same_assembly']) if sa else '',
                   b_different_partners_cross_same='%s/%s' % (sb['diff_partners_cross'], sb['diff_partners_same_assembly']) if sb else '',
                   status='reproduced' if A['core_key'] == B['core_key'] and rel == 'different_context_candidate' else 'changed')
        rows.append(res)
    cols = ['candidate', 'batch', 'subset', 'old_core_key', 'sample', 'status', 'new_core_key_a', 'new_core_key_b',
            'same_core', 'relation_k3', 'occurrence_a', 'occurrence_b', 'a_different_partners_cross_same',
            'b_different_partners_cross_same']
    pc.write_tsv(a.repo_out / ('regression_%s.tsv' % tag), rows, cols)
    return dict(Counter(r['status'] for r in rows))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest='mode', required=True)
    s = sub.add_parser('manifest')
    s.add_argument('--catalog', type=Path, required=True)
    s.add_argument('--genomes', type=Path, required=True)
    s.add_argument('--out', type=Path, required=True)
    for name in ('genes', 'families', 'occurrences'):
        s = sub.add_parser(name)
        s.add_argument('--manifest', type=Path, required=True)
        s.add_argument('--work', type=Path, required=True)
        s.add_argument('--batch', type=int, required=True)
        if name == 'genes':
            s.add_argument('--task', type=int, default=0)
            s.add_argument('--ntasks', type=int, default=1)
        if name == 'families':
            s.add_argument('--pool', action='store_true', help='cluster every pool protein once (pool_proteins)')
            s.add_argument('--mmseqs', default='/global/scratch/users/kh36969/longread_callout/env/s4/bin/mmseqs')
            s.add_argument('--threads', type=int, default=8)
    s = sub.add_parser('pool_proteins')
    s.add_argument('--manifest', type=Path, required=True)
    s.add_argument('--work', type=Path, required=True)
    s = sub.add_parser('compare')
    s.add_argument('--work', type=Path, required=True)
    s.add_argument('--batches', required=True)
    s.add_argument('--shard', type=int, required=True)
    s.add_argument('--max-sig-pairs', type=int, default=None)
    s = sub.add_parser('summarize')
    s.add_argument('--work', type=Path, required=True)
    s.add_argument('--batches', required=True)
    s.add_argument('--repo-out', type=Path, required=True)
    s.add_argument('--regression', type=Path)
    a = p.parse_args()
    dict(manifest=cmd_manifest, pool_proteins=cmd_pool_proteins, genes=cmd_genes, families=cmd_families, occurrences=cmd_occurrences,
         compare=cmd_compare, summarize=cmd_summarize)[a.mode](a)


if __name__ == '__main__':
    main()
