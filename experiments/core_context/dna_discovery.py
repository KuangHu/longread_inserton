"""DNA-first discovery entry for U-X-D candidates (discovery_source=dna). No protein, gene, product or gap input.

    position A:  U_A - X  - D_A
    position B:  U_B - X' - D_B        X / X' similar DNA, U and D different (decided only by the two checks)

Input: verified assembly FASTA (FNA) and the source manifest (sample, assembly, fasta, fna_sha256, batch, family).
Every FASTA is checked against its manifest sha256 before it is read; a mismatch is recorded and the assembly is
not used. Nothing here reads a protein family, a CDS annotation, a product name, a gap list, S6 or read evidence.

Stages (each resumable; one SLURM array task per unit; outputs under --work)
  plan         split the assemblies of the given batches into seed chunks (same family together, then by
               species cluster, so that related genomes share a chunk) and genome parts for the occurrence search.
  seeds        per seed chunk: all-vs-all local DNA alignment of the chunk's assemblies (minimap2 -x asm20, self
               diagonal skipped, same-contig repeats kept). A record whose spans on both sides are SEED_MIN..SEED_MAX
               bp is a seed: a repeated / homologous DNA segment. Seeds are only seeds: their ends are where the
               aligner stopped, not boundaries. Near-identical seed loci are merged, sequences clustered
               (greedy containment, mmseqs easy-cluster) into chunk representatives.
  library      new chunk representatives that already match a library sequence are dropped; the rest are clustered
               into a new library shard (L<k>_<n>). Earlier shards are never renumbered.
  occurrences  per (library shard, genome part): every library sequence is searched in every assembly of the part
               (minimap2 -x asm20, all primary and secondary chains kept, no copy cap). A hit is an occurrence of the
               library sequence, oriented by the hit strand (one orientation per library sequence). The occurrence
               window (seed +/- W bp, forward strand) is stored with it. New shards are searched in old parts and old
               shards in new parts, so cores shared between batches are found.
  pairs        per library-sequence shard: occurrences of one library sequence are paired (same assembly and across
               assemblies). Pair-specific candidate X: the two occurrence windows in the library orientation are
               aligned (blastn, as core_review.blast) and the shared block around the seeds is chained and extended
               (core_review.chain / extensions, unchanged). An end is resolved only if it is covered in both copies
               and the shared block does not reach a window or contig end; then X = seed + shared extension on each
               side (never cut to a length). Otherwise the pair is unresolved and is not checked.
               High-copy library sequences: all pairs up to FULL_PAIRS per sequence; above that each occurrence gets
               DESIGN_PARTNERS partners (same-assembly copies first, then other assemblies by a fixed hash order) and
               the pairs not formed are counted per sequence (pairs_deferred), never dropped silently.
  collect      resolved pairs of all shards, de-duplicated by the candidate intervals (one pair of X intervals found
               through several library sequences is kept once, the others counted), written as an ucd_batch input
               table (candidate_insertion only). ucd_batch then runs the two checks unchanged.

Seeds and library sequences are discovery devices. A library hit is not an insertion, and the extended X is a
candidate interval for the checks, not an insertion call. Discovery scope: a segment is found if it is a seed in at
least one chunk (two copies inside one chunk, or one copy repeated in one assembly) or if it is similar to a library
sequence; a core whose copies never meet in one chunk and that matches no library sequence is not found.
"""
import argparse
import gzip
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
from collections import Counter, defaultdict
from pathlib import Path

import core_review as cr
import ucd_batch as ub
import ucd_compare as uc

VERSION = 'dna_discovery_v0.2.2'
# v0.2.2: collect sends only pairs with 700 <= X < 2500 bp on both sides to the checks (others: out_of_target)
# v0.2.1: alignment keys sharded by subject window (one blastn call per subject with all its queries); per-key
#         results unchanged (one subject per call in both)
# v0.2.0: enumerate / align replace pairs: canonical alignment keys (grid windows, relative orientation, A-B = B-A)
#         computed once and cached; one subject per blastn call (result independent of batch composition);
#         high-copy pairs enumerated completely in deterministic chunks (pair_chunk), deferred chunks resumable.
# v0.1.1: minimap2 -c (v0.1 filtered identity on approximate PAF match counts)
# v0.1.2: occurrence search with minimap2 -P (v0.1.1 used -N 1000, which caps secondary chains per query contig and
#         silently dropped copies of high-copy segments); seed / library sequences clustered by greedy containment
#         (mmseqs easy-cluster, cluster-mode 2, member coverage), so a segment contained in a longer one is not kept
#         as a redundant library sequence (linclust compared each sequence only with one k-mer centre)
ENV = Path('/global/scratch/users/kh36969/longread_callout/env')
MINIMAP2 = str(ENV / 'msi' / 'bin' / 'minimap2')
MMSEQS = str(ENV / 's4' / 'bin' / 'mmseqs')
PARAMS = dict(
    seed_preset='asm20',          # minimap2 preset of the all-vs-all seed search and of the occurrence search
    seed_min=100, seed_max=10000,  # bp of both spans of a seed record (longer records are synteny, not seeds)
    seed_identity=0.80,           # matches / alignment length (base-level, minimap2 -c) of a seed record and an occurrence hit
    # seed and occurrence searches use minimap2 -P: all chains kept, -N / -p have no effect (no copy cap)
    merge_overlap=0.8,            # seed loci on one contig overlapping >= this of both are one locus
    cluster_identity=0.90, cluster_coverage=0.80,  # greedy containment clustering of seed / library sequences:
                                                   # a sequence covered >= coverage at >= identity by a longer one
                                                   # (either strand) is a member, the longest is the representative
    occ_min=100,                  # bp of an occurrence hit
    window=cr.W,                  # bp each side of the seed for the pair extension (core_review.W)
    chain_gap=cr.CHAIN_GAP, edge=cr.EDGE,
    full_pairs=20000,             # v0.1.2 pairs stage only: all pairs up to this many, design above (legacy)
    design_partners=8,            # partners per occurrence of the cover design (v0.2: priority chunk 0 only)
    pair_chunk=20000,             # v0.2: pairs per enumeration chunk (a chunk size, not a cap)
    window_grid=1000,             # v0.2: alignment windows widened outward to this grid (>= W each side), so copies
                                  #       seen through several library sequences share one alignment
    key_shards=128,               # v0.2: alignment-key shards (one align task each)
    target_min=700, target_max=2500,  # v0.2.2: candidate X length of both sides before the two checks (collect);
                                  #       seeds keep 100..10,000 bp, X is never cut to fit
    buckets=64,                   # occurrence files are split by library sequence into this many buckets
)
KIND = 'candidate_insertion'


def digest_params():
    return hashlib.sha256(json.dumps(PARAMS, sort_keys=True).encode()).hexdigest()[:16]


def rows(path):
    return uc.read_tsv(path)


def write(path, rs, cols):
    path = Path(path)
    tmp = path.with_name('.tmp_' + path.name)  # keeps the .gz suffix, so compression matches the final name
    uc.write_tsv(tmp, rs, cols)
    os.replace(tmp, path)


def manifest(path, batches=None):
    ms = [m for m in rows(path) if batches is None or int(m['batch']) in batches]
    return ms


def verified_fasta(m):
    """Sequences of a manifest assembly after its sha256 check, or None."""
    if ub.sha_file(m['fasta']) != m['fna_sha256']:
        return None
    return uc.read_fasta(m['fasta'])


def run(cmd, **kw):
    return subprocess.run(cmd, check=True, **kw)


# ---------------------------------------------------------------- plan

def cmd_plan(a):
    w = Path(a.work)
    w.mkdir(parents=True, exist_ok=True)
    batches = [int(x) for x in a.batches.split(',')]
    ms = manifest(a.manifest, set(batches))
    done = set()
    if (w / 'chunks.tsv').exists():
        done = {s for r in rows(w / 'chunks.tsv') for s in r['samples'].split(';')}
    new = [m for m in ms if m['sample'] not in done]
    new.sort(key=lambda m: (m.get('family', ''), m['species_cluster'], m['sample']))
    old_chunks = rows(w / 'chunks.tsv') if (w / 'chunks.tsv').exists() else []
    old_parts = rows(w / 'parts.tsv') if (w / 'parts.tsv').exists() else []
    tag = a.tag
    chunks, parts = [], []
    fam = defaultdict(list)
    for m in new:
        fam[m.get('family', '')].append(m)
    for f in sorted(fam):  # a seed chunk never mixes families (smaller last chunk of a family is kept)
        ms_f = fam[f]
        for i in range(0, len(ms_f), a.chunk_size):
            part = ms_f[i:i + a.chunk_size]
            chunks.append(dict(chunk='%s_s%03d' % (tag, len(chunks)), tag=tag, family=f, n=len(part),
                               bp=sum(int(m['total_bp']) for m in part), samples=';'.join(m['sample'] for m in part)))
    for i in range(0, len(new), a.part_size):
        part = new[i:i + a.part_size]
        parts.append(dict(part='%s_p%03d' % (tag, len(parts)), tag=tag, n=len(part),
                          samples=';'.join(m['sample'] for m in part)))
    if any(c['tag'] == tag for c in old_chunks):
        sys.exit('tag %s already planned' % tag)
    write(w / 'chunks.tsv', old_chunks + chunks, ['chunk', 'tag', 'family', 'n', 'bp', 'samples'])
    write(w / 'parts.tsv', old_parts + parts, ['part', 'tag', 'n', 'samples'])
    info = dict(version=VERSION, tag=tag, batches=batches, manifest=a.manifest, manifest_sha256=ub.sha_file(a.manifest),
                assemblies=len(ms), already_planned=len(ms) - len(new), new=len(new), seed_chunks=len(chunks),
                genome_parts=len(parts), chunk_size=a.chunk_size, part_size=a.part_size, params=PARAMS)
    (w / ('plan_%s.json' % tag)).write_text(json.dumps(info, indent=1) + '\n')
    print(json.dumps({k: v for k, v in info.items() if k != 'params'}, indent=1))


def unit(w, table, key, idx, tag):
    us = [r for r in rows(w / table) if r['tag'] == tag]
    if idx >= len(us):
        return None
    return us[idx]


# ---------------------------------------------------------------- seeds

def merge_loci(segs, frac):
    """Greedy merge of [lo, hi) on one sequence: overlap >= frac of both -> one locus (outer extent)."""
    out = []
    for lo, hi in sorted(segs):
        if out:
            plo, phi = out[-1]
            ov = min(hi, phi) - max(lo, plo)
            if ov >= frac * (hi - lo) and ov >= frac * (phi - plo):
                out[-1] = (min(lo, plo), max(hi, phi))
                continue
        out.append((lo, hi))
    return out


def seed_records(paf, p):
    """Seed segments from minimap2 PAF lines: (name, lo, hi) for both sides of each qualifying record."""
    for line in paf:
        f = line.split('\t')
        qn, ql, qs, qe, tn, tl, ts, te, nm, al = f[0], int(f[1]), int(f[2]), int(f[3]), f[5], int(f[6]), int(f[7]), \
            int(f[8]), int(f[9]), int(f[10])
        if not (p['seed_min'] <= qe - qs <= p['seed_max'] and p['seed_min'] <= te - ts <= p['seed_max']):
            continue
        if al == 0 or nm / al < p['seed_identity']:
            continue
        yield qn, qs, qe
        yield tn, ts, te


def cluster(fa, out_prefix, tmp, threads, p):
    run([MMSEQS, 'easy-cluster', str(fa), str(out_prefix), str(tmp), '--min-seq-id', str(p['cluster_identity']),
         '-c', str(p['cluster_coverage']), '--cov-mode', '1', '--cluster-mode', '2', '--threads', str(threads)],
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def cmd_seeds(a):
    w = Path(a.work)
    c = unit(w, 'chunks.tsv', 'chunk', a.index, a.tag)
    if c is None:
        return print('no chunk %d' % a.index)
    out = w / 'seeds' / c['chunk']
    if (out / 'done.json').exists():
        return print('done', out)
    out.mkdir(parents=True, exist_ok=True)
    man = {m['sample']: m for m in rows(a.manifest)}
    tmp = Path(tempfile.mkdtemp(prefix='dseed_', dir=os.environ.get('TMPDIR')))
    t0 = time.time()
    try:
        names, bad = {}, []
        fa = tmp / 'chunk.fa'
        with open(fa, 'w') as f:
            for i, smp in enumerate(c['samples'].split(';')):
                g = verified_fasta(man[smp])
                if g is None:
                    bad.append(smp)
                    continue
                for ctg, s in g.items():
                    n = 'g%d_%d' % (i, len(names))
                    names[n] = (smp, ctg)
                    f.write('>%s\n%s\n' % (n, s))
        genomes = {}
        paf = tmp / 'ava.paf'
        t1 = time.time()
        with open(paf, 'w') as fo:
            run([MINIMAP2, '-x', PARAMS['seed_preset'], '-c', '-DP', '--dual=no',
                 '-t', str(a.threads), str(fa), str(fa)], stdout=fo, stderr=subprocess.DEVNULL)
        t_map = time.time() - t1
        segs, nrec = defaultdict(list), 0
        with open(paf) as f:
            for n, lo, hi in seed_records(f, PARAMS):
                segs[n].append((lo, hi))
                nrec += 1
        seqs = uc.read_fasta(fa)
        loci = []
        with open(tmp / 'seeds.fa', 'w') as f:
            for n in sorted(segs):
                smp, ctg = names[n]
                for lo, hi in merge_loci(segs[n], PARAMS['merge_overlap']):
                    sid = '%s_%d' % (c['chunk'], len(loci))
                    loci.append(dict(seed=sid, sample=smp, contig=ctg, lo=lo, hi=hi))
                    f.write('>%s\n%s\n' % (sid, seqs[n][lo:hi]))
        del seqs
        reps = 0
        if loci:
            cluster(tmp / 'seeds.fa', tmp / 'cl', tmp / 'mm', a.threads, PARAMS)
            shutil.copy(tmp / 'cl_rep_seq.fasta', out / 'reps.fa')
            rep_of = {}
            for line in open(tmp / 'cl_cluster.tsv'):
                r, m = line.split()
                rep_of[m] = r
            for x in loci:
                x['rep'] = rep_of.get(x['seed'], '')
            reps = len(set(rep_of.values()))
        else:
            (out / 'reps.fa').write_text('')
        write(out / 'loci.tsv.gz', loci, ['seed', 'sample', 'contig', 'lo', 'hi', 'rep'])
        info = dict(chunk=c['chunk'], assemblies=c['n'], bp=c['bp'], sha256_mismatch=bad, seed_records=nrec // 2,
                    seed_loci=len(loci), representatives=reps, minimap2_s=round(t_map, 1),
                    wall_s=round(time.time() - t0, 1), params=digest_params(), host=os.uname()[1])
        (out / 'done.json').write_text(json.dumps(info, indent=1) + '\n')
        print(json.dumps(info))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


# ---------------------------------------------------------------- library

def read_fa(path):
    return uc.read_fasta(path) if Path(path).exists() and Path(path).stat().st_size else {}


def cmd_library(a):
    w = Path(a.work)
    ld = w / 'library'
    ld.mkdir(parents=True, exist_ok=True)
    shards = sorted(int(p.stem.split('_')[1]) for p in ld.glob('lib_*.fa'))
    k = (max(shards) + 1) if shards else 0
    chunks = [r['chunk'] for r in rows(w / 'chunks.tsv') if r['tag'] == a.tag]
    miss = [c for c in chunks if not (w / 'seeds' / c / 'done.json').exists()]
    if miss:
        sys.exit('seed chunks not done: %s' % miss[:5])
    tmp = Path(tempfile.mkdtemp(prefix='dlib_', dir=os.environ.get('TMPDIR')))
    try:
        new = {}
        for c in chunks:
            new.update(read_fa(w / 'seeds' / c / 'reps.fa'))
        with open(tmp / 'new.fa', 'w') as f:
            for n, s in new.items():
                f.write('>%s\n%s\n' % (n, s))
        known = 0
        if shards:
            with open(tmp / 'old.fa', 'w') as f:
                for s in shards:
                    f.write(open(ld / ('lib_%d.fa' % s)).read())
            run([MMSEQS, 'easy-search', str(tmp / 'new.fa'), str(tmp / 'old.fa'), str(tmp / 'hits.m8'), str(tmp / 'mm'),
                 '--search-type', '3', '--threads', str(a.threads), '--format-output',
                 'query,target,fident,qcov,tcov'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            hit = set()
            for line in open(tmp / 'hits.m8'):
                q, t, fid, qc, tc = line.split('\t')
                if float(fid) >= PARAMS['cluster_identity'] and float(qc) >= PARAMS['cluster_coverage']:
                    hit.add(q)
            known = len(hit)
            with open(tmp / 'rest.fa', 'w') as f:
                for n, s in new.items():
                    if n not in hit:
                        f.write('>%s\n%s\n' % (n, s))
            src = tmp / 'rest.fa'
        else:
            src = tmp / 'new.fa'
        reps = {}
        if src.stat().st_size:
            cluster(src, tmp / 'cl', tmp / 'mm2', a.threads, PARAMS)
            reps = uc.read_fasta(tmp / 'cl_rep_seq.fasta')
        names = {}
        with open(ld / ('lib_%d.fa.tmp' % k), 'w') as f:
            for i, n in enumerate(sorted(reps)):
                names[n] = 'L%d_%d' % (k, i)
                f.write('>%s\n%s\n' % (names[n], reps[n]))
        os.replace(ld / ('lib_%d.fa.tmp' % k), ld / ('lib_%d.fa' % k))
        write(ld / ('lib_%d_sources.tsv.gz' % k), [dict(library_id=v, chunk_rep=n) for n, v in names.items()],
              ['library_id', 'chunk_rep'])
        info = dict(shard=k, tag=a.tag, chunk_representatives=len(new), matched_existing_library=known,
                    new_library_sequences=len(reps), bp=sum(map(len, reps.values())), params=digest_params())
        (ld / ('lib_%d.json' % k)).write_text(json.dumps(info, indent=1) + '\n')
        print(json.dumps(info))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


# ---------------------------------------------------------------- occurrences

OCC_COLS = ['library_id', 'sample', 'contig', 'contig_len', 'lo', 'hi', 'strand', 'identity', 'lib_lo', 'lib_hi',
            'lib_len', 'win_lo', 'win_hi', 'window']


def occ_tasks(w):
    shards = sorted(int(p.stem.split('_')[1]) for p in (w / 'library').glob('lib_*.fa'))
    parts = [r['part'] for r in rows(w / 'parts.tsv')]
    return [(k, p) for k in shards for p in parts]


def cmd_occurrences(a):
    w = Path(a.work)
    tasks = occ_tasks(w)
    if a.list:  # array indices of the tasks not done (task order is fixed: shard-major, then part)
        return print(','.join(str(i) for i, t in enumerate(tasks)
                              if not (w / 'occ' / ('lib%d' % t[0]) / t[1] / 'done.json').exists()))
    if a.index >= len(tasks):
        return print('no task %d' % a.index)
    k, part = tasks[a.index]
    out = w / 'occ' / ('lib%d' % k) / part
    if (out / 'done.json').exists():
        return print('done', out)
    out.mkdir(parents=True, exist_ok=True)
    man = {m['sample']: m for m in rows(a.manifest)}
    samples = [r for r in rows(w / 'parts.tsv') if r['part'] == part][0]['samples'].split(';')
    lib = w / 'library' / ('lib_%d.fa' % k)
    tmp = Path(tempfile.mkdtemp(prefix='docc_', dir=os.environ.get('TMPDIR')))
    t0 = time.time()
    n, bad, per = 0, [], Counter()
    try:
        B = PARAMS['buckets']
        tmpf = {b: Path(str(out / ('b%02d.tsv.gz' % b)) + '.tmp') for b in range(B)}
        fos = {b: gzip.open(tmpf[b], 'wt') for b in range(B)}
        for fo in fos.values():
            fo.write('\t'.join(OCC_COLS) + '\n')
        if lib.stat().st_size:
            run([MINIMAP2, '-x', PARAMS['seed_preset'], '-d', str(tmp / 'lib.mmi'), str(lib)],
                stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            W = PARAMS['window']
            if True:
                for smp in samples:
                    g = verified_fasta(man[smp])
                    if g is None:
                        bad.append(smp)
                        continue
                    fa = tmp / 'g.fa'
                    with open(fa, 'w') as f:
                        for ctg, s in g.items():
                            f.write('>%s\n%s\n' % (ctg, s))
                    p = run([MINIMAP2, '-x', PARAMS['seed_preset'], '-c', '-P',
                             '-t', str(a.threads), str(tmp / 'lib.mmi'), str(fa)], capture_output=True, text=True)
                    for line in p.stdout.splitlines():
                        f = line.split('\t')
                        ctg, cl, qs, qe, st, lid, ll, ts, te, nm, al = f[0], int(f[1]), int(f[2]), int(f[3]), f[4], \
                            f[5], int(f[6]), int(f[7]), int(f[8]), int(f[9]), int(f[10])
                        if qe - qs < PARAMS['occ_min'] or al == 0 or nm / al < PARAMS['seed_identity']:
                            continue
                        lo, hi = max(0, qs - W), min(cl, qe + W)
                        fos[shard_of(lid, B)].write('\t'.join(map(str, (lid, smp, ctg, cl, qs, qe, st, round(nm / al, 4),
                                                                        ts, te, ll, lo, hi, g[ctg][lo:hi]))) + '\n')
                        n += 1
                        per[smp] += 1
        for b in range(B):
            fos[b].close()
            os.replace(tmpf[b], out / ('b%02d.tsv.gz' % b))
        info = dict(library_shard=k, part=part, assemblies=len(samples), sha256_mismatch=bad, hits=n,
                    assemblies_with_hits=len(per), wall_s=round(time.time() - t0, 1), params=digest_params(),
                    host=os.uname()[1])
        (out / 'done.json').write_text(json.dumps(info, indent=1) + '\n')
        print(json.dumps(info))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


# ---------------------------------------------------------------- pairs

def shard_of(lid, n):
    return int(hashlib.sha1(lid.encode()).hexdigest()[:8], 16) % n


def merge_occurrences(occs):
    """Hits of one library sequence on one contig that overlap are one occurrence (outer extent, best identity);
    a hit on the other strand overlapping it is kept apart (inverted copies are separate occurrences)."""
    out = []
    for o in sorted(occs, key=lambda o: (o['sample'], o['contig'], o['strand'], int(o['lo']))):
        if out:
            p = out[-1]
            if (p['sample'], p['contig'], p['strand']) == (o['sample'], o['contig'], o['strand']) and \
                    int(o['lo']) < int(p['hi']):
                if int(o['hi']) > int(p['hi']):
                    # extend; window must still contain the seed: rebuild from the two windows (contiguous, overlapping)
                    a_lo, a_hi, b_lo, b_hi = int(p['win_lo']), int(p['win_hi']), int(o['win_lo']), int(o['win_hi'])
                    if b_lo <= a_hi:
                        p['window'] = p['window'] + o['window'][a_hi - b_lo:] if b_hi > a_hi else p['window']
                        p['win_hi'] = str(max(a_hi, b_hi))
                    p['hi'] = o['hi']
                p['identity'] = str(max(float(p['identity']), float(o['identity'])))
                p['merged_hits'] = p.get('merged_hits', 1) + 1
                continue
        out.append(dict(o))
    return out


def oriented(o):
    """Window in library orientation and the seed's [c0, c1) in it."""
    wl, wh, lo, hi = int(o['win_lo']), int(o['win_hi']), int(o['lo']), int(o['hi'])
    s = o['window']
    if o['strand'] == '+':
        return s, lo - wl, hi - wl
    return cr.revcomp(s), wh - hi, wh - lo


def occ_id(o):
    return 'dna|%s|%s|%s-%s|%s' % (o['sample'], o['contig'], o['lo'], o['hi'], o['strand'])


def design(occs, p):
    """Pairs (i, j) of one library sequence: all, or the occurrence-cover design for high-copy sequences."""
    n = len(occs)
    total = n * (n - 1) // 2
    if total <= p['full_pairs']:
        return [(i, j) for i in range(n) for j in range(i + 1, n)], total, 'all_pairs'
    K = p['design_partners']
    rank = sorted(range(n), key=lambda i: hashlib.sha1(occ_id(occs[i]).encode()).hexdigest())
    pos = {i: r for r, i in enumerate(rank)}
    by_smp = defaultdict(list)
    for i, o in enumerate(occs):
        by_smp[o['sample']].append(i)
    chosen = set()
    for i, o in enumerate(occs):
        same = [j for j in by_smp[o['sample']] if j != i][:K // 2]
        need = K - len(same)
        cross = []
        r = pos[i]
        step = 1
        while len(cross) < need and step < n:
            j = rank[(r + step) % n]
            if occs[j]['sample'] != o['sample']:
                cross.append(j)
            step += 1
        for j in same + cross:
            chosen.add((min(i, j), max(i, j)))
    return sorted(chosen), total, 'occurrence_cover_design'


def extend_pair(A, B, hsps, p):
    """Candidate X of one pair from the HSPs of A's vs B's oriented windows (core_review chain / extensions)."""
    wa, a0, a1 = A['_w']
    wb, b0, b1 = B['_w']
    la, lb = A.get('_len', None), B.get('_len', None)
    la = len(wa) if la is None else la
    lb = len(wb) if lb is None else lb
    pieces = cr.chain([(h['qs'], h['qe'], h['ss'], h['se']) for h in hsps], a0, a1, b0, b1, p['chain_gap'])
    if not pieces:
        return dict(status='no_shared_block_at_seed')
    e = cr.extensions(pieces, a0, a1, la, b0, b1, lb, p['edge'])
    res = {}
    for end in ('up', 'down'):
        cov = e[end + '_cover']
        if not cov.startswith('covered'):
            res[end] = cov
        elif e[end + '_reaches_end']:
            res[end] = 'shared_to_window_or_contig_end'
        else:
            res[end] = 'resolved'
    return dict(status='resolved' if res['up'] == res['down'] == 'resolved' else 'unresolved', up=res['up'],
                down=res['down'], a_up=e['a_up'], a_down=e['a_down'], b_up=e['b_up'], b_down=e['b_down'],
                chain_pieces=len(pieces))


def genomic(o, up, down):
    """Genome [lo, hi) of the seed extended by up / down bp (library orientation)."""
    lo, hi = int(o['lo']), int(o['hi'])
    return (lo - up, hi + down) if o['strand'] == '+' else (lo - down, hi + up)


PAIR_COLS = ['pair_id', 'library_id', 'pair_scope', 'design', 'a_occurrence', 'b_occurrence', 'a_sample', 'b_sample',
             'a_contig', 'b_contig', 'a_strand', 'b_strand', 'a_seed', 'b_seed', 'seed_identity_a', 'seed_identity_b',
             'status', 'boundary_up', 'boundary_down', 'a_ext_up', 'a_ext_down', 'b_ext_up', 'b_ext_down',
             'a_x', 'b_x', 'x_len_a', 'x_len_b', 'chain_pieces']


def cmd_pairs(a):
    w = Path(a.work)
    B = PARAMS['buckets']
    bucket, sub = a.index % B, a.index // B
    out = w / 'pairs' / ('shard_%04d' % a.index)
    if (out / 'done.json').exists() and not a.force:
        return print('done', out)
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    occ = defaultdict(list)
    for f in sorted((w / 'occ').glob('lib*/*/b%02d.tsv.gz' % bucket)):
        with gzip.open(f, 'rt') as fh:
            h = fh.readline().rstrip('\n').split('\t')
            for line in fh:
                v = line.rstrip('\n').split('\t')
                if a.subs == 1 or shard_of('sub|' + v[0], a.subs) == sub:
                    occ[v[0]].append(dict(zip(h, v)))
    P, budget = [], []
    tmp = Path(tempfile.mkdtemp(prefix='dpair_', dir=os.environ.get('TMPDIR')))
    try:
        for lid in sorted(occ):
            os_ = merge_occurrences(occ[lid])
            os_.sort(key=occ_id)
            for o in os_:
                o['_w'] = oriented(o)
                o['_id'] = occ_id(o)
            pairs, total, how = design(os_, PARAMS)
            budget.append(dict(library_id=lid, occurrences=len(os_), hits=len(occ[lid]), pairs_total=total,
                               pairs_formed=len(pairs), pairs_deferred=total - len(pairs), design=how,
                               assemblies=len({o['sample'] for o in os_})))
            if not pairs:
                continue
            need = sorted({i for pr in pairs for i in pr})
            if how == 'all_pairs':
                hits = cr.blast('blastn', {str(i): os_[i]['_w'][0] for i in need}, 'blastn')
                get = lambda i, j: hits.get((str(i), str(j)), [])
            else:
                get = None
            for i, j in pairs:
                A, B = os_[i], os_[j]
                if get is not None:
                    hs = get(i, j)
                else:
                    hs = cr.blast('blastn', {'0': A['_w'][0], '1': B['_w'][0]}, 'blastn').get(('0', '1'), [])
                e = extend_pair(A, B, hs, PARAMS)
                r = dict(library_id=lid, pair_scope='same_assembly' if A['sample'] == B['sample'] else 'cross_assembly',
                         design=how, a_occurrence=A['_id'], b_occurrence=B['_id'], a_sample=A['sample'],
                         b_sample=B['sample'], a_contig=A['contig'], b_contig=B['contig'], a_strand=A['strand'],
                         b_strand=B['strand'], a_seed='%s-%s' % (A['lo'], A['hi']), b_seed='%s-%s' % (B['lo'], B['hi']),
                         seed_identity_a=A['identity'], seed_identity_b=B['identity'], status=e['status'],
                         boundary_up=e.get('up', e['status']), boundary_down=e.get('down', e['status']),
                         chain_pieces=e.get('chain_pieces', ''))
                for s, O in (('a', A), ('b', B)):
                    up, down = e.get(s + '_up'), e.get(s + '_down')
                    r[s + '_ext_up'] = '' if up is None else up
                    r[s + '_ext_down'] = '' if down is None else down
                    if e['status'] == 'resolved':
                        lo, hi = genomic(O, up, down)
                        r[s + '_x'] = '%d-%d' % (lo, hi)
                        r['x_len_' + s] = hi - lo
                key = (r['a_sample'], r['a_contig'], r.get('a_x') or r['a_seed'], r['a_strand'], r['b_sample'],
                       r['b_contig'], r.get('b_x') or r['b_seed'], r['b_strand'])
                r['pair_id'] = 'D' + hashlib.sha1('|'.join(key).encode()).hexdigest()[:14]
                P.append(r)
        write(out / 'pairs.tsv.gz', P, PAIR_COLS)
        write(out / 'budget.tsv.gz', budget, ['library_id', 'occurrences', 'hits', 'assemblies', 'pairs_total',
                                              'pairs_formed', 'pairs_deferred', 'design'])
        st = Counter(r['status'] for r in P)
        info = dict(shard=a.index, bucket=bucket, sub=sub, subs=a.subs, library_sequences=len(occ), pairs=len(P), status=dict(st),
                    pairs_deferred=sum(b['pairs_deferred'] for b in budget),
                    designed_sequences=sum(b['design'] != 'all_pairs' for b in budget),
                    wall_s=round(time.time() - t0, 1), params=digest_params(), host=os.uname()[1])
        (out / 'done.json').write_text(json.dumps(info, indent=1) + '\n')
        print(json.dumps(info))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


# ---------------------------------------------------------------- enumerate / align (v0.2: alignment de-duplicated)

SRC_SIDE = ['sample', 'contig', 'clen', 'lo', 'hi', 'strand', 'identity', 'wlo', 'whi']
SRC_COLS = ['key', 'library_id', 'chunk', 'design'] + ['a_' + c for c in SRC_SIDE] + ['b_' + c for c in SRC_SIDE]


def grid_window(lo, hi, clen, p):
    """Alignment window of an occurrence: seed +/- W, widened outward to the window grid, clipped to the contig."""
    W, G = p['window'], p['window_grid']
    return (max(0, lo - W) // G) * G, min(clen, -(-(hi + W) // G) * G)


def align_key(A, B):
    """Canonical alignment key of a source pair: the two windows (sample, contig, lo, hi), first = smaller, and the
    relative orientation of the two copies. A-B and B-A, and a pair seen in either global orientation, have one key."""
    wa = (A['sample'], A['contig'], int(A['wlo']), int(A['whi']))
    wb = (B['sample'], B['contig'], int(B['wlo']), int(B['whi']))
    rel = '+' if A['strand'] == B['strand'] else '-'
    first, second = (wa, wb) if wa <= wb else (wb, wa)
    k = '%s|%s|%d|%d||%s|%s|%d|%d||%s' % (first + second + (rel,))
    return hashlib.sha1(k.encode()).hexdigest()[:20], first, second, rel


def key_shard(second, rel, nk):
    """Shard of an alignment key: by its subject window (second window and relative orientation), so that all keys
    of one subject are aligned in one blastn call (v0.2.1; v0.2.0 sharded by key hash: ~1.3 queries per call)."""
    return int(hashlib.sha1(('%s|%s|%d|%d|%s' % (second + (rel,))).encode()).hexdigest()[:8], 16) % nk


def pair_chunks(n, occs, p):
    """Deterministic chunks of the n*(n-1)/2 pairs of one library sequence. Up to pair_chunk pairs: one chunk. Above:
    chunk 0 = the occurrence-cover design (DESIGN_PARTNERS partners per occurrence; priority only), then all other
    pairs in (i, j) order in chunks of pair_chunk. Every pair is in exactly one chunk; the chunk size is not a cap."""
    total = n * (n - 1) // 2
    C = p['pair_chunk']
    if total <= C:
        yield 0, 'all_pairs', ((i, j) for i in range(n) for j in range(i + 1, n)), total
        return
    first, _, _ = design(occs, p)
    fs = set(first)
    yield 0, 'cover_design', iter(first), len(first)
    rest = total - len(first)
    nchunks = -(-rest // C)

    def gen(c):
        k = 0
        lo, hi = c * C, (c + 1) * C
        for i in range(n):
            for j in range(i + 1, n):
                if (i, j) in fs:
                    continue
                if lo <= k < hi:
                    yield i, j
                k += 1
                if k >= hi:
                    return
    for c in range(nchunks):
        yield c + 1, 'remaining', gen(c), min(C, rest - c * C)


def cmd_enumerate(a):
    """Source pairs of every library sequence of one occurrence bucket, written to key shards (no alignment).
    --max-pairs bounds this run; chunks beyond it are recorded as deferred and enumerated by a later --gen."""
    w = Path(a.work)
    B, NK = PARAMS['buckets'], PARAMS['key_shards']
    bucket = a.index
    if bucket >= B:
        return print('no bucket %d' % bucket)
    ed = w / 'enum' / ('g%d' % a.gen) / ('b%02d' % bucket)
    if (ed / 'done.json').exists():
        return print('done', ed)
    done_chunks = set()
    for g in range(a.gen):
        pf = w / 'enum' / ('g%d' % g) / ('b%02d' % bucket) / 'progress.tsv.gz'
        if not pf.exists():
            sys.exit('generation %d of bucket %d not done' % (g, bucket))
        for r in rows(pf):
            if r['status'] == 'enumerated':  # a chunk is the same only for the same occurrence count of its sequence
                done_chunks.add((r['library_id'], int(r['chunk']), int(r['occurrences'])))
    ed.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    occ = defaultdict(list)
    for f in sorted((w / 'occ').glob('lib*/*/b%02d.tsv.gz' % bucket)):
        with gzip.open(f, 'rt') as fh:
            h = fh.readline().rstrip('\n').split('\t')
            for line in fh:
                v = line.rstrip('\n').split('\t')
                d = dict(zip(h, v))
                d.pop('window')
                occ[v[0]].append(d)
    tmpf = {k: ed / ('.tmp_k%03d.tsv.gz' % k) for k in range(NK)}
    fos = {k: gzip.open(tmpf[k], 'wt', compresslevel=3) for k in range(NK)}
    for fo in fos.values():
        fo.write('\t'.join(SRC_COLS) + '\n')
    prog, n_enum, n_def, n_skip, budget = [], 0, 0, 0, a.max_pairs
    try:
        for lid in sorted(occ):
            os_ = merge_occurrences([dict(o, window='') for o in occ[lid]])
            os_.sort(key=occ_id)
            for o in os_:
                o['wlo'], o['whi'] = grid_window(int(o['lo']), int(o['hi']), int(o['contig_len']), PARAMS)
            n = len(os_)
            for c, how, it, size in pair_chunks(n, os_, PARAMS):
                if (lid, c, n) in done_chunks:
                    n_skip += size
                    continue
                if budget and n_enum + size > budget and n_enum > 0:
                    prog.append(dict(library_id=lid, chunk=c, design=how, pairs=size, occurrences=n, status='deferred'))
                    n_def += size
                    continue
                m = 0
                for i, j in it:
                    A, Bo = os_[i], os_[j]
                    key, _, second, rel = align_key(A, Bo)
                    vals = [key, lid, c, how]
                    for O in (A, Bo):
                        vals += [O['sample'], O['contig'], O['contig_len'], O['lo'], O['hi'], O['strand'],
                                 O['identity'], O['wlo'], O['whi']]
                    fos[key_shard(second, rel, NK)].write('\t'.join(map(str, vals)) + '\n')
                    m += 1
                prog.append(dict(library_id=lid, chunk=c, design=how, pairs=m, occurrences=n, status='enumerated'))
                n_enum += m
        for k in range(NK):
            fos[k].close()
            os.replace(tmpf[k], ed / ('k%03d.tsv.gz' % k))
        write(ed / 'progress.tsv.gz', prog, ['library_id', 'chunk', 'design', 'pairs', 'occurrences', 'status'])
        info = dict(bucket=bucket, gen=a.gen, library_sequences=len(occ), pairs_enumerated=n_enum,
                    pairs_deferred=n_def, pairs_done_earlier=n_skip, max_pairs=a.max_pairs,
                    chunks_deferred=sum(r['status'] == 'deferred' for r in prog),
                    wall_s=round(time.time() - t0, 1), params=digest_params(), host=os.uname()[1])
        (ed / 'done.json').write_text(json.dumps(info, indent=1) + '\n')
        print(json.dumps(info))
    finally:
        for fo in fos.values():
            if not fo.closed:
                fo.close()


def _blast_group(job):
    """One subject window against all its query windows (blastn as core_review.blast: task blastn, plus strand, dust
    off, evalue 1e-5). One subject per call, so every (query, subject) result equals a separate pairwise search."""
    subj, queries = job
    with tempfile.NamedTemporaryFile('w', suffix='.fa') as q, tempfile.NamedTemporaryFile('w', suffix='.fa') as s:
        for k, seq in queries:
            q.write('>%s\n%s\n' % (k, seq))
        s.write('>s\n%s\n' % subj)
        q.flush()
        s.flush()
        cmd = ['blastn', '-query', q.name, '-subject', s.name, '-outfmt', cr.OUTFMT, '-evalue', '1e-5',
               '-max_target_seqs', '500', '-task', 'blastn', '-strand', 'plus', '-dust', 'no']
        out = subprocess.run(cmd, check=True, capture_output=True, text=True).stdout
    res = {k: [] for k, _ in queries}
    for line in out.splitlines():
        f = line.split('\t')
        res[f[0]].append((int(f[4]), int(f[5]), int(f[6]), int(f[7])))
    return res


def to_source(hsps, L1, L2, a_is_first, flip):
    """HSPs of the canonical frame (first window forward, second in relative orientation) in a source's frame
    (A and B each in its library orientation): swap if A is the second window, reverse-complement both if needed."""
    out = []
    for qs, qe, ss, se in hsps:
        if not a_is_first:
            qs, qe, ss, se = ss, se, qs, qe
        out.append((qs, qe, ss, se))
    La, Lb = (L1, L2) if a_is_first else (L2, L1)
    if flip:
        out = [(La - qe + 1, La - qs + 1, Lb - se + 1, Lb - ss + 1) for qs, qe, ss, se in out]
    return [dict(qs=a, qe=b, ss=c, se=d) for a, b, c, d in out]


def seed_in_window(o):
    """[c0, c1) of the seed in the occurrence's grid window, in library orientation."""
    wl, wh, lo, hi = int(o['wlo']), int(o['whi']), int(o['lo']), int(o['hi'])
    return (lo - wl, hi - wl) if o['strand'] == '+' else (wh - hi, wh - lo)


ALIGN_COLS = PAIR_COLS + ['align_key', 'chunk', 'alignment']


def cmd_align(a):
    """Per key shard: every distinct alignment key once (cached across generations), then the extension of every
    source pair from the HSPs in its own frame. Each source keeps its own seeds, X and library sequence."""
    w = Path(a.work)
    NK = PARAMS['key_shards']
    s = a.index
    if s >= NK:
        return print('no shard %d' % s)
    out = w / 'align' / ('g%d' % a.gen) / ('k%03d' % s)
    if (out / 'done.json').exists():
        return print('done', out)
    nb = PARAMS['buckets']
    miss = [b for b in range(nb) if not (w / 'enum' / ('g%d' % a.gen) / ('b%02d' % b) / 'done.json').exists()]
    if miss:
        sys.exit('enumerate g%d not done for buckets %s' % (a.gen, miss[:5]))
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    srcs = []
    for b in range(nb):
        srcs.extend(rows(w / 'enum' / ('g%d' % a.gen) / ('b%02d' % b) / ('k%03d.tsv.gz' % s)))
    cache_dir = w / 'align' / 'cache' / ('k%03d' % s)
    cache = {}
    for f in sorted(cache_dir.glob('hsps_g*.tsv.gz')):
        for r in rows(f):
            cache[r['key']] = (int(r['qlen']), int(r['slen']),
                               [tuple(map(int, x.split(','))) for x in r['hsps'].split(';') if x])
    keys = {}
    for r in srcs:
        if r['key'] in cache or r['key'] in keys:
            continue
        A = {k: r['a_' + k] for k in SRC_SIDE}
        Bo = {k: r['b_' + k] for k in SRC_SIDE}
        key, first, second, rel = align_key(A, Bo)
        if key != r['key']:
            sys.exit('key mismatch %s %s' % (key, r['key']))
        keys[key] = (first, second, rel)
    # window sequences, one verified genome at a time
    man = {m['sample']: m for m in rows(a.manifest)}
    need = defaultdict(set)
    for first, second, rel in keys.values():
        need[first[0]].add(first)
        need[second[0]].add(second)
    seqs, bad = {}, []
    for smp in sorted(need):
        g = verified_fasta(man[smp])
        if g is None:
            bad.append(smp)
            continue
        for wnd in need[smp]:
            seqs[wnd] = g[wnd[1]][wnd[2]:wnd[3]]
        del g
    if bad:
        sys.exit('sha256 mismatch: %s' % bad)
    groups = defaultdict(list)
    for key, (first, second, rel) in keys.items():
        groups[(second, rel)].append((key, seqs[first]))
    jobs = []
    for (second, rel), qs in sorted(groups.items()):
        subj = seqs[second] if rel == '+' else cr.revcomp(seqs[second])
        jobs.append((subj, sorted(qs)))
    t1 = time.time()
    new = {}
    if jobs:
        import multiprocessing as mp
        with mp.Pool(a.threads) as pool:
            for res in pool.imap_unordered(_blast_group, jobs, chunksize=1):
                new.update(res)
    t_blast = time.time() - t1
    cache_dir.mkdir(parents=True, exist_ok=True)
    crow = []
    for key, (first, second, rel) in keys.items():
        hs = new.get(key, [])
        L1, L2 = first[3] - first[2], second[3] - second[2]
        cache[key] = (L1, L2, hs)
        crow.append(dict(key=key, first='%s|%s|%d-%d' % first, second='%s|%s|%d-%d' % second, rel=rel, qlen=L1,
                         slen=L2, hsps=';'.join('%d,%d,%d,%d' % h for h in hs)))
    write(cache_dir / ('hsps_g%d.tsv.gz' % a.gen), crow, ['key', 'first', 'second', 'rel', 'qlen', 'slen', 'hsps'])
    P = []
    for r in srcs:
        A = {k: r['a_' + k] for k in SRC_SIDE}
        Bo = {k: r['b_' + k] for k in SRC_SIDE}
        _, first, second, rel = align_key(A, Bo)
        L1, L2, hs = cache[r['key']]
        wa = (A['sample'], A['contig'], int(A['wlo']), int(A['whi']))
        a_first = wa == first
        if first == second:  # both copies in one window: A is the query by occurrence order
            a_first = True
        o_first = A['strand'] if a_first else Bo['strand']
        flip = o_first != '+'
        if not a_first:
            flip = A['strand'] != rel
        hsps = to_source(hs, L1, L2, a_first, flip)
        for O in (A, Bo):
            c0, c1 = seed_in_window(O)
            O['_w'] = (None, c0, c1)
            O['_len'] = int(O['whi']) - int(O['wlo'])
            O['contig_len'] = O['clen']
            O['_id'] = occ_id(O)
        e = extend_pair(A, Bo, hsps, PARAMS)
        pr = dict(library_id=r['library_id'], design=r['design'], chunk=r['chunk'], align_key=r['key'],
                  alignment='computed' if r['key'] in keys else 'cached',
                  pair_scope='same_assembly' if A['sample'] == Bo['sample'] else 'cross_assembly',
                  a_occurrence=A['_id'], b_occurrence=Bo['_id'], a_sample=A['sample'], b_sample=Bo['sample'],
                  a_contig=A['contig'], b_contig=Bo['contig'], a_strand=A['strand'], b_strand=Bo['strand'],
                  a_seed='%s-%s' % (A['lo'], A['hi']), b_seed='%s-%s' % (Bo['lo'], Bo['hi']),
                  seed_identity_a=A['identity'], seed_identity_b=Bo['identity'], status=e['status'],
                  boundary_up=e.get('up', e['status']), boundary_down=e.get('down', e['status']),
                  chain_pieces=e.get('chain_pieces', ''))
        for sd, O in (('a', A), ('b', Bo)):
            up, down = e.get(sd + '_up'), e.get(sd + '_down')
            pr[sd + '_ext_up'] = '' if up is None else up
            pr[sd + '_ext_down'] = '' if down is None else down
            if e['status'] == 'resolved':
                lo, hi = genomic(O, up, down)
                pr[sd + '_x'] = '%d-%d' % (lo, hi)
                pr['x_len_' + sd] = hi - lo
        k2 = (pr['a_sample'], pr['a_contig'], pr.get('a_x') or pr['a_seed'], pr['a_strand'], pr['b_sample'],
              pr['b_contig'], pr.get('b_x') or pr['b_seed'], pr['b_strand'])
        pr['pair_id'] = 'D' + hashlib.sha1('|'.join(k2).encode()).hexdigest()[:14]
        P.append(pr)
    write(out / 'pairs.tsv.gz', P, ALIGN_COLS)
    info = dict(shard=s, gen=a.gen, sources=len(srcs), distinct_keys=len({r['key'] for r in srcs}),
                keys_computed=len(keys), keys_cached=len({r['key'] for r in srcs} - set(keys)),
                blast_calls=len(jobs), windows=len(seqs), blast_s=round(t_blast, 1),
                status=dict(Counter(r['status'] for r in P)), wall_s=round(time.time() - t0, 1),
                params=digest_params(), host=os.uname()[1])
    (out / 'done.json').write_text(json.dumps(info, indent=1) + '\n')
    print(json.dumps(info))


# ---------------------------------------------------------------- collect

def length_class(r):
    """Target length of a pair's candidate X (both sides): in_target / out_of_target; unresolved when not known."""
    if r['status'] != 'resolved':
        return 'unresolved'
    lo, hi = PARAMS['target_min'], PARAMS['target_max']
    ok = all(lo <= int(r['x_len_' + s]) < hi for s in 'ab')
    return 'in_target' if ok else 'out_of_target'


def cmd_collect(a):
    """ucd_batch input table of the resolved pairs (candidate_insertion), one row per distinct pair of X intervals."""
    w = Path(a.work)
    if a.gens:  # v0.2: align shards of the given generations; deferred = enumerate chunks never enumerated
        gens = [int(g) for g in a.gens.split(',')]
        NK, nb = PARAMS['key_shards'], PARAMS['buckets']
        units = [w / 'align' / ('g%d' % g) / ('k%03d' % k) for g in gens for k in range(NK)]
        miss = [str(u) for u in units if not (u / 'done.json').exists()]
        if miss:
            sys.exit('align shards not done: %d (first %s)' % (len(miss), miss[:3]))
        enumerated = set()
        deferred_rows = {}
        for g in gens:
            for b in range(nb):
                for r in rows(w / 'enum' / ('g%d' % g) / ('b%02d' % b) / 'progress.tsv.gz'):
                    k = (r['library_id'], r['chunk'])
                    if r['status'] == 'enumerated':
                        enumerated.add(k)
                    else:
                        deferred_rows[k] = int(r['pairs'])
        pending = {k: v for k, v in deferred_rows.items() if k not in enumerated}
        shard_files = [(u / 'pairs.tsv.gz', None) for u in units]
        deferred_info = dict(chunks_pending=len(pending), pairs_pending=sum(pending.values()),
                             libraries_pending=len({k[0] for k in pending}))
    else:
        S = PARAMS['buckets'] * a.subs
        miss = [i for i in range(S) if not (w / 'pairs' / ('shard_%04d' % i) / 'done.json').exists()]
        if miss:
            sys.exit('pair shards not done: %d (first %s)' % (len(miss), miss[:5]))
        shard_files = [(w / 'pairs' / ('shard_%04d' % i) / 'pairs.tsv.gz', w / 'pairs' / ('shard_%04d' % i))
                       for i in range(S)]
        deferred_info = None
    man = {m['sample']: m for m in rows(a.manifest)}
    seen, dup, st, scope, n_all = {}, Counter(), Counter(), Counter(), 0
    out_target = set()
    deferred = 0
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    ci = {c: i for i, c in enumerate(ub.INPUT_COLS)}
    with gzip.open(out / 'dna_pairs_all.tsv.gz', 'wt') as fa:
        fa.write('\t'.join(PAIR_COLS + ['length_class']) + '\n')
        alignment = Counter()
        lclass = Counter()  # distinct resolved pairs by target length of both X
        for pf, d in shard_files:
            if d is not None:
                deferred += json.loads((d / 'done.json').read_text())['pairs_deferred']
            for r in rows(pf):
                alignment[r.get('alignment', '')] += 1
                n_all += 1
                lc = length_class(r)
                fa.write('\t'.join([str(r.get(c, '')) for c in PAIR_COLS] + [lc]) + '\n')
                st[r['status']] += 1
                if r['status'] != 'resolved':
                    continue
                if r['pair_id'] in seen or r['pair_id'] in out_target:
                    dup[r['pair_id']] += 1
                    continue
                lclass[lc] += 1
                if lc != 'in_target':  # length known and outside the target: not sent to the two checks
                    out_target.add(r['pair_id'])
                    continue
                scope[r['pair_scope']] += 1
                x = dict(pair_id=r['pair_id'], structure='dna_segment', core_key=r['library_id'],
                         pair_scope=r['pair_scope'], copy_class='', s6_grouping='', core_match_class='',
                         core_dna_identity='', core_dna_coverage_min='', protein_identity_min='', chunk='')
                for s in 'ab':
                    lo, hi = r[s + '_x'].split('-')
                    slo, shi = r[s + '_seed'].split('-')
                    x.update({s + '_occurrence': r[s + '_occurrence'], s + '_sample': r[s + '_sample'],
                              s + '_assembly': man[r[s + '_sample']]['assembly'], s + '_contig': r[s + '_contig'],
                              s + '_start0': slo, s + '_end0': shi, s + '_strand': r[s + '_strand'],
                              s + '_ext_up': r[s + '_ext_up'], s + '_ext_down': r[s + '_ext_down'],
                              s + '_elem_lo': lo, s + '_elem_hi': hi, s + '_boundary_up': 'resolved',
                              s + '_boundary_down': 'resolved', s + '_estimated_span': int(hi) - int(lo),
                              s + '_input_error': ''})
                if x['a_occurrence'] > x['b_occurrence']:  # canonical order of ucd_batch (A = smaller occurrence)
                    sw = {}
                    for f in ub.SIDE:
                        sw['a_' + f], sw['b_' + f] = x['b_' + f], x['a_' + f]
                    x.update(sw)
                seen[r['pair_id']] = [x.get(c, '') for c in ub.INPUT_COLS]
    with gzip.open(out / 'ucd_input.tsv.gz', 'wt') as f:
        f.write('\t'.join(ub.INPUT_COLS) + '\n')
        for pid in sorted(seen):
            f.write('\t'.join(str(v) for v in seen[pid]) + '\n')
    info = dict(version=VERSION, discovery_source='dna', params=PARAMS, pair_rows=n_all, status=dict(st),
                resolved_distinct=len(seen) + len(out_target), resolved_duplicate_rows=sum(dup.values()),
                target_length=dict(min=PARAMS['target_min'], max=PARAMS['target_max'], rule='both X, min <= len < max'),
                resolved_distinct_by_length=dict(lclass), sent_to_checks=len(seen), out_of_target=len(out_target),
                scope_sent_to_checks=dict(scope),
                pairs_deferred_high_copy=deferred if deferred_info is None else deferred_info['pairs_pending'],
                deferred=deferred_info, gens=a.gens, source_alignment=dict(alignment),
                ucd_input=str(out / 'ucd_input.tsv.gz'),
                ucd_input_sha256=ub.sha_file(out / 'ucd_input.tsv.gz'),
                dna_pairs_all_sha256=ub.sha_file(out / 'dna_pairs_all.tsv.gz'))
    (out / 'collect.json').write_text(json.dumps(info, indent=1) + '\n')
    print(json.dumps({k: v for k, v in info.items() if k != 'params'}, indent=1))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = p.add_subparsers(dest='cmd')
    sp.required = True
    for name in ('plan', 'seeds', 'library', 'occurrences', 'pairs', 'enumerate', 'align', 'collect'):
        x = sp.add_parser(name)
        x.add_argument('--work', required=True)
        x.add_argument('--manifest', required=True)
        x.add_argument('--threads', type=int, default=int(os.environ.get('SLURM_CPUS_PER_TASK', 4)))
        x.add_argument('--index', type=int, default=int(os.environ.get('SLURM_ARRAY_TASK_ID', 0)))
        if name in ('plan', 'seeds', 'library'):
            x.add_argument('--tag', required=True)
        if name == 'plan':
            x.add_argument('--batches', required=True)
            x.add_argument('--chunk-size', type=int, default=100)
            x.add_argument('--part-size', type=int, default=100)
        if name == 'occurrences':
            x.add_argument('--list', action='store_true', help='print the array indices of tasks not done')
        if name == 'pairs':
            x.add_argument('--force', action='store_true')
        if name in ('pairs', 'collect'):
            x.add_argument('--subs', type=int, default=1, help='pair tasks per bucket (array size buckets x subs)')
        if name in ('enumerate', 'align'):
            x.add_argument('--gen', type=int, default=0, help='generation (0 = first; later = resumed deferred chunks / '
                                                           'new batches)')
        if name == 'enumerate':
            x.add_argument('--max-pairs', type=int, default=0, help='pairs enumerated by this task (0 = all)')
        if name == 'collect':
            x.add_argument('--out', required=True)
            x.add_argument('--gens', default='', help='v0.2: generations of align shards to collect (e.g. 0,1)')
    a = p.parse_args()
    dict(plan=cmd_plan, seeds=cmd_seeds, library=cmd_library, occurrences=cmd_occurrences, pairs=cmd_pairs,
         enumerate=cmd_enumerate, align=cmd_align, collect=cmd_collect)[a.cmd](a)


if __name__ == '__main__':
    main()
