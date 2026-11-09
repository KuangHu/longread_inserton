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

VERSION = 'dna_discovery_v0.1.2'
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
    full_pairs=20000,             # all pairs of a library sequence up to this many
    design_partners=8,            # partners per occurrence above it
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
    pieces = cr.chain([(h['qs'], h['qe'], h['ss'], h['se']) for h in hsps], a0, a1, b0, b1, p['chain_gap'])
    if not pieces:
        return dict(status='no_shared_block_at_seed')
    e = cr.extensions(pieces, a0, a1, len(wa), b0, b1, len(wb), p['edge'])
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


# ---------------------------------------------------------------- collect

def cmd_collect(a):
    """ucd_batch input table of the resolved pairs (candidate_insertion), one row per distinct pair of X intervals."""
    w = Path(a.work)
    S = PARAMS['buckets'] * a.subs
    miss = [i for i in range(S) if not (w / 'pairs' / ('shard_%04d' % i) / 'done.json').exists()]
    if miss:
        sys.exit('pair shards not done: %d (first %s)' % (len(miss), miss[:5]))
    man = {m['sample']: m for m in rows(a.manifest)}
    seen, dup, st, scope, n_all = {}, Counter(), Counter(), Counter(), 0
    deferred = 0
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    ci = {c: i for i, c in enumerate(ub.INPUT_COLS)}
    with gzip.open(out / 'dna_pairs_all.tsv.gz', 'wt') as fa:
        fa.write('\t'.join(PAIR_COLS) + '\n')
        for i in range(S):
            d = w / 'pairs' / ('shard_%04d' % i)
            deferred += json.loads((d / 'done.json').read_text())['pairs_deferred']
            for r in rows(d / 'pairs.tsv.gz'):
                n_all += 1
                fa.write('\t'.join(str(r.get(c, '')) for c in PAIR_COLS) + '\n')
                st[r['status']] += 1
                if r['status'] != 'resolved':
                    continue
                if r['pair_id'] in seen:
                    dup[r['pair_id']] += 1
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
                resolved_distinct=len(seen), resolved_duplicate_rows=sum(dup.values()), scope=dict(scope),
                pairs_deferred_high_copy=deferred, ucd_input=str(out / 'ucd_input.tsv.gz'),
                ucd_input_sha256=ub.sha_file(out / 'ucd_input.tsv.gz'),
                dna_pairs_all_sha256=ub.sha_file(out / 'dna_pairs_all.tsv.gz'))
    (out / 'collect.json').write_text(json.dumps(info, indent=1) + '\n')
    print(json.dumps({k: v for k, v in info.items() if k != 'params'}, indent=1))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = p.add_subparsers(dest='cmd')
    sp.required = True
    for name in ('plan', 'seeds', 'library', 'occurrences', 'pairs', 'collect'):
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
        if name == 'collect':
            x.add_argument('--out', required=True)
    a = p.parse_args()
    dict(plan=cmd_plan, seeds=cmd_seeds, library=cmd_library, occurrences=cmd_occurrences, pairs=cmd_pairs,
         collect=cmd_collect)[a.cmd](a)


if __name__ == '__main__':
    main()
