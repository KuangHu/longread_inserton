"""U-insertion-D comparison: is the same insertion found between different contexts?

For a pair of occurrences

    position A:  U_A - X  - D_A
    position B:  U_B - X' - D_B

the module compares DNA directly and decides whether X and X' are similar while U and D both differ. It requires no
empty site, says nothing about mobility, and does not use gap callers, S6 site groups or read-check results as
admission criteria. Read-support and other labels of earlier steps are carried as columns for viewing only.

Commands
  prepare   build the input pair table (one row per occurrence pair, both sides' coordinates, interval of each kind,
            labels) from existing assembly / candidate tables; no annotation, download or read search.
  run       compare every pair for each interval kind; writes ucd_pairs.tsv, counts, and the sequence / HSP evidence.
  figures   draw U-insertion-D figures of a few representative pairs, with the measured identity and coverage.

Interval kinds (kept apart in every row and every count; a CDS is never called a complete insertion):
  cds_core              the core coding span (catalog CDS) of each occurrence
  candidate_insertion   coding span plus the pair's own shared extensions with both boundaries resolved (the
                        boundary review's element interval); pair-specific, and only where it exists

Sequences (all in the insertion's orientation): X = interval; U = the CONTEXT bp before X; D = the CONTEXT bp after X.
A minus-strand occurrence is reverse-complemented, so U is always 5' of X as read along X.

Comparison (one blastn per pair and interval kind; thresholds in CONFIG, written to every output)
  1. X_A vs X_B, both strands. Reliable HSPs: identity >= ins_identity and length >= min_hsp. Coverage of each
     sequence = union of reliable HSPs on the better strand (non-redundant), identity = alignment-length weighted.
     insertion_similar: both coverages >= ins_coverage. Some reliable HSP but less coverage: insertion_partial.
     None: insertion_dissimilar. If X' matches X on the minus strand, B is flipped (U_B' = rc(D_B), D_B' = rc(U_B))
     so that both contexts are read in the same insertion orientation.
  2. Contexts, after that orientation: U_A/U_B and D_A/D_B on the plus strand (same side), and as checks for
     swapped or rearranged context U_A/U_B and D_A/D_B on the minus strand and U_A/D_B, D_A/U_B on both strands.
     Each comparison: union of reliable HSPs (identity >= ctx_identity, length >= min_hsp) on each flank.
       similar    both flanks covered >= ctx_similar_cov
       different  both flanks covered <= ctx_different_cov
       uncertain  anything between (not forced into either)
  3. Context usable: a flank with fewer than min_usable_bp bases that are not N and not in low-complexity windows
     (3-mer entropy < lc_entropy bits in lc_window bp windows) is insufficient (contig edge, N runs, repeats of
     low complexity). Missing context is never counted as different.
  4. Same assembly: if a flank of one occurrence overlaps a flank or the insertion of the other on the same contig
     (the two reuse one physical stretch), the pair is shared_physical_context whatever the sequences say.

Verdict (first that applies; every applicable reason is also listed):
  insertion_dissimilar / insertion_partial
  shared_physical_context
  insufficient_context
  swapped_or_rearranged_context    an inverted same-side or a cross-side comparison is similar
  same_context                     U similar and D similar
  only_one_side_different          one side similar, the other different
  ucd_supported                    U different, D different, and every inverted / cross comparison different
  uncertain_context                anything else (a side or a cross check uncertain)

ucd_supported is a structural pattern only (similar insertion, different contexts). Conserved paralogs that meet the
rule stay in, with their labels; nothing is excluded by product name.

Pair order is canonical (A = the smaller occurrence id), so swapping the input sides or shuffling rows does not
change a result; reverse-complementing a side's sequences does not change the structural verdict (tests).
"""
import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import random
import shutil
import subprocess
import sys
import tempfile
import time
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

VERSION = 'ucd_compare_v0.1.1'  # v0.1.1: rc() complements every IUPAC code (v0.1 left R/Y/K/M/... uncomplemented)
HERE = Path(__file__).resolve().parent
REPO = HERE / 'scale_comamonadaceae_v0.1'
SCRATCH = Path('/global/scratch/users/kh36969/longread_callout/core_context')
CONFIG = dict(
    context=1000,             # bp of U and of D
    ins_identity=90.0, ins_coverage=0.90,
    ctx_identity=80.0, ctx_similar_cov=0.80, ctx_different_cov=0.20,
    min_hsp=30,               # bp; shorter HSPs are not reliable
    evalue=1e-5,
    min_usable_bp=800,        # of CONTEXT; fewer usable bases -> insufficient_context
    lc_window=64, lc_step=32, lc_entropy=3.0,
)
KINDS = ('cds_core', 'candidate_insertion')
_FWD = 'ACGTURYSWKMBDHVN'
_REV = 'TGCAAYRSWMKVHDBN'
COMP = str.maketrans(_FWD + _FWD.lower(), _REV + _REV.lower())  # full IUPAC, case kept
SIDE = ('occurrence', 'sample', 'contig', 'start0', 'end0', 'strand', 'elem_lo', 'elem_hi', 'boundary_up',
        'boundary_down', 'estimated_span', 'position_class')
LABELS = ('pair_id', 'source', 'case', 'tier', 'structure', 'core_key', 'product', 'product_hint', 'control_class',
          'pair_scope', 'core_match_class', 'core_dna_identity', 'pair_class', 'reference_context_class')
OUTFMT = '6 qseqid sseqid pident length nident qstart qend sstart send qlen slen evalue bitscore sstrand'


# ---------------------------------------------------------------- io

def read_tsv(path):
    op = gzip.open if str(path).endswith('.gz') else open
    with op(path, 'rt') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def write_tsv(path, rows, cols):
    op = gzip.open if str(path).endswith('.gz') else open
    with op(path, 'wt') as f:
        w = csv.DictWriter(f, cols, delimiter='\t', lineterminator='\n', extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)


def read_fasta(path):
    op = gzip.open if str(path).endswith('.gz') else open
    seqs, name, buf = {}, None, []
    with op(path, 'rt') as f:
        for line in f:
            if line.startswith('>'):
                if name is not None:
                    seqs[name] = ''.join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name is not None:
        seqs[name] = ''.join(buf)
    return seqs


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def rc(s):
    return s.translate(COMP)[::-1]


# ---------------------------------------------------------------- sequences

def sides_of(contig_seq, lo, hi, strand, L):
    """X, U, D in the insertion's orientation and their genomic intervals (contig coordinates, no wrap-around)."""
    n = len(contig_seq)
    left, right = (max(0, lo - L), lo), (hi, min(n, hi + L))
    x = contig_seq[lo:hi].upper()
    l_, r_ = contig_seq[left[0]:left[1]].upper(), contig_seq[right[0]:right[1]].upper()
    if strand == '+':
        return dict(X=x, U=l_, D=r_, X_iv=(lo, hi), U_iv=left, D_iv=right)
    return dict(X=rc(x), U=rc(r_), D=rc(l_), X_iv=(lo, hi), U_iv=right, D_iv=left)


def usable_bp(s, cfg):
    """Bases that are not N and not inside a low-complexity window (3-mer entropy < lc_entropy bits)."""
    bad = [c not in 'ACGT' for c in s]
    w, st = cfg['lc_window'], cfg['lc_step']
    for i in range(0, max(1, len(s) - w + 1), st):
        win = s[i:i + w]
        k = Counter(win[j:j + 3] for j in range(len(win) - 2) if 'N' not in win[j:j + 3])
        tot = sum(k.values())
        ent = -sum(v / tot * math.log2(v / tot) for v in k.values()) if tot else 0.0
        if ent < cfg['lc_entropy']:
            for j in range(i, min(len(s), i + w)):
                bad[j] = True
    return len(s) - sum(bad)


# ---------------------------------------------------------------- alignment

def blast_pair(seqs_a, seqs_b, tmp, cfg):
    """All HSPs of the sequences of A (query) against those of B (subject), both strands."""
    qa, sb = Path(tmp) / 'a.fa', Path(tmp) / 'b.fa'
    qa.write_text(''.join('>%s\n%s\n' % (k, v) for k, v in seqs_a.items() if v))
    sb.write_text(''.join('>%s\n%s\n' % (k, v) for k, v in seqs_b.items() if v))
    if not qa.stat().st_size or not sb.stat().st_size:
        return []
    out = subprocess.run(['blastn', '-task', 'blastn', '-query', str(qa), '-subject', str(sb), '-evalue',
                          str(cfg['evalue']), '-dust', 'no', '-outfmt', OUTFMT], stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE, universal_newlines=True, check=True).stdout
    hsps = []
    for line in out.splitlines():
        f = line.split('\t')
        q, s, pid, ln, nid, qs, qe, ss, se, ql, sl, ev, bs, sst = f
        ss, se = int(ss), int(se)
        hsps.append(dict(q=q, s=s, pident=float(pid), length=int(ln), nident=int(nid), qs=int(qs) - 1, qe=int(qe),
                         ss=min(ss, se) - 1, se=max(ss, se), qlen=int(ql), slen=int(sl), evalue=float(ev),
                         bitscore=float(bs), strand='plus' if sst == 'plus' else 'minus'))
    return hsps


def union_len(ivs):
    tot, cur = 0, None
    for a, b in sorted(ivs):
        if cur is None or a > cur[1]:
            if cur:
                tot += cur[1] - cur[0]
            cur = [a, b]
        else:
            cur[1] = max(cur[1], b)
    return tot + (cur[1] - cur[0] if cur else 0)


def measure(hsps, q, s, strand, min_id, cfg, qlen, slen):
    """Non-redundant coverage of both sequences by the reliable HSPs of one (q, s, strand) comparison."""
    h = [x for x in hsps if x['q'] == q and x['s'] == s and x['strand'] == strand
         and x['pident'] >= min_id and x['length'] >= cfg['min_hsp']]
    if not h or not qlen or not slen:
        return dict(identity=0.0, cov_a=0.0, cov_b=0.0, hsps=len(h))
    return dict(identity=round(100.0 * sum(x['nident'] for x in h) / sum(x['length'] for x in h), 2),
                cov_a=round(union_len([(x['qs'], x['qe']) for x in h]) / qlen, 3),
                cov_b=round(union_len([(x['ss'], x['se']) for x in h]) / slen, 3), hsps=len(h))


def ctx_class(m, cfg):
    if m['cov_a'] >= cfg['ctx_similar_cov'] and m['cov_b'] >= cfg['ctx_similar_cov']:
        return 'similar'
    if m['cov_a'] <= cfg['ctx_different_cov'] and m['cov_b'] <= cfg['ctx_different_cov']:
        return 'different'
    return 'uncertain'


def overlaps(a, b):
    return a[0] < b[1] and b[0] < a[1]


def physical_overlap(A, B):
    """Same contig of the same assembly, and a flank of one overlaps a flank or the insertion of the other."""
    if not (A.get('sample') and A.get('sample') == B.get('sample') and A.get('contig') == B.get('contig')):
        return []
    out = []
    for p, q in (('U', 'U'), ('U', 'D'), ('D', 'U'), ('D', 'D'), ('U', 'X'), ('D', 'X'), ('X', 'U'), ('X', 'D'),
                 ('X', 'X')):
        if overlaps(A[p + '_iv'], B[q + '_iv']):
            out.append('%s_A~%s_B' % (p, q))
    return out


def compare(A, B, cfg=CONFIG, tmp=None, keep=None):
    """A, B: dict(X, U, D, and optionally sample, contig, X_iv, U_iv, D_iv, id). Returns the result row (and fills
    `keep` with the oriented sequences and HSPs as evidence)."""
    if (A.get('id') or '') > (B.get('id') or ''):
        A, B = B, A
    own = tmp is None
    tmp = tmp or tempfile.mkdtemp(prefix='ucd_')
    try:
        h = blast_pair({'X': A['X']}, {'X': B['X']}, tmp, cfg)
        ins = {st: measure(h, 'X', 'X', st, cfg['ins_identity'], cfg, len(A['X']), len(B['X']))
               for st in ('plus', 'minus')}
        ori = max(('plus', 'minus'), key=lambda st: (min(ins[st]['cov_a'], ins[st]['cov_b']), st == 'plus'))
        mi = ins[ori]
        if mi['cov_a'] >= cfg['ins_coverage'] and mi['cov_b'] >= cfg['ins_coverage']:
            ins_cls = 'insertion_similar'
        elif mi['hsps']:
            ins_cls = 'insertion_partial'
        else:
            ins_cls = 'insertion_dissimilar'
        Bo = dict(B) if ori == 'plus' else dict(B, X=rc(B['X']), U=rc(B['D']), D=rc(B['U']),
                                                  U_iv=B.get('D_iv'), D_iv=B.get('U_iv'))
        hc = blast_pair({'U': A['U'], 'D': A['D']}, {'U': Bo['U'], 'D': Bo['D']}, tmp, cfg)
        L = {('A', k): len(A[k]) for k in 'UD'}
        L.update({('B', k): len(Bo[k]) for k in 'UD'})
        cmp = {}
        for qa, sb, st, name in (('U', 'U', 'plus', 'UU'), ('D', 'D', 'plus', 'DD'), ('U', 'U', 'minus', 'UU_inv'),
                                 ('D', 'D', 'minus', 'DD_inv'), ('U', 'D', 'plus', 'UD'), ('U', 'D', 'minus', 'UD_inv'),
                                 ('D', 'U', 'plus', 'DU'), ('D', 'U', 'minus', 'DU_inv')):
            m = measure(hc, qa, sb, st, cfg['ctx_identity'], cfg, L[('A', qa)], L[('B', sb)])
            m['class'] = ctx_class(m, cfg)
            cmp[name] = m
        usable = {('A', k): usable_bp(A[k], cfg) for k in 'UD'}
        usable.update({('B', k): usable_bp(Bo[k], cfg) for k in 'UD'})
        short = ['%s_%s' % k for k, v in sorted(usable.items()) if v < cfg['min_usable_bp']]
        phys = physical_overlap(A, Bo) if A.get('U_iv') and Bo.get('U_iv') else []
        cross = [k for k in ('UU_inv', 'DD_inv', 'UD', 'UD_inv', 'DU', 'DU_inv') if cmp[k]['class'] == 'similar']
        cross_unc = [k for k in ('UU_inv', 'DD_inv', 'UD', 'UD_inv', 'DU', 'DU_inv') if cmp[k]['class'] == 'uncertain']
        u, d = cmp['UU']['class'], cmp['DD']['class']
        reasons = []
        if ins_cls != 'insertion_similar':
            reasons.append(ins_cls)
        if phys:
            reasons.append('shared_physical_context:' + ','.join(phys))
        if short:
            reasons.append('insufficient_context:' + ','.join(short))
        if cross:
            reasons.append('cross_or_inverted_similar:' + ','.join(cross))
        if cross_unc:
            reasons.append('cross_or_inverted_uncertain:' + ','.join(cross_unc))
        reasons.append('upstream_%s' % u)
        reasons.append('downstream_%s' % d)
        if ins_cls != 'insertion_similar':
            verdict = ins_cls
        elif phys:
            verdict = 'shared_physical_context'
        elif short:
            verdict = 'insufficient_context'
        elif cross:
            verdict = 'swapped_or_rearranged_context'
        elif u == d == 'similar':
            verdict = 'same_context'
        elif {u, d} == {'similar', 'different'}:
            verdict = 'only_one_side_different'
        elif u == d == 'different' and not cross_unc:
            verdict = 'ucd_supported'
        else:
            verdict = 'uncertain_context'
        row = dict(occurrence_a=A.get('id', ''), occurrence_b=B.get('id', ''), b_flipped=ori == 'minus',
                   x_len_a=len(A['X']), x_len_b=len(B['X']), insertion_class=ins_cls,
                   insertion_identity=mi['identity'], insertion_cov_a=mi['cov_a'], insertion_cov_b=mi['cov_b'],
                   insertion_orientation=ori,
                   **{'%s_len_%s' % (k, s): L[(s.upper(), k)] for k in 'UD' for s in 'ab'},
                   **{'%s_usable_%s' % (k, s): usable[(s.upper(), k)] for k in 'UD' for s in 'ab'},
                   **{'%s_%s' % (c, f): cmp[c][f] for c in cmp for f in ('class', 'identity', 'cov_a', 'cov_b')},
                   upstream=u, downstream=d, verdict=verdict, reasons=';'.join(reasons))
        if keep is not None:
            keep.update(seqs={'A_U': A['U'], 'A_X': A['X'], 'A_D': A['D'], 'B_U': Bo['U'], 'B_X': Bo['X'],
                              'B_D': Bo['D']}, hsps=[dict(x, part='insertion') for x in h] +
                        [dict(x, part='context') for x in hc])
        return row
    finally:
        if own:
            shutil.rmtree(tmp, ignore_errors=True)


# ---------------------------------------------------------------- prepare

def manifest(path):
    return {r['sample']: r for r in read_tsv(path)}


def cmd_prepare(a):
    """Input pairs: (1) every read-check round-5 pair whose two samples are both in --samples (dev / validation);
    (2) the named regression cases (--cases), built from the boundary-review pair shards and occurrence shards (the
    same element() rule as ucd_readcheck.pairs); labels of earlier steps are copied for viewing only."""
    samples = set()
    for p in a.samples.split(','):
        samples |= {r['sample'] for r in read_tsv(p)}
    r5 = read_tsv(a.round5)
    ref = {}
    if a.reference and Path(a.reference).exists():
        for r in read_tsv(a.reference):
            ref[frozenset((r['occurrence_a'], r['occurrence_b']))] = r['context_class']
    rows = []
    for r in r5:
        if r['a_sample'] in samples and r['b_sample'] in samples:
            x = {k: r.get(k, '') for k in LABELS}
            x.update(source='dev_validation_samples', case='')
            x['reference_context_class'] = ref.get(frozenset((r['a_occurrence'], r['b_occurrence'])), '')
            for s in 'ab':
                x.update({'%s_%s' % (s, f): r.get('%s_%s' % (s, f), '') for f in SIDE})
            rows.append(x)
    have = {frozenset((r['a_occurrence'], r['b_occurrence'])): r for r in rows}
    r5ix = {frozenset((r['a_occurrence'], r['b_occurrence'])): r for r in r5}
    cases = read_tsv(a.cases) if a.cases else []
    need = {frozenset((c['occurrence_a'], c['occurrence_b'])): c for c in cases}
    found, occ_need = {}, set()
    work = Path(a.work)
    for k in range(64):
        p = work / 'review' / 'boundary' / a.tag / ('shard_%03d' % k) / 'pairs.tsv.gz'
        for r in read_tsv(p):
            key = frozenset((r['occurrence_a'], r['occurrence_b']))
            c = need.get(key)
            if c and r['core_key'] in (c['new_core_key_a'], c['new_core_key_b']) and key not in found:
                found[key] = r
                occ_need |= key
    occ = {}
    for k in range(64):
        for r in read_tsv(work / 'review' / 's6' / a.tag / ('shard_%03d' % k) / 'occurrences.tsv.gz'):
            if r['occurrence_id'] in occ_need:
                occ[r['occurrence_id']] = r
    bnd = {}
    for scope in ('same_assembly', 'cross_assembly'):
        p = REPO / ('boundary_occurrences_%s_%s.tsv.gz' % (scope, a.tag))
        for r in read_tsv(p):
            if r['occurrence_id'] in occ_need:
                bnd[(scope, r['structure'], r['occurrence_id'])] = r
    missing = []
    for key, c in need.items():
        if key in have:
            have[key]['case'] = c['candidate']
            continue
        r = found.get(key)
        if r is None:
            missing.append(c['candidate'])
            continue
        x = dict(source='named_case', case=c['candidate'], structure=r['structure'], core_key=r['core_key'],
                 pair_scope=r['pair_scope'], core_match_class=r['core_match_class'],
                 core_dna_identity=r['core_dna_identity'], tier='', product='', product_hint='', control_class='')
        old = r5ix.get(key)
        x['pair_id'] = old['pair_id'] if old else 'NC_' + c['candidate']
        x['pair_class'] = old['pair_class'] if old else 'not_in_read_check'
        x['reference_context_class'] = ref.get(key, '')
        for s in 'ab':
            o = occ[r['occurrence_' + s]]
            st, en = int(o['start0']), int(o['end0'])
            if r[s + '_up'] == '' or r[s + '_down'] == '':  # extension not determined: no candidate interval
                lo = hi = ''
            else:
                up, down = int(r[s + '_up']), int(r[s + '_down'])
                lo, hi = (st - up, en + down) if o['strand'] == '+' else (st - down, en + up)
            b = bnd.get((r['pair_scope'], r['structure'], o['occurrence_id']), {})
            x.update({'%s_%s' % (s, f): v for f, v in dict(
                occurrence=o['occurrence_id'], sample=o['sample'], contig=o['contig'], start0=st, end0=en,
                strand=o['strand'], elem_lo=lo, elem_hi=hi, boundary_up=b.get('up_status', ''),
                boundary_down=b.get('down_status', ''), estimated_span=b.get('estimated_span', ''),
                position_class=(old or {}).get(s + '_position_class', 'not_in_read_check')).items()})
        rows.append(x)
    cols = list(LABELS) + ['%s_%s' % (s, f) for s in 'ab' for f in SIDE]
    rows.sort(key=lambda r: (r['source'], r['pair_id']))
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    write_tsv(out / 'input_pairs.tsv.gz', rows, cols)
    info = dict(version=VERSION, pairs=len(rows), by_source=dict(Counter(r['source'] for r in rows)),
                named_cases_found=sorted(r['case'] for r in rows if r['case']), named_cases_missing=sorted(missing),
                samples=len(samples), inputs={str(p): sha256(p) for p in [Path(a.round5)] + [Path(x) for x in
                                                                                       a.samples.split(',')]})
    (out / 'prepare.json').write_text(json.dumps(info, indent=1) + '\n')
    print(json.dumps(info, indent=1))


# ---------------------------------------------------------------- run

_G = {}


def _init(man, cfg):
    _G.update(man=man, cfg=cfg, genomes={})


def genome(sample):
    g = _G['genomes']
    if sample not in g:
        if len(g) > 8:
            g.pop(next(iter(g)))
        g[sample] = read_fasta(_G['man'][sample]['fasta'])
    return g[sample]


def interval(r, s, kind):
    if kind == 'cds_core':
        return int(r[s + '_start0']), int(r[s + '_end0'])
    if r.get(s + '_elem_lo', '') == '' or r.get(s + '_boundary_up') != 'resolved' or \
            r.get(s + '_boundary_down') != 'resolved':
        return None
    return int(r[s + '_elem_lo']), int(r[s + '_elem_hi'])


def side(r, s, kind):
    iv = interval(r, s, kind)
    if iv is None:
        return None
    seq = genome(r[s + '_sample'])[r[s + '_contig']]
    d = sides_of(seq, iv[0], iv[1], r[s + '_strand'], _G['cfg']['context'])
    d.update(id=r[s + '_occurrence'], sample=r[s + '_sample'], contig=r[s + '_contig'], contig_len=len(seq))
    return d


def work_one(r):
    out, ev = [], []
    t0 = time.time()
    tmp = tempfile.mkdtemp(prefix='ucd_', dir=os.environ.get('TMPDIR'))
    try:
        for kind in KINDS:
            A, B = side(r, 'a', kind), side(r, 'b', kind)
            base = {k: r.get(k, '') for k in LABELS}
            if A is None or B is None:
                out.append(dict(base, interval_kind=kind, pair_key='%s|%s' % (r['pair_id'], kind),
                                occurrence_a=r['a_occurrence'], occurrence_b=r['b_occurrence'],
                                a_sample=r['a_sample'], b_sample=r['b_sample'], verdict='interval_not_available'))
                continue
            keep = {}
            x = compare(A, B, _G['cfg'], tmp, keep)
            sw = x['occurrence_a'] != r['a_occurrence']
            P, Q = ('b', 'a') if sw else ('a', 'b')
            loc = {}
            for s, src, D in (('a', P, A if not sw else B), ('b', Q, B if not sw else A)):
                loc.update({'%s_%s' % (s, f): r.get('%s_%s' % (src, f), '') for f in
                            ('sample', 'contig', 'strand', 'boundary_up', 'boundary_down', 'estimated_span',
                             'position_class')})
                loc['%s_x_interval' % s] = '%d-%d' % D['X_iv']
                loc['%s_u_interval' % s] = '%d-%d' % D['U_iv']
                loc['%s_d_interval' % s] = '%d-%d' % D['D_iv']
            key = '%s|%s' % (r['pair_id'], kind)
            out.append(dict(base, interval_kind=kind, pair_key=key, **loc, **x))
            ev.append(dict(key=key, seqs=keep['seqs'], hsps=keep['hsps']))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    return out, ev, time.time() - t0


RESULT_COLS = (['pair_key', 'interval_kind'] + list(LABELS) +
               ['occurrence_a', 'a_sample', 'a_contig', 'a_strand', 'a_x_interval', 'a_u_interval', 'a_d_interval',
                'occurrence_b', 'b_sample', 'b_contig', 'b_strand', 'b_x_interval', 'b_u_interval', 'b_d_interval',
                'b_flipped', 'x_len_a', 'x_len_b', 'insertion_class', 'insertion_identity', 'insertion_cov_a',
                'insertion_cov_b', 'insertion_orientation', 'U_len_a', 'U_len_b', 'D_len_a', 'D_len_b', 'U_usable_a',
                'U_usable_b', 'D_usable_a', 'D_usable_b'] +
               ['%s_%s' % (c, f) for c in ('UU', 'DD', 'UU_inv', 'DD_inv', 'UD', 'UD_inv', 'DU', 'DU_inv')
                for f in ('class', 'identity', 'cov_a', 'cov_b')] +
               ['upstream', 'downstream', 'verdict', 'reasons', 'a_boundary_up', 'a_boundary_down',
                'a_estimated_span', 'b_boundary_up', 'b_boundary_down', 'b_estimated_span', 'a_position_class',
                'b_position_class'])


def counts(rows):
    """Per interval kind and verdict: pairs, distinct occurrences, samples, cores (pairs are not events)."""
    out = {}
    for kind in KINDS:
        rk = [r for r in rows if r['interval_kind'] == kind]
        by = defaultdict(list)
        for r in rk:
            by[r['verdict']].append(r)
        out[kind] = {v: dict(pairs=len(x), occurrences=len({o for r in x for o in (r.get('occurrence_a') or
                                                                                    r.get('a_occurrence', ''),
                                                                                    r.get('occurrence_b') or
                                                                                    r.get('b_occurrence', ''))}),
                             samples=len({r.get('a_sample', '') for r in x} | {r.get('b_sample', '') for r in x}),
                             cores=len({(r['structure'], r['core_key']) for r in x}))
                     for v, x in sorted(by.items())}
    return out


def cmd_run(a):
    cfg = dict(CONFIG)
    if a.config:
        cfg.update(json.loads(Path(a.config).read_text()))
    out, work = Path(a.out), Path(a.work)
    work.mkdir(parents=True, exist_ok=True)
    pairs = read_tsv(out / 'input_pairs.tsv.gz')
    man = manifest(a.manifest)
    miss = sorted({r[s + '_sample'] for r in pairs for s in 'ab'} - set(man))
    if miss:
        sys.exit('samples without an assembly in the manifest: %s' % ','.join(miss))
    pairs.sort(key=lambda r: (r['a_sample'], r['b_sample'], r['pair_id']))  # genome cache locality
    t0 = time.time()
    rows, secs = [], []
    seq_f = gzip.open(work / 'evidence_sequences.fa.gz', 'wt')
    hsp_f = gzip.open(work / 'evidence_hsps.tsv.gz', 'wt')
    hcols = ['pair_key', 'part', 'q', 's', 'strand', 'pident', 'length', 'nident', 'qs', 'qe', 'ss', 'se', 'qlen',
             'slen', 'evalue', 'bitscore']
    hw = csv.DictWriter(hsp_f, hcols, delimiter='\t', lineterminator='\n', extrasaction='ignore')
    hw.writeheader()
    with ProcessPoolExecutor(a.threads, initializer=_init, initargs=(man, cfg)) as ex:
        for res, ev, dt in ex.map(work_one, pairs, chunksize=4):
            rows += res
            secs.append(dt)
            for e in ev:
                for k, v in e['seqs'].items():
                    seq_f.write('>%s|%s\n%s\n' % (e['key'], k, v))
                for h in e['hsps']:
                    hw.writerow(dict(h, pair_key=e['key']))
    seq_f.close()
    hsp_f.close()
    wall = time.time() - t0
    rows.sort(key=lambda r: (r['interval_kind'], r['source'], r.get('pair_key') or r['pair_id']))
    write_tsv(out / 'ucd_pairs.tsv.gz', rows, RESULT_COLS)
    summ = dict(
        version=VERSION, config=cfg, input_pairs=len(pairs), comparisons=sum(r['verdict'] != 'interval_not_available'
                                                                             for r in rows),
        by_kind_and_verdict=counts(rows),
        named_cases={r['case']: {k: x['verdict'] for k in KINDS for x in rows if x['case'] == r['case'] and
                                 x['interval_kind'] == k} for r in pairs if r['case']},
        verdict_by_reference_context_class={k: dict(sorted(Counter('%s | %s' % (r['reference_context_class'] or '-',
                                                                                  r['verdict'])
                                                                    for r in rows if r['interval_kind'] == k).items()))
                                            for k in KINDS},
        timing=dict(wall_s=round(wall, 1), threads=a.threads, cpu_s_per_pair=round(sum(secs) / max(1, len(secs)), 3)),
        evidence={p.name: dict(path=str(p), sha256=sha256(p), bytes=p.stat().st_size)
                  for p in (work / 'evidence_sequences.fa.gz', work / 'evidence_hsps.tsv.gz')},
        input_pairs_sha256=sha256(out / 'input_pairs.tsv.gz'))
    (out / 'summary.json').write_text(json.dumps(summ, indent=1) + '\n')
    print(json.dumps({k: v for k, v in summ.items() if k not in ('named_cases', 'verdict_by_reference_context_class')},
                     indent=1))


# ---------------------------------------------------------------- figures

COLORS = dict(U='#6b8fb8', X='#d9822b', D='#5aa36b', cov='#1f2d3d', cross='#c0392b', text='#222', grid='#999')


def svg_pair(r, hsps, title):
    """Two tracks (A over B), each U | X | D in insertion orientation, with the reliable HSP blocks of the same-side
    and the best cross comparison drawn on the flanks and the measured identity / coverage written next to them."""
    W, H, M = 980, 320, 40
    lens = [int(r['U_len_a']), int(r['x_len_a']), int(r['D_len_a']), int(r['U_len_b']), int(r['x_len_b']),
            int(r['D_len_b'])]
    total = max(lens[0] + lens[1] + lens[2], lens[3] + lens[4] + lens[5])
    sc = (W - 2 * M) / total
    ys = {'A': 70, 'B': 170}
    o = ['<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" font-family="sans-serif" font-size="11">'
         % (W, H), '<rect width="100%" height="100%" fill="white"/>',
         '<text x="%d" y="18" font-size="13" fill="%s">%s</text>' % (M, COLORS['text'], title)]
    x0 = {}
    for t, (lu, lx, ld) in (('A', lens[:3]), ('B', lens[3:])):
        y, x = ys[t], M
        for part, ln in (('U', lu), ('X', lx), ('D', ld)):
            o.append('<rect x="%.1f" y="%d" width="%.1f" height="18" fill="%s" opacity="0.55"/>'
                     % (x, y, ln * sc, COLORS[part]))
            o.append('<text x="%.1f" y="%d" fill="%s">%s_%s %d bp</text>' % (x + 2, y - 4, COLORS['text'], part, t, ln))
            x0[(t, part)] = x
            x += ln * sc
        lab = r['occurrence_a'] if t == 'A' else r['occurrence_b']
        o.append('<text x="%d" y="%d" fill="%s">%s %s</text>' % (M, y + 34, COLORS['grid'], lab,
                                                                  '(flipped)' if t == 'B' and r['b_flipped'] in
                                                                  (True, 'True') else ''))
    for h in hsps:
        if h['part'] != 'context' or float(h['pident']) < CONFIG['ctx_identity'] or int(h['length']) < CONFIG['min_hsp']:
            continue
        same = h['q'] == h['s'] and h['strand'] == 'plus'
        for t, part, a_, b_ in (('A', h['q'], int(h['qs']), int(h['qe'])), ('B', h['s'], int(h['ss']), int(h['se']))):
            o.append('<rect x="%.1f" y="%d" width="%.1f" height="6" fill="%s"/>'
                     % (x0[(t, part)] + a_ * sc, ys[t] + (19 if same else 26), max(1, (b_ - a_) * sc),
                        COLORS['cov'] if same else COLORS['cross']))
    txt = ['X/X′: identity %s%%, coverage %s / %s (%s)' % (r['insertion_identity'], r['insertion_cov_a'],
                                                                r['insertion_cov_b'], r['insertion_class']),
           'U/U: %s (id %s%%, cov %s / %s)   D/D: %s (id %s%%, cov %s / %s)' % (
               r['UU_class'], r['UU_identity'], r['UU_cov_a'], r['UU_cov_b'], r['DD_class'], r['DD_identity'],
               r['DD_cov_a'], r['DD_cov_b']),
           'inverted same side: ' + ', '.join('%s %s (cov %s / %s)' % (c, r[c + '_class'], r[c + '_cov_a'],
                                                                       r[c + '_cov_b']) for c in ('UU_inv', 'DD_inv')),
           'cross side: ' + ', '.join('%s %s (cov %s / %s)' % (c, r[c + '_class'], r[c + '_cov_a'], r[c + '_cov_b'])
                                     for c in ('UD', 'UD_inv', 'DU', 'DU_inv')),
           'verdict: %s   [%s]' % (r['verdict'], r['interval_kind'])]
    for i, t in enumerate(txt):
        o.append('<text x="%d" y="%d" fill="%s">%s</text>' % (M, 222 + i * 15, COLORS['text'],
                                                              t.replace('&', '&amp;')))
    o.append('<text x="%d" y="%d" fill="%s">blocks under a flank: reliable HSPs (≥%d%%, ≥%d bp); dark = same side, '
             'same strand; red = inverted or cross side</text>' % (M, H - 6, COLORS['grid'], CONFIG['ctx_identity'], CONFIG['min_hsp']))
    o.append('</svg>')
    return '\n'.join(o)


def cmd_figures(a):
    out, work = Path(a.out), Path(a.work)
    rows = read_tsv(out / 'ucd_pairs.tsv.gz')
    want = [('ucd_supported', 'clear pass'), ('same_context', 'same context'),
            ('only_one_side_different', 'one side different'), ('swapped_or_rearranged_context', 'orientation / swapped'),
            ('insufficient_context', 'insufficient context')]
    pick = []
    for v, title in want:
        for kind in ('candidate_insertion', 'cds_core'):
            c = sorted((r for r in rows if r['verdict'] == v and r['interval_kind'] == kind),
                       key=lambda r: (r['case'] == '', -float(r['insertion_identity'] or 0), r['pair_key']))
            if c:
                pick.append((c[0], title))
                break
    keys = {r['pair_key'] for r, _ in pick}
    hs = defaultdict(list)
    with gzip.open(work / 'evidence_hsps.tsv.gz', 'rt') as f:
        for h in csv.DictReader(f, delimiter='\t'):
            if h['pair_key'] in keys:
                hs[h['pair_key']].append(h)
    fd = out / 'figures'
    fd.mkdir(exist_ok=True)
    made = []
    for r, title in pick:
        name = '%s_%s.svg' % (r['verdict'], r['pair_key'].replace('|', '_'))
        (fd / name).write_text(svg_pair(r, hs[r['pair_key']], '%s: %s %s' % (title, r['pair_key'],
                                                                            r['case'] or r['product'][:50])))
        made.append(dict(file='figures/' + name, verdict=r['verdict'], pair_key=r['pair_key'], case=r['case']))
    (fd / 'index.json').write_text(json.dumps(made, indent=1) + '\n')
    print(json.dumps(made, indent=1))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = p.add_subparsers(dest='cmd', required=True)
    x = sp.add_parser('prepare')
    x.add_argument('--out', required=True)
    x.add_argument('--samples', required=True, help='comma-separated samples.tsv files (sample column)')
    x.add_argument('--round5', default=str(REPO / 'readcheck_round5' / 'pairs.tsv.gz'))
    x.add_argument('--cases', default=str(REPO / 'regression_b1-2-3.tsv'))
    x.add_argument('--reference', default=str(REPO / 's6_ucd_structural_pairs_b1-2-3.tsv.gz'),
                   help='earlier context classes, copied for comparison only')
    x.add_argument('--work', default=str(SCRATCH / 'scale_comamonadaceae_v0.1'))
    x.add_argument('--tag', default='b1-2-3')
    x = sp.add_parser('run')
    x.add_argument('--out', required=True)
    x.add_argument('--work', required=True, help='scratch dir for the sequence / HSP evidence')
    x.add_argument('--manifest', default=str(REPO / 'manifest.tsv'))
    x.add_argument('--config', default='', help='JSON overriding CONFIG keys')
    x.add_argument('--threads', type=int, default=int(os.environ.get('SLURM_CPUS_PER_TASK', 4)))
    x = sp.add_parser('figures')
    x.add_argument('--out', required=True)
    x.add_argument('--work', required=True)
    a = p.parse_args()
    dict(prepare=cmd_prepare, run=cmd_run, figures=cmd_figures)[a.cmd](a)


if __name__ == '__main__':
    main()
