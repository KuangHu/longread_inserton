"""Check 1 (source and extraction) and Check 2 (independent U-X-D comparison) used by ucd_batch.py.

Both checks are written separately from the ucd_compare.py screen (v0.1, frozen) and do not call its extraction or
coverage code. Check 2 shares only the frozen rule definitions: the CONFIG thresholds and the low-complexity rule
(ucd_compare.usable_bp).

Check 1  source and extraction
  The assembly FASTA is verified against the manifest sha256 (the checksum recorded when it was downloaded and
  verified). It is then indexed from its bytes: for each record, the byte offset, the length and the line widths.
  Every interval is read again with seek(), and the minus strand is reverse-complemented with a full IUPAC table.
  Coordinates are 0-based half-open on the forward strand of the contig (COORDS).
  The expected U / X / D intervals are recomputed from the input coordinates. Intervals and sequences must equal
  those the screen used, base by base.
  - Errors are reported, never repaired: FASTA checksum mismatch, contig missing, interval out of bounds or empty,
    irregular line widths, interval or sequence different from the screen.
  - A flank shortened by a contig edge is not an error. It is recorded and later judged by the usable-context rule.

Check 2  independent U-X-D comparison of the Check-1 sequences
  One blastn per pair and interval kind: A's X, U, D against B's X, U, D, on both strands.
  - Insertion: the reliable HSPs (>= ins_identity, >= min_hsp) on one strand are chained collinearly:
    - each next HSP starts and ends later in both sequences (on the minus strand, earlier in the subject);
    - an overlap of at most max_chain_overlap bp with the previous HSP is allowed, and overlapping bases are counted
      once.
    The best chain (largest identical bases) gives the coverage of X and of X'. No position is used twice and no
    out-of-order or inverted fragment is added.
      insertion_similar                    chain covers >= ins_coverage of both
      rearranged_or_fragmented_insertion   the chain does not, but the union of all reliable HSPs on both strands
                                           covers >= ins_coverage of both (same pieces, other order/orientation)
      insertion_partial                    some reliable HSP, less coverage
      insertion_dissimilar                 no reliable HSP
  - Orientation: if the best chain is on the minus strand, B is read reverse-complemented
    (U_B' = rc D_B, D_B' = rc U_B).
  - Contexts: the frozen v0.1 rule. Non-redundant union of reliable HSPs (>= ctx_identity, >= min_hsp) on each
    flank, for U/U and D/D (same side), the same inverted, and U_A/D_B', D_A/U_B' on both strands:
    similar both >= ctx_similar_cov; different both <= ctx_different_cov; otherwise uncertain.
  - Usable context (< min_usable_bp non-N, non-low-complexity bp: insufficient), and flanks that reuse one physical
    stretch of the same contig, as in v0.1.
  - Verdict, first that applies:
    insertion class
    > shared_physical_context
    > insufficient_context
    > swapped_or_rearranged_context
    > same_context
    > only_one_side_different
    > ucd_supported (no cross or inverted comparison similar or uncertain)
    > uncertain_context
"""
import hashlib
import subprocess
from pathlib import Path

import ucd_compare as uc

CHECK_VERSION = 'ucd_checks_v0.1'
CHECK_CONFIG = dict(max_chain_overlap=20)
COORDS = ('0-based half-open, forward strand of the contig as in the source FASTA; U/X/D sequences are given in the '
          'insertion orientation (minus strand = reverse complement of the forward interval)')
IUPAC = dict(zip('ACGTURYSWKMBDHVN', 'TGCAAYRSWMKVHDBN'))
FMT = '6 qseqid sseqid pident length nident qstart qend sstart send qlen slen evalue bitscore sstrand'
PAIRS8 = (('UU', 'U', 'U', 'plus'), ('DD', 'D', 'D', 'plus'), ('UU_inv', 'U', 'U', 'minus'),
          ('DD_inv', 'D', 'D', 'minus'), ('UD', 'U', 'D', 'plus'), ('UD_inv', 'U', 'D', 'minus'),
          ('DU', 'D', 'U', 'plus'), ('DU_inv', 'D', 'U', 'minus'))
CROSS = ('UU_inv', 'DD_inv', 'UD', 'UD_inv', 'DU', 'DU_inv')


# ---------------------------------------------------------------- Check 1

def revcomp(s):
    return ''.join(IUPAC.get(c, '?') for c in reversed(s))


class Source:
    """One assembly FASTA: sha256 verified against the manifest, and a byte index of its records."""

    def __init__(self, path, expect_sha256):
        self.path, self.error, self.index = str(path), '', {}
        try:
            data = Path(path).read_bytes()
        except OSError as e:
            self.error = 'fasta_unreadable:%s' % e.__class__.__name__
            return
        self.sha256 = hashlib.sha256(data).hexdigest()
        if self.sha256 != expect_sha256:
            self.error = 'fasta_sha256_mismatch'
            return
        name, pos, rec = None, 0, None
        for line in data.split(b'\n'):
            nxt = pos + len(line) + 1
            if line.startswith(b'>'):
                name = line[1:].split()[0].decode()
                rec = self.index[name] = dict(offset=nxt, length=0, widths=[], line_bytes=None)
            elif rec is not None and not line.rstrip(b'\r'):
                rec['blank'] = True
            elif rec is not None:
                if rec.get('blank'):
                    rec['bad'] = True
                body = line.rstrip(b'\r')
                rec['widths'].append(len(body))
                rec['line_bytes'] = rec['line_bytes'] or len(line) + 1
                rec['length'] += len(body)
            pos = nxt
        for r in self.index.values():
            w = r.pop('widths')
            r['line_bases'] = w[0] if w else 0
            r['regular'] = all(x == w[0] for x in w[:-1]) and (not w or w[-1] <= w[0]) and not r.pop('bad', False)
            r.pop('blank', None)

    def fetch(self, contig, lo, hi):
        """Forward-strand bases [lo, hi) read from the file by offset; raises ValueError naming the problem."""
        r = self.index.get(contig)
        if r is None:
            raise ValueError('contig_missing')
        if not r['regular']:
            raise ValueError('irregular_line_width')
        if not (0 <= lo < hi <= r['length']):
            raise ValueError('interval_out_of_bounds')
        lb, LB = r['line_bases'], r['line_bytes']
        a = r['offset'] + (lo // lb) * LB + lo % lb
        b = r['offset'] + ((hi - 1) // lb) * LB + (hi - 1) % lb + 1
        with open(self.path, 'rb') as f:
            f.seek(a)
            raw = f.read(b - a)
        s = raw.replace(b'\n', b'').replace(b'\r', b'').decode().upper()
        if len(s) != hi - lo:
            raise ValueError('short_read_from_source')
        return s


def check1_side(src, contig, lo, hi, strand, L, screen=None):
    """Re-extract U, X, D of one occurrence from the verified source and compare them with the screen's."""
    out = dict(status='fail', reasons=[], contig_len='', edge='')
    if src.error:
        out['reasons'].append(src.error)
        return out
    r = src.index.get(contig)
    if r is None:
        out['reasons'].append('contig_missing')
        return out
    n = r['length']
    out['contig_len'] = n
    if strand not in '+-' or not strand:
        out['reasons'].append('strand_invalid')
        return out
    if not (0 <= lo < hi <= n):
        out['reasons'].append('interval_out_of_bounds')
        return out
    left, right = (max(0, lo - L), lo), (hi, min(n, hi + L))
    ivs = dict(X=(lo, hi), U=left, D=right) if strand == '+' else dict(X=(lo, hi), U=right, D=left)
    seqs = {}
    try:
        for k, (a, b) in ivs.items():
            s = src.fetch(contig, a, b) if b > a else ''
            seqs[k] = s if strand == '+' else revcomp(s)
    except ValueError as e:
        out['reasons'].append(str(e))
        return out
    edge = [k for k in 'UD' if ivs[k][1] - ivs[k][0] < L]
    out['edge'] = ','.join('%s_at_contig_edge' % k for k in edge)
    if any('?' in s for s in seqs.values()):
        out['reasons'].append('unknown_base_symbol')
    if screen is not None:
        for k in 'XUD':
            if tuple(screen[k + '_iv']) != ivs[k]:
                out['reasons'].append('screen_interval_differs:%s' % k)
            elif screen[k] != seqs[k]:
                i = next((j for j, (p, q) in enumerate(zip(screen[k], seqs[k])) if p != q),
                         min(len(screen[k]), len(seqs[k])))
                out['reasons'].append('screen_sequence_differs:%s@%d' % (k, i))
    out.update(ivs=ivs, seqs=seqs, sha={k: hashlib.sha256(v.encode()).hexdigest()[:16] for k, v in seqs.items()})
    if not out['reasons']:
        out['status'] = 'pass'
    return out


# ---------------------------------------------------------------- Check 2

def blast_all(A, B, tmp, cfg):
    q, s = Path(tmp) / 'q.fa', Path(tmp) / 's.fa'
    q.write_text(''.join('>%s\n%s\n' % (k, A[k]) for k in 'XUD' if A[k]))
    s.write_text(''.join('>%s\n%s\n' % (k, B[k]) for k in 'XUD' if B[k]))
    if not q.stat().st_size or not s.stat().st_size:
        return []
    txt = subprocess.run(['blastn', '-task', 'blastn', '-query', str(q), '-subject', str(s), '-evalue',
                          str(cfg['evalue']), '-dust', 'no', '-outfmt', FMT], stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE, universal_newlines=True, check=True).stdout
    out = []
    for line in txt.splitlines():
        f = line.split('\t')
        ss, se = int(f[7]), int(f[8])
        out.append(dict(q=f[0], s=f[1], pident=float(f[2]), length=int(f[3]), nident=int(f[4]), qs=int(f[5]) - 1,
                        qe=int(f[6]), ss=min(ss, se) - 1, se=max(ss, se), qlen=int(f[9]), slen=int(f[10]),
                        evalue=float(f[11]), bitscore=float(f[12]), strand='plus' if f[13] == 'plus' else 'minus'))
    return out


def merged_len(ivs):
    tot, end = 0, None
    for a, b in sorted(ivs):
        if end is None or a > end:
            tot += b - a
            end = b
        elif b > end:
            tot += b - end
            end = b
    return tot


def reliable(h, q, s, strand, min_id, cfg):
    return [x for x in h if x['q'] == q and x['s'] == s and x['strand'] == strand and x['pident'] >= min_id
            and x['length'] >= cfg['min_hsp']]


def best_chain(hs, strand, max_ov):
    """Collinear chain of HSPs with the most identical bases (O(n^2); n is small)."""
    hs = sorted(hs, key=lambda x: (x['qs'], x['qe']))

    def follows(p, n):
        if not (n['qs'] > p['qs'] and n['qe'] > p['qe'] and p['qe'] - n['qs'] <= max_ov):
            return False
        if strand == 'plus':
            return n['ss'] > p['ss'] and n['se'] > p['se'] and p['se'] - n['ss'] <= max_ov
        return n['ss'] < p['ss'] and n['se'] < p['se'] and n['se'] - p['ss'] <= max_ov

    best, prev = [], []
    for i, h in enumerate(hs):
        sc, pv = h['nident'], -1
        for j in range(i):
            if follows(hs[j], h) and best[j] + h['nident'] > sc:
                sc, pv = best[j] + h['nident'], j
        best.append(sc)
        prev.append(pv)
    if not hs:
        return []
    i = max(range(len(hs)), key=lambda k: best[k])
    chain = []
    while i >= 0:
        chain.append(hs[i])
        i = prev[i]
    return chain[::-1]


def cover(hs, qlen, slen):
    if not hs or not qlen or not slen:
        return dict(identity=0.0, cov_a=0.0, cov_b=0.0, hsps=len(hs))
    return dict(identity=round(100.0 * sum(x['nident'] for x in hs) / sum(x['length'] for x in hs), 2),
                cov_a=round(merged_len([(x['qs'], x['qe']) for x in hs]) / qlen, 3),
                cov_b=round(merged_len([(x['ss'], x['se']) for x in hs]) / slen, 3), hsps=len(hs))


def ctx_class(m, cfg):
    if min(m['cov_a'], m['cov_b']) >= cfg['ctx_similar_cov']:
        return 'similar'
    if max(m['cov_a'], m['cov_b']) <= cfg['ctx_different_cov']:
        return 'different'
    return 'uncertain'


def check_ucd(A, B, cfg, ccfg, tmp):
    """A, B: dict(id, X, U, D, sample, contig, ivs) from Check 1, each in its own insertion orientation."""
    if A['id'] > B['id']:
        A, B = B, A
    h = blast_all(A, B, tmp, cfg)
    la, lb = {k: len(A[k]) for k in 'XUD'}, {k: len(B[k]) for k in 'XUD'}
    ins = {}
    for st in ('plus', 'minus'):
        ch = best_chain(reliable(h, 'X', 'X', st, cfg['ins_identity'], cfg), st, ccfg['max_chain_overlap'])
        ins[st] = dict(cover(ch, la['X'], lb['X']), chain=ch)
    ori = max(('plus', 'minus'), key=lambda st: (min(ins[st]['cov_a'], ins[st]['cov_b']), st == 'plus'))
    mi = ins[ori]
    allx = reliable(h, 'X', 'X', 'plus', cfg['ins_identity'], cfg) + \
        reliable(h, 'X', 'X', 'minus', cfg['ins_identity'], cfg)
    ux = cover(allx, la['X'], lb['X'])
    if min(mi['cov_a'], mi['cov_b']) >= cfg['ins_coverage']:
        ins_cls = 'insertion_similar'
    elif min(ux['cov_a'], ux['cov_b']) >= cfg['ins_coverage']:
        ins_cls = 'rearranged_or_fragmented_insertion'
    elif allx:
        ins_cls = 'insertion_partial'
    else:
        ins_cls = 'insertion_dissimilar'
    flip = ori == 'minus'
    side_b = {'U': 'D', 'D': 'U'} if flip else {'U': 'U', 'D': 'D'}
    other = {'plus': 'minus', 'minus': 'plus'}
    cmp = {}
    for name, qa, sb, st in PAIRS8:
        sk, sst = side_b[sb], (other[st] if flip else st)
        m = cover(reliable(h, qa, sk, sst, cfg['ctx_identity'], cfg), la[qa], lb[sk])
        m['class'] = ctx_class(m, cfg)
        cmp[name] = m
    Bo = dict(U=revcomp(B['D']) if flip else B['U'], D=revcomp(B['U']) if flip else B['D'])
    biv = dict(X=B['ivs']['X'], U=B['ivs'][side_b['U']], D=B['ivs'][side_b['D']])
    usable = dict(U_A=uc.usable_bp(A['U'], cfg), D_A=uc.usable_bp(A['D'], cfg), U_B=uc.usable_bp(Bo['U'], cfg),
                  D_B=uc.usable_bp(Bo['D'], cfg))
    short = [k for k in ('D_A', 'D_B', 'U_A', 'U_B') if usable[k] < cfg['min_usable_bp']]
    phys = []
    if A.get('sample') and A['sample'] == B.get('sample') and A['contig'] == B['contig']:
        for p in 'UDX':
            for q in 'UDX':
                a, b = A['ivs'][p], biv[q]
                if a[0] < b[1] and b[0] < a[1] and a[1] > a[0] and b[1] > b[0]:
                    phys.append('%s_A~%s_B' % (p, q))
    cross = [k for k in CROSS if cmp[k]['class'] == 'similar']
    cross_unc = [k for k in CROSS if cmp[k]['class'] == 'uncertain']
    u, d = cmp['UU']['class'], cmp['DD']['class']
    reasons = ([ins_cls] if ins_cls != 'insertion_similar' else []) + \
        (['shared_physical_context:' + ','.join(phys)] if phys else []) + \
        (['insufficient_context:' + ','.join(short)] if short else []) + \
        (['cross_or_inverted_similar:' + ','.join(cross)] if cross else []) + \
        (['cross_or_inverted_uncertain:' + ','.join(cross_unc)] if cross_unc else []) + \
        ['upstream_' + u, 'downstream_' + d]
    if ins_cls != 'insertion_similar':
        v = ins_cls
    elif phys:
        v = 'shared_physical_context'
    elif short:
        v = 'insufficient_context'
    elif cross:
        v = 'swapped_or_rearranged_context'
    elif u == d == 'similar':
        v = 'same_context'
    elif {u, d} == {'similar', 'different'}:
        v = 'only_one_side_different'
    elif u == d == 'different' and not cross_unc:
        v = 'ucd_supported'
    else:
        v = 'uncertain_context'
    row = dict(chk_occurrence_a=A['id'], chk_occurrence_b=B['id'], chk_b_flipped=flip, chk_insertion_class=ins_cls,
               chk_insertion_identity=mi['identity'], chk_insertion_cov_a=mi['cov_a'],
               chk_insertion_cov_b=mi['cov_b'], chk_insertion_chain_hsps=mi['hsps'], chk_insertion_orientation=ori,
               chk_insertion_union_cov_a=ux['cov_a'], chk_insertion_union_cov_b=ux['cov_b'],
               **{'chk_%s_usable' % k: v_ for k, v_ in usable.items()},
               **{'chk_%s_%s' % (c, f): cmp[c][f] for c in cmp for f in ('class', 'identity', 'cov_a', 'cov_b')},
               chk_upstream=u, chk_downstream=d, ucd_check=v, chk_reasons=';'.join(reasons))
    chain_ids = {id(x) for x in mi['chain']}
    hsps = [dict(x, in_insertion_chain=id(x) in chain_ids) for x in h]
    return row, hsps
