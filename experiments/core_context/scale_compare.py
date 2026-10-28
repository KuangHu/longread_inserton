"""Exact, resumable pair accounting for protein-entry cores at scale (replaces the MAX_PAIRS=20000 truncation of
protein_context.compare for the scale-up; the frozen subset outputs are unchanged). Evidence stays assembly_only.

The K=3 pair relation of protein_context.relation depends only on each occurrence's context signature per side: the
set of neighbour labels, whether every neighbour is labelled, and whether the side has K neighbours. Occurrences are
collapsed to unique signatures with per-assembly multiplicities, and every unordered signature pair (including a
signature with itself) is related once and weighted by the number of occurrence pairs it stands for. Every occurrence
pair of a core is therefore counted exactly once, split into same-assembly and cross-assembly pairs, with no cap and
no dependence on input order.

Work is split into blocks of signature rows (at most BLOCK_PAIRS signature pairs each). A finished block is appended
to the shard's progress.jsonl with the digest of the core's inputs; a rerun skips finished blocks of unchanged cores
and recomputes cores whose inputs changed (e.g. a new batch). --max-sig-pairs bounds one run. Per core the output
records pairs_total, pairs_accounted and pairs_not_compared; a pair in an unfinished block is never counted as
similar, different or absent.

Per signature (not per occurrence pair) the output keeps what is needed to name partners per occurrence, because one
signature can hold occurrences of several assemblies:
  partner_occ        different_context_candidate partner occurrences of any member (sum over partner signatures)
  partner_by_sample  {assembly of a member: partner occurrences in that assembly}, so each occurrence gets its own
                     same-assembly count (partner_by_sample[own]) and cross-assembly count (partner_occ - that)
  witnesses          up to two (partner signature, assembly) with distinct assemblies, smallest first: whatever the
                     member's assembly, one of them is in another assembly whenever a cross partner exists
  same_first         {assembly of a member: smallest partner signature holding that assembly}
All three merge exactly across blocks (sums, per-key minima, the two smallest distinct-assembly witnesses).
(v1 kept one cross and one same partner signature per signature; a member whose assembly did not fit that partner
got no named partner although it had one. Pair counts were not affected.) Context classes are unions of
similar_context links (provisional, as in protein entry; not S6 site groups).
"""
import hashlib
import json
from collections import Counter, defaultdict

RELATIONS = ['similar_context', 'different_context_candidate', 'one_side_shared_unresolved', 'insufficient_context_labels',
             'edge_or_missing_context']
K = 3
BLOCK_PAIRS = 2_000_000
RECORD_VERSION = 2


def side_sig(labels):
    """labels: list of neighbour labels (None = unlabelled), nearest first."""
    return (frozenset(l for l in labels if l), len(labels) >= K, all(labels))


def relation_sig(x, y):
    """Same rule as protein_context.relation on (up side signature, down side signature) pairs."""
    (xu, xd), (yu, yd) = x, y
    if not (xu[1] and xd[1] and yu[1] and yd[1]):
        return 'edge_or_missing_context'
    shared = [bool(xu[0] & yu[0]), bool(xd[0] & yd[0])]
    if all(shared):
        return 'similar_context'
    labelled = [xu[2] and yu[2], xd[2] and yd[2]]
    if any(not sh and not lb for sh, lb in zip(shared, labelled)):
        return 'insufficient_context_labels'
    if not any(shared):
        return 'different_context_candidate'
    return 'one_side_shared_unresolved'


def sig_repr(s):
    return json.dumps([[sorted(side[0]), side[1], side[2]] for side in s])


def collapse(occs):
    """occs: [(occurrence_id, sample, up labels, down labels)]. Returns signatures (sorted by repr), per-signature
    sample multiplicities, occurrence -> signature index, and a digest of the inputs."""
    by = defaultdict(list)
    for oid, smp, up, down in occs:
        by[sig_repr((side_sig(up), side_sig(down)))].append((oid, smp, up, down))
    keys = sorted(by)
    sigs = [(side_sig(by[k][0][2]), side_sig(by[k][0][3])) for k in keys]
    mult = [Counter(x[1] for x in by[k]) for k in keys]
    index = {x[0]: i for i, k in enumerate(keys) for x in by[k]}
    h = hashlib.sha256()
    for k in keys:
        h.update(k.encode()); h.update(json.dumps(sorted((x[0], x[1]) for x in by[k])).encode())
    return sigs, mult, index, h.hexdigest()


def blocks(mult, limit=BLOCK_PAIRS):
    """Row ranges [i0, i1) of the upper signature triangle with at most `limit` signature pairs each."""
    u, out, i0, cur = len(mult), [], 0, 0
    for i in range(u):
        cur += u - i
        if cur >= limit:
            out.append((i0, i + 1)); i0, cur = i + 1, 0
    if i0 < u:
        out.append((i0, u))
    return out


def pairs_of(mult, i, j):
    """(occurrence pairs, same-assembly occurrence pairs) represented by signature pair (i, j)."""
    if i == j:
        n = sum(mult[i].values())
        return n * (n - 1) // 2, sum(m * (m - 1) // 2 for m in mult[i].values())
    a, b = (mult[i], mult[j]) if len(mult[i]) <= len(mult[j]) else (mult[j], mult[i])
    return sum(mult[i].values()) * sum(mult[j].values()), sum(m * b.get(s, 0) for s, m in a.items())


def run_block(sigs, mult, i0, i1):
    """Relations of rows [i0, i1) against every j >= i. Returns a JSON-serialisable block record."""
    total, same = Counter(), Counter()
    parent = {}

    def root(x):
        while parent.get(x, x) != x:
            parent[x] = parent.get(parent[x], parent[x]); x = parent[x]
        return x
    unions, diff = [], {}
    for i in range(i0, i1):
        for j in range(i, len(sigs)):
            n, ss = pairs_of(mult, i, j)
            if not n:
                continue
            r = relation_sig(sigs[i], sigs[j])
            total[r] += n; same[r] += ss
            if r == 'similar_context' and i != j:
                ri, rj = root(i), root(j)
                if ri != rj:
                    parent[ri] = rj; unions.append((i, j))
            elif r == 'different_context_candidate':
                for x, y in ((i, j), (j, i)) if i != j else ((i, i),):
                    d = diff.setdefault(x, new_partner())
                    add_partner(d, x, y, mult)
    return dict(v=RECORD_VERSION, i0=i0, i1=i1, total=dict(total), same=dict(same), unions=unions,
                diff={str(k): v for k, v in diff.items()}, pairs=sum(total.values()))


def new_partner():
    return dict(n=0, ss=0, partner_occ=0, partner_by_sample={}, witnesses=[], same_first={})


def keep_witnesses(ws):
    """The two smallest (signature, assembly) entries with distinct assemblies."""
    out, seen = [], set()
    for j, smp in sorted(ws):
        if smp not in seen:
            out.append([j, smp]); seen.add(smp)
        if len(out) == 2:
            break
    return out


def add_partner(d, x, y, mult):
    """Signature y is a different-context partner signature of x (x == y: pairs inside one signature)."""
    n, ss = pairs_of(mult, x, y)
    d['n'] += n; d['ss'] += ss
    d['partner_occ'] += sum(mult[y].values())
    for smp, m in mult[x].items():
        k = mult[y].get(smp, 0) - (1 if x == y else 0)
        if k > 0:
            d['partner_by_sample'][smp] = d['partner_by_sample'].get(smp, 0) + k
            if smp not in d['same_first'] or y < d['same_first'][smp]:
                d['same_first'][smp] = y
    d['witnesses'] = keep_witnesses(d['witnesses'] + [[y, smp] for smp in sorted(mult[y])[:2]])


def merge_partner(d, e):
    d['n'] += e['n']; d['ss'] += e['ss']; d['partner_occ'] += e['partner_occ']
    for smp, k in e['partner_by_sample'].items():
        d['partner_by_sample'][smp] = d['partner_by_sample'].get(smp, 0) + k
    for smp, j in e['same_first'].items():
        if smp not in d['same_first'] or j < d['same_first'][smp]:
            d['same_first'][smp] = j
    d['witnesses'] = keep_witnesses(d['witnesses'] + e['witnesses'])


def occurrence_partners(row, smp):
    """Per occurrence of assembly smp in a signature: (cross partners, same-assembly partners, cross witness
    (signature, assembly) or None, same-assembly partner signature or None)."""
    same_n = row['partner_by_sample'].get(smp, 0)
    cross_n = row['partner_occ'] - same_n
    wit = next(((j, s) for j, s in row['witnesses'] if s != smp), None) if cross_n > 0 else None
    return cross_n, same_n, wit, row['same_first'].get(smp) if same_n > 0 else None


def name_partners(occs, index, rows):
    """occs: [(occurrence_id, assembly)] of one core; index: occurrence -> signature; rows: finalize() rows.
    Returns {occurrence: dict(cross partners, same-assembly partners, a named partner of each kind or '')}."""
    members = defaultdict(list)
    for oid, smp in sorted(occs):
        members[index[oid]].append((oid, smp))
    out = {}
    for oid, smp in occs:
        cross_n, same_n, wit, same_sig = occurrence_partners(rows[index[oid]], smp)
        cross = next(x for x, s in members[wit[0]] if s == wit[1]) if wit else ''
        same = next(x for x, s in members[same_sig] if s == smp and x != oid) if same_sig is not None else ''
        out[oid] = dict(diff_partners_cross=cross_n, diff_partners_same_assembly=same_n,
                        diff_partner_cross=cross, diff_partner_same=same)
    return out


def finalize(sigs, mult, records):
    """Merge finished block records of one core. Returns (core summary, per-signature rows)."""
    n_occ = sum(sum(m.values()) for m in mult)
    total_pairs = n_occ * (n_occ - 1) // 2
    total, same = Counter(), Counter()
    parent = list(range(len(sigs)))

    def root(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    diff = defaultdict(new_partner)
    for rec in records:
        total.update(rec['total']); same.update(rec['same'])
        for i, j in rec['unions']:
            ri, rj = root(i), root(j)
            if ri != rj:
                parent[ri] = rj
        for k, e in rec['diff'].items():
            merge_partner(diff[int(k)], e)
    accounted = sum(total.values())
    done_rows = sum(r['i1'] - r['i0'] for r in records)
    summary = dict(occurrences=n_occ, assemblies=len(set().union(*[set(m) for m in mult])) if mult else 0,
                   signatures=len(sigs), pairs_total=total_pairs, pairs_accounted=accounted,
                   pairs_not_compared=total_pairs - accounted, complete=done_rows == len(sigs) and accounted == total_pairs,
                   **{r: total[r] for r in RELATIONS}, **{r + '_same_assembly': same[r] for r in RELATIONS},
                   context_classes_provisional=len({root(i) for i in range(len(sigs))}) if done_rows == len(sigs) else '')
    rows = []
    for i in range(len(sigs)):
        d = diff[i] if i in diff else new_partner()
        rows.append(dict(sig=i, context_class=root(i), diff_pairs=d['n'], diff_pairs_same_assembly=d['ss'],
                         partner_occ=d['partner_occ'], partner_by_sample=d['partner_by_sample'],
                         witnesses=d['witnesses'], same_first=d['same_first']))
    return summary, rows
