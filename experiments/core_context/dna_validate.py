"""Validation of the DNA entry on a set of real FNA: the two checks on the DNA-entry pairs, and a regression comparison
with the protein entry's strict pairs between the same assemblies. The protein-entry pairs are a reference, not the
complete truth; the DNA entry is not required to reproduce their number, and its own new candidates are kept.

v0.2 correspondence rules (fixed; not tuned to raise recovery). No family, core or pair ID is used.
  side correspondence of a protein X (sample, contig, [s, e)) and a DNA X (same sample and contig):
    same_interval        |delta start| <= TOL and |delta end| <= TOL bp
    boundary_differs     otherwise, reciprocal overlap >= MIN_OVERLAP of both intervals (same physical locus)
    uncertain            overlap > 0 but below MIN_OVERLAP
  pair correspondence: both sides correspond (A-B or B-A) and the relative orientation of the two copies (same /
    opposite strand) agrees; same_interval if both sides are, else boundary_differs; a pair whose sides only overlap
    (uncertain) or whose relative orientation disagrees is correspondence_uncertain.

Pair-level class of each protein strict pair (first that applies):
  recovered_same_interval / recovered_boundary_differs   a corresponding DNA row passes both checks, both X < 2500
  not_recovered:length_changed           a corresponding DNA row passes both checks, an X >= 2500 bp
  not_recovered:checks_not_passed        corresponding DNA rows checked, none passes (reason: ucd_check)
  not_recovered:correspondence_uncertain DNA rows overlap both X, but only as `uncertain` or with other orientation
  one_side_recovered:<reason>            one X is a side of a DNA strict-pass row (with another partner); the
                                         pair itself is classed by the reasons below
  not_recovered:boundary_unresolved      a DNA pair whose seeds overlap both X is unresolved (reason: ends)
  not_recovered:pair_deferred            both X have occurrences of one library sequence whose pairs were deferred
  not_recovered:occurrences_not_paired   both X have occurrences of one library sequence, nothing deferred (a gap)
  not_recovered:occurrence_not_detected  occurrences on one / both sides, but of no common library sequence; or a
                                         seed locus overlaps X but no occurrence
  not_recovered:no_seed                  neither side overlaps any seed locus or occurrence
Occurrence level: every distinct X of the protein strict pairs: in a DNA strict-pass row / in a checked DNA row / in a
DNA occurrence (any library hit overlapping it) / in a seed locus only / none.
"""
import argparse
import gzip
import json
from collections import Counter, defaultdict
from pathlib import Path

import ucd_compare as uc

VERSION = 'dna_validate_v0.2.1'  # v0.2.1: reason of the pair for one_side_recovered
MIN_OVERLAP = 0.5
TOL = 10
STRICT_MAX = 2500


def iv(s):
    a, b = s.split('-')
    return int(a), int(b)


def side_corr(x, y):
    if abs(x[0] - y[0]) <= TOL and abs(x[1] - y[1]) <= TOL:
        return 'same_interval'
    ov = min(x[1], y[1]) - max(x[0], y[0])
    if ov <= 0:
        return None
    if ov >= MIN_OVERLAP * (x[1] - x[0]) and ov >= MIN_OVERLAP * (y[1] - y[0]):
        return 'boundary_differs'
    return 'uncertain'


def stream(path, keep=lambda r: True):
    with gzip.open(path, 'rt') as f:
        h = f.readline().rstrip('\n').split('\t')
        for line in f:
            r = dict(zip(h, line.rstrip('\n').split('\t')))
            if keep(r):
                yield r


def strict(r):
    return r['final_status'] == 'pass' and 0 < int(r['x_len_a'] or 0) < STRICT_MAX and \
        0 < int(r['x_len_b'] or 0) < STRICT_MAX


def lkey(r):
    return tuple(sorted(((r['a_sample'], r['a_contig']), (r['b_sample'], r['b_contig']))))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--manifest', required=True)
    p.add_argument('--dna-work', required=True)
    p.add_argument('--dna-collect', required=True)
    p.add_argument('--dna-ucd', required=True, help='ucd_batch out dir of the DNA pairs (completeness.json)')
    p.add_argument('--dna-rows', required=True, help='ucd_batch rows_candidate_insertion of the DNA pairs')
    p.add_argument('--protein-final', required=True)
    p.add_argument('--out', required=True)
    a = p.parse_args()
    S = {m['sample'] for m in uc.read_tsv(a.manifest)}
    inset = lambda r: r['a_sample'] in S and r['b_sample'] in S
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    work = Path(a.dna_work)

    # ---------------- DNA run
    keep_cols = ('pair_id', 'pair_scope', 'final_status', 'source_check', 'ucd_check', 'x_len_a', 'x_len_b',
                 'a_sample', 'a_contig', 'a_strand', 'a_x_interval', 'b_sample', 'b_contig', 'b_strand',
                 'b_x_interval')
    dna = [{k: r[k] for k in keep_cols} for r in stream(a.dna_rows)]
    by_loc = defaultdict(list)
    for r in dna:
        by_loc[lkey(r)].append(r)
    st = [r for r in dna if strict(r)]
    socc = {(r[s + '_sample'], r[s + '_contig'], r[s + '_x_interval']) for r in st for s in 'ab'}
    side_strict = defaultdict(list)  # (sample, contig) -> X intervals of strict-pass sides
    for x in socc:
        side_strict[(x[0], x[1])].append(iv(x[2]))
    coll = json.loads((Path(a.dna_collect) / 'collect.json').read_text())
    comp = json.loads((Path(a.dna_ucd) / 'completeness.json').read_text())
    D = dict(collect_pair_rows=coll['pair_rows'], collect_status=coll['status'],
             resolved_distinct_pairs_input=coll['resolved_distinct'],
             resolved_duplicate_rows=coll['resolved_duplicate_rows'],
             pairs_deferred_rows_high_copy=coll['pairs_deferred_high_copy'],
             checked_rows=len(dna), completeness=comp,
             source_check=Counter(r['source_check'] for r in dna), ucd_check=Counter(r['ucd_check'] for r in dna),
             final=Counter(r['final_status'] for r in dna),
             checked_by_scope=Counter(r['pair_scope'] for r in dna),
             strict_pass_pairs=len(st), strict_pass_by_scope=Counter(r['pair_scope'] for r in st),
             strict_pass_by_strands=Counter(r['a_strand'] + r['b_strand'] for r in st),
             strict_pass_occurrences=len(socc),
             strict_pass_assemblies=len({x[0] for x in socc}),
             pass_not_strict=sum(r['final_status'] == 'pass' and not strict(r) for r in dna))

    # unresolved / no-block pair rows indexed by locus key (seeds), from the pair table
    unres = defaultdict(list)
    for r in stream(Path(a.dna_collect) / 'dna_pairs_all.tsv.gz', lambda r: r['status'] != 'resolved'):
        unres[lkey(r)].append({k: r[k] for k in ('a_sample', 'a_contig', 'a_seed', 'b_sample', 'b_contig',
                                                  'b_seed', 'status', 'boundary_up', 'boundary_down')})
    # occurrences, seed loci, deferred library sequences
    occ = defaultdict(list)
    for f in work.glob('occ/lib*/*/b*.tsv.gz'):
        with gzip.open(f, 'rt') as fh:
            fh.readline()
            for line in fh:
                v = line.split('\t', 6)
                occ[(v[1], v[2])].append((int(v[4]), int(v[5]), v[0]))
    seeds = defaultdict(list)
    for f in work.glob('seeds/*/loci.tsv.gz'):
        for r in stream(f):
            seeds[(r['sample'], r['contig'])].append((int(r['lo']), int(r['hi'])))
    deferred = {}
    for f in work.glob('pairs/shard_*/budget.tsv.gz'):  # v0.1.2 pairs stage
        for r in stream(f):
            if int(r['pairs_deferred']):
                deferred[r['library_id']] = int(r['pairs_deferred'])
    done_chunks, pend = set(), {}
    for f in work.glob('enum/g*/b*/progress.tsv.gz'):  # v0.2 enumerate: chunks never enumerated
        for r in stream(f):
            if r['status'] == 'enumerated':
                done_chunks.add((r['library_id'], r['chunk']))
            else:
                pend[(r['library_id'], r['chunk'])] = int(r['pairs'])
    for (l, c), n in pend.items():
        if (l, c) not in done_chunks:
            deferred[l] = deferred.get(l, 0) + n

    def ovl(lst, x):
        return [t for t in lst if min(t[1], x[1]) > max(t[0], x[0])]

    def pair_corr(pr, r):
        P = [(pr['a_sample'], pr['a_contig'], iv(pr['a_x_interval'])),
             (pr['b_sample'], pr['b_contig'], iv(pr['b_x_interval']))]
        Q = [(r['a_sample'], r['a_contig'], iv(r['a_x_interval'])) if r['a_x_interval'] else None,
             (r['b_sample'], r['b_contig'], iv(r['b_x_interval'])) if r['b_x_interval'] else None]
        if None in Q:
            return None
        best = None
        for q in (Q, Q[::-1]):
            if all(P[i][:2] == q[i][:2] for i in (0, 1)):
                c = [side_corr(P[i][2], q[i][2]) for i in (0, 1)]
                if None in c:
                    continue
                if 'uncertain' in c:
                    k = 'uncertain'
                else:
                    k = 'same_interval' if c == ['same_interval'] * 2 else 'boundary_differs'
                rank = {'same_interval': 0, 'boundary_differs': 1, 'uncertain': 2}[k]
                if best is None or rank < best[0]:
                    best = (rank, k)
        if best is None:
            return None
        rel_p = pr['a_strand'] == pr['b_strand']
        rel_d = r['a_strand'] == r['b_strand']
        if best[1] != 'uncertain' and rel_p != rel_d:
            return 'orientation_differs'
        return best[1]

    def classify(pr):
        rows = by_loc.get(lkey(pr), [])
        cs = [(pair_corr(pr, r), r) for r in rows]
        good = [(c, r) for c, r in cs if c in ('same_interval', 'boundary_differs')]
        if good:
            sp = [(c, r) for c, r in good if strict(r)]
            if sp:
                c = 'same_interval' if any(c == 'same_interval' for c, _ in sp) else 'boundary_differs'
                return 'recovered_' + c, [r['pair_id'] for _, r in sp]
            if any(r['final_status'] == 'pass' for _, r in good):
                return 'not_recovered:length_changed', [r['pair_id'] for _, r in good]
            return ('not_recovered:checks_not_passed:' + sorted(r['ucd_check'] if r['source_check'] == 'pass'
                                                                else 'source_' + r['source_check'] for _, r in good)[0],
                    [r['pair_id'] for _, r in good])
        if any(c in ('uncertain', 'orientation_differs') for c, _ in cs):
            return 'not_recovered:correspondence_uncertain', [r['pair_id'] for c, r in cs if c]
        A = (pr['a_sample'], pr['a_contig'], iv(pr['a_x_interval']))
        B = (pr['b_sample'], pr['b_contig'], iv(pr['b_x_interval']))
        one = [s for s, X in (('a', A), ('b', B)) if any(side_corr(X[2], y) in ('same_interval', 'boundary_differs')
                                                          for y in side_strict.get(X[:2], []))]

        def seed_ov(r, s, X):
            return r[s + '_sample'] == X[0] and r[s + '_contig'] == X[1] and bool(ovl([iv(r[s + '_seed'])], X[2]))
        un = [r for r in unres.get(lkey(pr), []) if (seed_ov(r, 'a', A) and seed_ov(r, 'b', B)) or
              (seed_ov(r, 'a', B) and seed_ov(r, 'b', A))]
        if un:
            u = un[0]
            why = 'boundary_unresolved:%s/%s' % (u['boundary_up'], u['boundary_down'])
        else:
            why = absent(A, B)
        return ('one_side_recovered:' if one else 'not_recovered:') + why, []

    def absent(A, B):
        la = {t[2] for t in ovl(occ.get(A[:2], []), A[2])}
        lb = {t[2] for t in ovl(occ.get(B[:2], []), B[2])}
        common = la & lb
        if common:
            if any(l in deferred for l in common):
                return 'pair_deferred'
            return 'occurrences_not_paired'
        if la and lb:
            return 'occurrence_not_detected:no_common_library_sequence'
        if la or lb:
            return 'occurrence_not_detected:one_side'
        if ovl(seeds.get(A[:2], []), A[2]) or ovl(seeds.get(B[:2], []), B[2]):
            return 'occurrence_not_detected:seed_only'
        return 'no_seed'

    prot = list(stream(Path(a.protein_final) / 'strict_lt2500_candidate_insertion.tsv.gz', inset))
    rec, rows = Counter(), []
    by_scope, by_strand = defaultdict(Counter), defaultdict(Counter)
    for pr in prot:
        c, ids = classify(pr)
        top = ':'.join(c.split(':')[:2])
        rec[c] += 1
        by_scope[pr['pair_scope']][top] += 1
        by_strand[pr['a_strand'] + pr['b_strand']][top] += 1
        rows.append(dict(protein_pair_id=pr['pair_id'], pair_scope=pr['pair_scope'], a_sample=pr['a_sample'],
                         a_contig=pr['a_contig'], a_x=pr['a_x_interval'], a_strand=pr['a_strand'],
                         b_sample=pr['b_sample'], b_contig=pr['b_contig'], b_x=pr['b_x_interval'],
                         b_strand=pr['b_strand'], dna_class=c, dna_pair_ids=';'.join(sorted(set(ids)))[:500]))
    cols = list(rows[0]) if rows else ['protein_pair_id']
    uc.write_tsv(out / 'recovery_pairs.tsv.gz', rows, cols)

    # occurrence level
    pocc = {(pr[s + '_sample'], pr[s + '_contig'], pr[s + '_x_interval']) for pr in prot for s in 'ab'}
    checked_side = defaultdict(list)
    for r in dna:
        for s in 'ab':
            if r[s + '_x_interval']:
                checked_side[(r[s + '_sample'], r[s + '_contig'])].append(iv(r[s + '_x_interval']))
    oc = Counter()
    for x in pocc:
        X = iv(x[2])
        good = lambda lst: any(side_corr(X, y) in ('same_interval', 'boundary_differs') for y in lst)
        if good(side_strict.get(x[:2], [])):
            oc['in_dna_strict_pass_row'] += 1
        elif good(checked_side.get(x[:2], [])):
            oc['in_checked_dna_row_only'] += 1
        elif ovl(occ.get(x[:2], []), X):
            oc['in_dna_occurrence_only'] += 1
        elif ovl(seeds.get(x[:2], []), X):
            oc['in_seed_locus_only'] += 1
        else:
            oc['none'] += 1

    # DNA strict-pass pairs without a protein-entry strict counterpart (DNA-only candidates, kept)
    pidx = defaultdict(list)
    for pr in prot:
        pidx[lkey(pr)].append(pr)
    dna_only = 0
    for r in st:
        if not any(pair_corr(pr, r) in ('same_interval', 'boundary_differs') for pr in pidx.get(lkey(r), [])):
            dna_only += 1

    R = dict(protein_strict_pairs=len(prot), classes=dict(rec.most_common()),
             top_classes=dict(Counter(':'.join(k.split(':')[:2]) for k in rec.elements()).most_common()),
             by_scope={k: dict(v) for k, v in by_scope.items()},
             by_strands={k: dict(v) for k, v in by_strand.items()},
             occurrences=dict(protein_strict_occurrences=len(pocc), **oc),
             dna_strict_pass_pairs_without_protein_strict_counterpart=dna_only)
    ctl = Counter()
    for pr in stream(Path(a.protein_final) / 'all_rows_candidate_insertion.tsv.gz',
                     lambda r: inset(r) and r['ucd_check'] in ('same_context', 'only_one_side_different')):
        c, _ = classify(pr)
        ctl['%s -> %s' % (pr['ucd_check'], ':'.join(c.split(':')[:2]))] += 1
    C = dict(protein_control_rows=sum(ctl.values()), dna_judgement=dict(ctl.most_common()))
    res = dict(version=VERSION, assemblies=len(S), rules=dict(tol_bp=TOL, min_overlap=MIN_OVERLAP,
                                                               strict_max_bp=STRICT_MAX),
               dna=D, recovery=R, controls=C)
    res = json.loads(json.dumps(res, default=lambda o: dict(o)))
    (out / 'validation.json').write_text(json.dumps(res, indent=1) + '\n')
    print(json.dumps(res, indent=1))


if __name__ == '__main__':
    main()
