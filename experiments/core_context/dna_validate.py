"""Validation of the DNA entry on a small set of real FNA: correctness of the two checks on DNA-entry pairs and recovery
of known U-X-D pairs. Not a requirement to reproduce the protein entry's candidates or counts.

Inputs: the DNA run (dna_discovery work dir and collect dir, ucd_batch rows of the DNA pairs) and the protein entry's
final rows (all_rows_candidate_insertion, strict table) restricted to the validation assemblies.

Reported:
  dna            candidate pairs, resolved / unresolved (with reasons), checked rows by source_check and ucd_check,
                 pass, strict (0 < both X < 2500 bp), by pair scope and strand combination.
  recovery       every protein-entry strict pair between validation assemblies, matched to DNA rows by the same two
                 assemblies and contigs and reciprocal overlap >= MIN_OVERLAP of both X intervals:
                   dna_strict_pass / dna_pass_not_strict / dna_checked_not_passed:<verdict> /
                   dna_pair_unresolved:<reason> (a DNA pair whose seeds overlap both X but whose boundary is open) /
                   occurrences_not_paired (both X overlap DNA occurrences of one library sequence, no pair formed) /
                   occurrence_one_side / not_found
                 counted overall, by scope and by strand combination.
  controls       protein-entry candidate rows between validation assemblies whose check was same_context or
                 only_one_side_different: how DNA rows at the same loci were judged (a DNA pass there is listed).
"""
import argparse
import gzip
import json
from collections import Counter, defaultdict
from pathlib import Path

import ucd_compare as uc

MIN_OVERLAP = 0.5


def iv(s):
    a, b = s.split('-')
    return int(a), int(b)


def recip(x, y):
    ov = min(x[1], y[1]) - max(x[0], y[0])
    return ov > 0 and ov >= MIN_OVERLAP * (x[1] - x[0]) and ov >= MIN_OVERLAP * (y[1] - y[0])


def stream(path, keep):
    with gzip.open(path, 'rt') as f:
        h = f.readline().rstrip('\n').split('\t')
        for line in f:
            r = dict(zip(h, line.rstrip('\n').split('\t')))
            if keep(r):
                yield r


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--manifest', required=True)
    p.add_argument('--dna-work', required=True)
    p.add_argument('--dna-collect', required=True)
    p.add_argument('--dna-rows', required=True, help='ucd_batch rows_candidate_insertion of the DNA pairs')
    p.add_argument('--protein-final', required=True)
    p.add_argument('--out', required=True)
    a = p.parse_args()
    S = {m['sample'] for m in uc.read_tsv(a.manifest)}
    inset = lambda r: r['a_sample'] in S and r['b_sample'] in S
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    # DNA rows (checked)
    dna = list(stream(a.dna_rows, lambda r: True))
    strict = lambda r: r['final_status'] == 'pass' and 0 < int(r['x_len_a'] or 0) < 2500 and 0 < int(r['x_len_b'] or 0) < 2500
    by_loc = defaultdict(list)
    for r in dna:
        k = tuple(sorted(((r['a_sample'], r['a_contig']), (r['b_sample'], r['b_contig']))))
        by_loc[k].append(r)
    D = dict(rows=len(dna), source_check=Counter(r['source_check'] for r in dna),
             ucd_check=Counter(r['ucd_check'] for r in dna), final=Counter(r['final_status'] for r in dna),
             strict_pass=sum(map(strict, dna)),
             strict_by_scope=Counter(r['pair_scope'] for r in dna if strict(r)),
             strict_by_strands=Counter('%s%s' % (r['a_strand'], r['b_strand']) for r in dna if strict(r)),
             source_check_failures=Counter(x for r in dna for x in r['source_reasons'].split(';') if x),
             screen_check_agree=Counter(r['screen_check_agree'] for r in dna))
    allp = list(stream(Path(a.dna_collect) / 'dna_pairs_all.tsv.gz', lambda r: True))
    D['candidate_pairs'] = len(allp)
    D['pair_status'] = Counter(r['status'] for r in allp)
    D['unresolved_ends'] = Counter('%s/%s' % (r['boundary_up'], r['boundary_down']) for r in allp
                                   if r['status'] != 'resolved')
    strict_occ = {(r[s + '_sample'], r[s + '_contig'], r[s + '_x_interval']) for r in dna if strict(r) for s in 'ab'}
    D['strict_pass_occurrences'] = len(strict_occ)
    pairs_loc = defaultdict(list)
    for r in allp:
        k = tuple(sorted(((r['a_sample'], r['a_contig']), (r['b_sample'], r['b_contig']))))
        pairs_loc[k].append(r)
    occ = defaultdict(list)
    for f in Path(a.dna_work).glob('occ/lib*/*/b*.tsv.gz'):
        for r in stream(f, lambda r: True):
            occ[(r['sample'], r['contig'])].append((int(r['lo']), int(r['hi']), r['library_id']))

    def side_match(r, s, x):
        return r[s + '_sample'] == x[0] and r[s + '_contig'] == x[1] and recip(iv(r[s + '_x_interval']), x[2])

    def classify(pr):
        A = (pr['a_sample'], pr['a_contig'], iv(pr['a_x_interval']))
        B = (pr['b_sample'], pr['b_contig'], iv(pr['b_x_interval']))
        k = tuple(sorted(((A[0], A[1]), (B[0], B[1]))))
        hits = [r for r in by_loc.get(k, []) if (side_match(r, 'a', A) and side_match(r, 'b', B)) or
                (side_match(r, 'a', B) and side_match(r, 'b', A))]
        if hits:
            if any(strict(r) for r in hits):
                return 'dna_strict_pass', hits
            if any(r['final_status'] == 'pass' for r in hits):
                return 'dna_pass_not_strict', hits
            return 'dna_checked_not_passed:' + sorted(r['ucd_check'] for r in hits)[0], hits

        def seed_ov(r, s, x):
            if r[s + '_sample'] != x[0] or r[s + '_contig'] != x[1]:
                return False
            y = iv(r[s + '_seed'])
            return min(y[1], x[2][1]) > max(y[0], x[2][0])
        un = [r for r in pairs_loc.get(k, []) if r['status'] != 'resolved' and
              ((seed_ov(r, 'a', A) and seed_ov(r, 'b', B)) or (seed_ov(r, 'a', B) and seed_ov(r, 'b', A)))]
        if un:
            return 'dna_pair_unresolved:%s/%s' % (un[0]['boundary_up'], un[0]['boundary_down']), []
        la = {l for lo, hi, l in occ.get((A[0], A[1]), []) if min(hi, A[2][1]) > max(lo, A[2][0])}
        lb = {l for lo, hi, l in occ.get((B[0], B[1]), []) if min(hi, B[2][1]) > max(lo, B[2][0])}
        if la & lb:
            return 'occurrences_not_paired', []
        if la or lb:
            return 'occurrence_one_side', []
        return 'not_found', []

    prot_strict = list(stream(Path(a.protein_final) / 'strict_lt2500_candidate_insertion.tsv.gz', inset))
    rec, rows = Counter(), []
    by_scope, by_strand = defaultdict(Counter), defaultdict(Counter)
    for pr in prot_strict:
        c, hits = classify(pr)
        top = c.split(':')[0]
        rec[c] += 1
        by_scope[pr['pair_scope']][top] += 1
        by_strand['%s%s' % (pr['a_strand'], pr['b_strand'])][top] += 1
        rows.append(dict(pair_id=pr['pair_id'], core_key=pr['core_key'], pair_scope=pr['pair_scope'],
                         a_sample=pr['a_sample'], a_contig=pr['a_contig'], a_x=pr['a_x_interval'],
                         a_strand=pr['a_strand'], b_sample=pr['b_sample'], b_contig=pr['b_contig'],
                         b_x=pr['b_x_interval'], b_strand=pr['b_strand'], dna_class=c,
                         dna_pair_ids=';'.join(sorted({r['pair_id'] for r in hits}))[:500]))
    uc.write_tsv(out / 'recovery.tsv.gz', rows, list(rows[0]) if rows else ['pair_id'])
    cores = defaultdict(set)
    for r in rows:
        cores[r['core_key']].add(r['dna_class'].split(':')[0])
    R = dict(protein_strict_pairs=len(prot_strict), classes=dict(rec.most_common()),
             by_scope={k: dict(v) for k, v in by_scope.items()}, by_strands={k: dict(v) for k, v in by_strand.items()},
             protein_cores=len(cores), cores_with_dna_strict_pass=sum('dna_strict_pass' in v for v in cores.values()),
             cores_with_any_dna_check=sum(bool(v & {'dna_strict_pass', 'dna_pass_not_strict', 'dna_checked_not_passed'})
                                          for v in cores.values()))

    ctl = Counter()
    for pr in stream(Path(a.protein_final) / 'all_rows_candidate_insertion.tsv.gz',
                     lambda r: inset(r) and r['ucd_check'] in ('same_context', 'only_one_side_different')):
        c, hits = classify(pr)
        ctl['%s -> %s' % (pr['ucd_check'], c)] += 1
    C = dict(protein_control_rows=sum(ctl.values()), dna_judgement=dict(ctl.most_common()))
    res = dict(assemblies=len(S), min_overlap=MIN_OVERLAP, dna=D, recovery=R, controls=C)
    res = json.loads(json.dumps(res, default=lambda o: dict(o)))
    (out / 'validation.json').write_text(json.dumps(res, indent=1) + '\n')
    print(json.dumps(res, indent=1))


if __name__ == '__main__':
    main()
