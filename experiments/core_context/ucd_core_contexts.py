"""Distinct-context count per core (summary after the two checks; discovery and checks unchanged).

Counted are distinct contexts, certified by pairwise passes, never pairs or raw occurrences:
  position   one physical site: X intervals of the checked rows on the same sample and contig whose reciprocal
             overlap is >= POS_OVERLAP are one position (several seeds, library sequences or boundary versions of one
             site count once).
  edge       two positions are linked when at least one row between them passes both checks with
             TARGET_MIN <= X < TARGET_MAX on both sides (the target table). Position pairs that have both a passing
             and a non-passing row (another boundary version) are counted and reported.
  core       a connected component of the edge graph (X / X' similar along every edge).
  contexts   the largest set of positions of a core that pass PAIRWISE (a clique): k representative positions with
             all k(k-1)/2 pairs passing. This prevents "A differs from B, C, D but B, C, D share one context" (B-C,
             B-D, C-D would fail as same_context and break the clique) and related strains with X between the same
             U / D (their pairs fail). Searched exactly up to CAP (reported as >= CAP beyond).
  variants   contexts_any (any positions), contexts_single_assembly (all positions in one assembly; copies within a
             genome) and contexts_distinct_assemblies (every position from a different assembly), reported apart.
Tiers (project rule, not a biological standard): 2-3 exploratory, >=4 main list, >=6 priority review, >=8 showcase.
A core below a tier whose library sequences still have deferred (not yet enumerated / aligned) pairs is
`insufficient_evidence`, not `below_tier`: unfinished pairs are not evidence of fewer contexts.
"""
import argparse
import glob
import gzip
import hashlib
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

TARGET_MIN, TARGET_MAX = 700, 2500
POS_OVERLAP = 0.5
CAP = 12
TIERS = (2, 4, 6, 8)
sys.setrecursionlimit(10000)


def stream(p):
    with gzip.open(p, 'rt') as f:
        h = f.readline().rstrip('\n').split('\t')
        for line in f:
            yield dict(zip(h, line.rstrip('\n').split('\t')))


def iv(s):
    a, b = s.split('-')
    return int(a), int(b)


def max_clique(nodes, adj, cap, distinct=None):
    """Largest clique (up to cap) among nodes: branch and bound with a greedy-colouring bound (Tomita style).
    distinct: node -> label; nodes with equal labels are treated as not adjacent (one per label in the clique)."""
    S = set(nodes)
    if distinct is None:
        nb = {v: adj[v] & S for v in S}
    else:
        nb = {v: {u for u in adj[v] & S if distinct[u] != distinct[v]} for v in S}
    best = []

    def colour_sort(P):
        classes = []
        for v in sorted(P, key=lambda x: -len(nb[x])):
            for c in classes:
                if not (nb[v] & c):
                    c.add(v)
                    break
            else:
                classes.append({v})
        order, cols = [], []
        for k, c in enumerate(classes, 1):
            for v in c:
                order.append(v)
                cols.append(k)
        return order, cols

    def expand(R, P):
        nonlocal best
        order, cols = colour_sort(P)
        for i in range(len(order) - 1, -1, -1):
            if len(R) + cols[i] <= len(best):
                return False
            v = order[i]
            R2 = R + [v]
            NP = P & nb[v]
            if NP:
                if expand(R2, NP):
                    return True
            elif len(R2) > len(best):
                best = R2
                if len(best) >= cap:
                    return True
            P = P - {v}
        return False

    if S:
        best = [next(iter(S))]
        expand([], set(S))
    return best[:cap]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--final', required=True, help='dir with all_rows_ / strict_700_2500_candidate_insertion.tsv.gz')
    ap.add_argument('--dna-work', default='', help='dna_discovery work dir (deferred chunks -> insufficient_evidence)')
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    fin, out = Path(a.final), Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    # positions from every checked row (so boundary versions of a site merge even if only one passes)
    sides = defaultdict(set)
    for r in stream(fin / 'all_rows_candidate_insertion.tsv.gz'):
        for s in 'ab':
            if r[s + '_x_interval']:
                sides[(r[s + '_sample'], r[s + '_contig'])].add(iv(r[s + '_x_interval']))
    pos_of, pos_info = {}, {}
    for (smp, ctg), ivs in sides.items():
        ivs = sorted(ivs)
        groups = []  # greedy merge in start order: reciprocal overlap with the group's first interval
        for x in ivs:
            for g in groups:
                y = g[0]
                ov = min(x[1], y[1]) - max(x[0], y[0])
                if ov > 0 and ov >= POS_OVERLAP * (x[1] - x[0]) and ov >= POS_OVERLAP * (y[1] - y[0]):
                    g.append(x)
                    break
            else:
                groups.append([x])
        for g in groups:
            pid = '%s|%s|%d-%d' % (smp, ctg, g[0][0], g[0][1])
            pos_info[pid] = dict(sample=smp, contig=ctg, intervals=len(g), lo=min(x[0] for x in g),
                                 hi=max(x[1] for x in g))
            for x in g:
                pos_of[(smp, ctg, x)] = pid

    def P(r, s):
        return pos_of[(r[s + '_sample'], r[s + '_contig'], iv(r[s + '_x_interval']))]

    def in_target(r):
        return all(TARGET_MIN <= int(r['x_len_' + s] or 0) < TARGET_MAX for s in 'ab')

    adj = defaultdict(set)
    edge_rows, fail_pairs, libs_of = defaultdict(list), set(), defaultdict(set)
    for r in stream(fin / 'all_rows_candidate_insertion.tsv.gz'):
        if not r['a_x_interval'] or not r['b_x_interval']:
            continue
        p, q = P(r, 'a'), P(r, 'b')
        if p == q:
            continue
        k = (min(p, q), max(p, q))
        if r['final_status'] == 'pass' and in_target(r):
            adj[p].add(q)
            adj[q].add(p)
            edge_rows[k].append(r['pair_id'])
            libs_of[p].add(r['core_key'])
            libs_of[q].add(r['core_key'])
        else:
            fail_pairs.add(k)
    conflict = sum(1 for k in edge_rows if k in fail_pairs)

    # deferred library sequences (not yet enumerated / aligned chunks)
    pending = set()
    if a.dna_work:
        done = set()
        for f in glob.glob(str(Path(a.dna_work) / 'enum' / 'g*' / 'b*' / 'progress.tsv.gz')):
            for r in stream(f):
                k = (r['library_id'], r['chunk'], r['occurrences'])
                if r['status'] == 'enumerated':
                    done.add(k)
                else:
                    pending.add(k)
        pending = {k[0] for k in pending - done}

    # components
    seen, comps = set(), []
    for n in sorted(adj):
        if n in seen:
            continue
        stack, comp = [n], []
        seen.add(n)
        while stack:
            v = stack.pop()
            comp.append(v)
            for u in adj[v]:
                if u not in seen:
                    seen.add(u)
                    stack.append(u)
        comps.append(sorted(comp))

    rows, summ = [], defaultdict(Counter)
    for comp in comps:
        smp = {n: pos_info[n]['sample'] for n in comp}
        best_any = max_clique(comp, adj, CAP)
        best_dist = max_clique(comp, adj, CAP, distinct=smp)
        by_s = defaultdict(list)
        for n in comp:
            by_s[smp[n]].append(n)
        best_single, single_s = [], ''
        for s_, ns in by_s.items():
            if len(ns) > len(best_single):
                c = max_clique(ns, adj, CAP)
                if len(c) > len(best_single):
                    best_single, single_s = c, s_
        n_any, n_single, n_dist = len(best_any), len(best_single), len(best_dist)
        libs = set().union(*(libs_of[n] for n in comp))
        unfinished = bool(libs & pending)
        tier = 'showcase_ge8' if n_any >= 8 else 'priority_ge6' if n_any >= 6 else 'main_ge4' if n_any >= 4 else \
            'exploratory_2_3'
        status = 'final' if not unfinished else 'insufficient_evidence'
        reps = best_any[:8]
        rep_pairs = []
        for i in range(min(4, len(reps))):
            for j in range(i + 1, min(4, len(reps))):
                k = (min(reps[i], reps[j]), max(reps[i], reps[j]))
                rep_pairs.append(sorted(edge_rows[k])[0])
        cid = 'K' + hashlib.sha1('|'.join(comp).encode()).hexdigest()[:12]
        rows.append(dict(core_id=cid, positions=len(comp), assemblies=len(set(smp.values())),
                         passing_position_pairs=sum(len(adj[n] & set(comp)) for n in comp) // 2,
                         contexts_any=('>=%d' % CAP) if n_any >= CAP else n_any,
                         contexts_single_assembly=('>=%d' % CAP) if n_single >= CAP else n_single,
                         single_assembly_sample=single_s,
                         contexts_distinct_assemblies=('>=%d' % CAP) if n_dist >= CAP else n_dist,
                         tier=tier, evidence_status=status, library_ids=';'.join(sorted(libs))[:300],
                         representative_positions=';'.join(reps),
                         representative_pair_ids_first4=';'.join(rep_pairs)))
        for name, n in (('any', n_any), ('single_assembly', n_single), ('distinct_assemblies', n_dist)):
            for t in TIERS:
                if n >= t:
                    summ[name]['ge%d' % t] += 1
                elif unfinished:
                    summ[name]['insufficient_evidence_below_%d' % t] += 1
    cols = list(rows[0]) if rows else ['core_id']
    import ucd_compare as uc
    uc.write_tsv(out / 'core_contexts.tsv.gz', rows, cols)
    uc.write_tsv(out / 'core_main_ge4.tsv.gz', [r for r in rows if r['tier'] != 'exploratory_2_3'], cols)
    res = dict(rules=dict(target='%d <= X < %d bp both sides, pass both checks' % (TARGET_MIN, TARGET_MAX),
                          position_overlap=POS_OVERLAP, contexts='max clique of pairwise-passing positions',
                          cap=CAP, tiers='2-3 exploratory, >=4 main, >=6 priority, >=8 showcase'),
               positions_total=len(pos_info), positions_in_pass_graph=len(adj), cores=len(rows),
               position_pairs_with_pass_and_fail_versions=conflict,
               cores_by_tier=dict(Counter(r['tier'] for r in rows)),
               cores_by_evidence=dict(Counter(r['evidence_status'] for r in rows)),
               cores_ge=dict((k, dict(v)) for k, v in summ.items()),
               pending_library_sequences=len(pending))
    (out / 'core_contexts_summary.json').write_text(json.dumps(res, indent=1) + '\n')
    print(json.dumps(res, indent=1))


if __name__ == '__main__':
    main()
