import itertools
import random
import unittest
from collections import Counter

import protein_context as pc
import scale_compare as sc


def brute(occs):
    """Every occurrence pair with protein_context.relation (the frozen subset rule)."""
    tot, same = Counter(), Counter()
    for (a, sa, ua, da), (b, sb, ub, db) in itertools.combinations(occs, 2):
        r = pc.relation(dict(up=[(l, True, '') for l in ua], down=[(l, True, '') for l in da]),
                        dict(up=[(l, True, '') for l in ub], down=[(l, True, '') for l in db]))
        tot[r] += 1; same[r] += sa == sb
    return tot, same


def run(occs, limit=sc.BLOCK_PAIRS, budget=None):
    sigs, mult, index, _ = sc.collapse(occs)
    recs, used = [], 0
    for i0, i1 in sc.blocks(mult, limit):
        if budget is not None and used >= budget:
            break
        rec = sc.run_block(sigs, mult, i0, i1)
        recs.append(rec); used += i1 - i0
    return sc.finalize(sigs, mult, recs), index


def random_occs(n, labels, samples, rnd, missing=0.05, short=0.05):
    out = []
    for k in range(n):
        side = lambda: [None if rnd.random() < missing else rnd.choice(labels)
                        for _ in range(2 if rnd.random() < short else 3)]
        out.append(('o%d' % k, rnd.choice(samples), side(), side()))
    return out


class EquivalenceTests(unittest.TestCase):
    def test_counts_equal_old_rule_on_every_pair(self):
        rnd = random.Random(3)
        for _ in range(20):
            occs = random_occs(rnd.randint(2, 60), ['L%d' % i for i in range(rnd.randint(2, 12))],
                               ['S%d' % i for i in range(rnd.randint(1, 6))], rnd)
            tot, same = brute(occs)
            (summary, _), _ = run(occs, limit=rnd.randint(1, 50))
            self.assertTrue(summary['complete'])
            for r in sc.RELATIONS:
                self.assertEqual(summary[r], tot[r], r)
                self.assertEqual(summary[r + '_same_assembly'], same[r], r)

    def test_above_old_cap_is_complete(self):
        # 260 occurrences = 33,670 pairs > MAX_PAIRS (20,000); the old compare kept only the first 20,000 in input order
        rnd = random.Random(5)
        occs = random_occs(260, ['L%d' % i for i in range(40)], ['S%d' % i for i in range(30)], rnd)
        self.assertGreater(260 * 259 // 2, pc.MAX_PAIRS)
        (summary, _), _ = run(occs, limit=500)
        self.assertEqual((summary['pairs_total'], summary['pairs_accounted'], summary['pairs_not_compared']),
                         (33670, 33670, 0))
        tot, _ = brute(occs)
        self.assertEqual({r: summary[r] for r in sc.RELATIONS}, {r: tot[r] for r in sc.RELATIONS})


class DiscoveryTests(unittest.TestCase):
    def test_single_copy_per_assembly_cross_assembly(self):
        # one copy in each of two assemblies, different contexts: must be found without any multi-copy assembly
        occs = [('a', 'S1', ['A', 'B', 'C'], ['D', 'E', 'F']), ('b', 'S2', ['P', 'Q', 'R'], ['X', 'Y', 'Z'])]
        (summary, rows), index = run(occs)
        self.assertEqual(summary['different_context_candidate'], 1)
        self.assertEqual(summary['different_context_candidate_same_assembly'], 0)
        got = sc.name_partners([('a', 'S1'), ('b', 'S2')], index, rows)['a']
        self.assertEqual((got['diff_partners_cross'], got['diff_partner_cross'], got['diff_partner_same']), (1, 'b', ''))

    def test_two_gene_core_contexts_are_compared_the_same_way(self):
        # two-gene occurrences carry their own up/down contexts; the comparator is structure-agnostic
        occs = [('t|A,B|0', 'S1', ['U1', 'U2', 'U3'], ['D1', 'D2', 'D3']),
                ('t|A,B|1', 'S1', ['V1', 'V2', 'V3'], ['E1', 'E2', 'E3']),
                ('t|A,B|2', 'S2', ['U1', 'W2', 'W3'], ['D1', 'F2', 'F3'])]
        (summary, _), _ = run(occs)
        self.assertEqual((summary['different_context_candidate'], summary['different_context_candidate_same_assembly'],
                          summary['similar_context']), (2, 1, 1))

    def test_identical_signatures_collapse_but_pairs_are_counted(self):
        occs = [('o%d' % i, 'S%d' % (i % 3), ['A', 'B', 'C'], ['D', 'E', 'F']) for i in range(9)]
        (summary, _), _ = run(occs)
        self.assertEqual(summary['signatures'], 1)
        self.assertEqual((summary['similar_context'], summary['similar_context_same_assembly']), (36, 9))


def brute_partners(occs):
    """Per occurrence: the sets of cross-assembly and same-assembly different_context_candidate partners."""
    out = {o[0]: (set(), set()) for o in occs}
    for (a, sa, ua, da), (b, sb, ub, db) in itertools.combinations(occs, 2):
        r = pc.relation(dict(up=[(l, True, '') for l in ua], down=[(l, True, '') for l in da]),
                        dict(up=[(l, True, '') for l in ub], down=[(l, True, '') for l in db]))
        if r == 'different_context_candidate':
            out[a][sa == sb].add(b); out[b][sa == sb].add(a)
    return out


def named(occs, limit=sc.BLOCK_PAIRS):
    (summary, rows), index = run(occs, limit=limit)
    return sc.name_partners([(o[0], o[1]) for o in occs], index, rows)


class PartnerTests(unittest.TestCase):
    def test_review_counterexample(self):
        # x1 (S1) and x2 (S2) share one signature X; Y holds only S2, Z only S1; both differ from X. v1 kept a single
        # cross and a single same-assembly partner signature for all of X (the smallest one, here Z), so x1 got no
        # cross partner (Z holds only its own assembly) and x2 no same-assembly partner, although both exist.
        X = (['A', 'B', 'C'], ['D', 'E', 'F'])
        occs = [('x1', 'S1') + X, ('x2', 'S2') + X,
                ('y', 'S2', ['P', 'Q', 'R'], ['U', 'V', 'W']), ('z', 'S1', ['K', 'L', 'M'], ['N', 'O', 'T'])]
        got = named(occs)
        self.assertEqual((got['x2']['diff_partners_cross'], got['x2']['diff_partners_same_assembly']), (1, 1))
        self.assertEqual((got['x2']['diff_partner_cross'], got['x2']['diff_partner_same']), ('z', 'y'))
        self.assertEqual((got['x1']['diff_partner_cross'], got['x1']['diff_partner_same']), ('y', 'z'))

    def test_every_occurrence_matches_brute_force(self):
        rnd = random.Random(21)
        for _ in range(40):
            occs = random_occs(rnd.randint(2, 50), ['L%d' % i for i in range(rnd.randint(2, 9))],
                               ['S%d' % i for i in range(rnd.randint(1, 5))], rnd, missing=0.02, short=0.02)
            want = brute_partners(occs)
            got = named(occs, limit=rnd.randint(1, 30))
            for oid, (cross, same) in want.items():
                g = got[oid]
                self.assertEqual((g['diff_partners_cross'], g['diff_partners_same_assembly']), (len(cross), len(same)))
                self.assertEqual(g['diff_partner_cross'] != '', bool(cross))
                self.assertEqual(g['diff_partner_same'] != '', bool(same))
                if cross:
                    self.assertIn(g['diff_partner_cross'], cross)
                if same:
                    self.assertIn(g['diff_partner_same'], same)


class ResumeTests(unittest.TestCase):
    def test_unfinished_blocks_are_not_compared_and_resume_matches(self):
        rnd = random.Random(9)
        occs = random_occs(120, ['L%d' % i for i in range(25)], ['S%d' % i for i in range(8)], rnd)
        (part, _), _ = run(occs, limit=200, budget=10)
        self.assertFalse(part['complete'])
        self.assertGreater(part['pairs_not_compared'], 0)
        self.assertEqual(part['pairs_accounted'] + part['pairs_not_compared'], part['pairs_total'])
        self.assertEqual(part['context_classes_provisional'], '')
        (full, _), _ = run(occs, limit=200)
        (one, _), _ = run(occs)
        self.assertTrue(full['complete'])
        self.assertEqual(full, one)


if __name__ == '__main__':
    unittest.main()
