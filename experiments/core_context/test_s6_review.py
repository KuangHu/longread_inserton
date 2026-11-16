import unittest

import s6_review as rv


def pr(verdict, up='different', down='different', a='S1', b='S2'):
    return dict(verdict=verdict, up_side=up, down_side=down, sample_a=a, sample_b=b)


class ContextClassTests(unittest.TestCase):
    def test_coordinates_and_context_are_separate(self):
        self.assertEqual(rv.context_class('different', 'different'), 'context_different_both_sides')
        self.assertEqual(rv.context_class('same', 'same'), 'context_repeated_both_sides')
        self.assertEqual(rv.context_class('same', 'different'), 'context_repeated_one_side')
        self.assertEqual(rv.context_class('different', 'unknown'), 'context_unknown')
        self.assertEqual(rv.context_class('same', 'unknown'), 'context_unknown')


class UnitRelationTests(unittest.TestCase):
    def test_all_cross_pairs_different(self):
        cross = [pr('different_position_candidate'), pr('different_position_same_assembly', a='S1', b='S1')]
        self.assertEqual(rv.unit_relation(cross), ('different', 'context_different_both_sides'))

    def test_repeated_context_is_still_different_coordinates(self):
        cross = [pr('different_position_candidate'), pr('different_position_same_assembly', 'same', 'same', 'S1', 'S1')]
        self.assertEqual(rv.unit_relation(cross), ('different', 'context_repeated_in_some_pair'))
        cross = [pr('different_position_same_assembly', 'unknown', 'different', 'S1', 'S1')]
        self.assertEqual(rv.unit_relation(cross), ('different', 'context_unknown_in_some_pair'))

    def test_two_groups_are_not_different_by_default(self):
        # one unresolved cross pair: two accepted groups, relation between them unknown
        self.assertEqual(rv.unit_relation([pr('different_position_candidate'), pr('unresolved', 'unknown')])[0], 'unknown')
        self.assertEqual(rv.unit_relation([])[0], 'unknown')

    def test_ambiguous_and_conflicting(self):
        self.assertEqual(rv.unit_relation([pr('different_position_candidate'), pr('same_position_candidate', 'same', 'same')])[0],
                         'ambiguous')
        self.assertEqual(rv.unit_relation([pr('conflicting_sides', 'same', 'different')])[0], 'conflicting')


class CliqueTests(unittest.TestCase):
    def test_largest_pairwise_different_set(self):
        # A-B, B-C, A-C different, D different only from A: 3, not 4
        adj = {'A': {'B', 'C', 'D'}, 'B': {'A', 'C'}, 'C': {'A', 'B'}, 'D': {'A'}}
        self.assertEqual(rv.max_clique(list(adj), adj), 3)
        # A-B and B-C different, A-C unknown: only 2 units are shown pairwise different
        adj = {'A': {'B'}, 'B': {'A', 'C'}, 'C': {'B'}}
        self.assertEqual(rv.max_clique(list(adj), adj), 2)
        self.assertEqual(rv.max_clique(['A', 'B'], {'A': set(), 'B': set()}), 1)


class FlankTests(unittest.TestCase):
    def test_covered_union(self):
        self.assertEqual(rv.covered([(0, 10), (5, 20), (30, 40)], 0, 35), 25)
        self.assertEqual(rv.covered([(0, 10)], 20, 30), 0)

    def test_outer_region_skips_the_core_end(self):
        self.assertEqual(rv.outer('up', 3000), (0, 2500))
        self.assertEqual(rv.outer('down', 3000), (500, 3000))
        self.assertEqual(rv.outer('up', 300), (0, 0))

    def test_flank_class(self):
        self.assertEqual(rv.flank_class(0.9, 2500), 'shared')
        self.assertEqual(rv.flank_class(0.0, 2500), 'not_shared')
        self.assertEqual(rv.flank_class(0.3, 2500), 'partial')
        self.assertEqual(rv.flank_class(0.0, 400), 'short_flank')

    def test_flank_check_blastn(self):
        import random
        rnd = random.Random(1)
        seq = lambda n: ''.join(rnd.choice('ACGT') for _ in range(n))
        shared_up, end = seq(2500), seq(500)
        # up: the outer 2,500 bp are identical; down: only the 500 bp next to the core (element end) is shared
        x = dict(up=shared_up + end, down=end + seq(2500))
        y = dict(up=shared_up + seq(500), down=end + seq(2500))
        fc = rv.flank_check(x, y)
        self.assertEqual(fc['up'][2], 'shared')
        self.assertEqual(fc['down'][2], 'not_shared')


if __name__ == '__main__':
    unittest.main()
