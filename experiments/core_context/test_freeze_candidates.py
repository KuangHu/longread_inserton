import unittest

import freeze_candidates as fc


def pair(cls='context_different_both_sides', st='single_gene'):
    return dict(structure=st, context_class=cls, occurrence_a='a', occurrence_b='b')


class EligibleTests(unittest.TestCase):
    def test_gene_order_rules_only(self):
        occ = {'a': dict(s6_status='site_group'), 'b': dict(s6_status='single_observation_site')}
        self.assertTrue(fc.eligible(pair(), occ, 0))
        for cls in ('context_repeated_both_sides', 'context_repeated_one_side', 'context_unknown'):
            self.assertFalse(fc.eligible(pair(cls), occ, 0))
        self.assertFalse(fc.eligible(pair(st='two_gene'), occ, 0))
        # interlocking windows of two nearby copies: different labels per side, but the same physical genes
        self.assertFalse(fc.eligible(pair(), occ, 4))
        occ['b']['s6_status'] = 'unresolved'
        self.assertFalse(fc.eligible(pair(), occ, 0))


class RankTests(unittest.TestCase):
    def row(self, res, dna, copies, key='k'):
        return dict(copies_both_ends_resolved=res, dna_v2_not_shared_both=dna, copies_in_assembly=copies, subset='s',
                    core_key=key, occurrence_a='a', occurrence_b='b')

    def test_order(self):
        rows = [self.row(0, True, 9, 'z'), self.row(2, False, 2, 'y'), self.row(2, True, 2, 'x'), self.row(1, True, 5, 'w')]
        self.assertEqual([r['core_key'] for r in sorted(rows, key=fc.rank)], ['x', 'y', 'w', 'z'])

    def test_unknown_boundaries_stay_eligible(self):
        # an unresolved boundary only ranks lower; it never removes the candidate
        self.assertLess(fc.rank(self.row(0, True, 2)), fc.rank(self.row(0, False, 2)))


if __name__ == '__main__':
    unittest.main()
