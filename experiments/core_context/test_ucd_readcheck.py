import unittest

import ucd_readcheck as u
from ucd_readcheck import rc3


class GeometryTests(unittest.TestCase):
    def test_element_follows_core_orientation(self):
        o = dict(start0='1000', end0='1600', strand='+')
        self.assertEqual(u.element(o, 50, 20), (950, 1620))
        o['strand'] = '-'  # upstream in core orientation is to the right on the contig
        self.assertEqual(u.element(o, 50, 20), (980, 1650))

    def test_alleles_and_anchor_positions_match_the_assembly(self):
        seq = ''.join('ACGT'[(i * 7 + i // 3) % 4] for i in range(9000))
        p = dict(test_id='P_x', contig='c', lo='3000', hi='4200')
        tests, fa, skipped = u.make_tests([p], {'c': seq})
        self.assertFalse(skipped)
        al = dict(fa)
        self.assertEqual(al['P_x|C|P'], seq[1000:6200])
        self.assertEqual(al['P_x|C|A'], seq[1000:3000] + seq[4200:6200])
        g = rc3.test_geometry(tests[0])['C']
        for allele in 'PA':
            lo, hi, _ = g['iv'][allele]
            el, er = g['expect'][allele]
            # allele position lo / hi - 1 sit at these assembly positions
            self.assertEqual(al['P_x|C|' + allele][lo], seq[el])
            self.assertEqual(al['P_x|C|' + allele][hi - 1], seq[er])
        self.assertEqual(g['limit'], 600)
        self.assertEqual(g['need'], 1200 + 2 * rc3.ANCHOR)

    def test_interval_crossing_contig_end_is_not_tested(self):
        tests, fa, skipped = u.make_tests([dict(test_id='P_y', contig='c', lo='-5', hi='300')], {'c': 'A' * 1000})
        self.assertEqual((tests, fa), ([], []))
        self.assertEqual(skipped[0]['reason'], 'element_interval_crosses_contig_end')


class ClassTests(unittest.TestCase):
    def test_tiers(self):
        self.assertEqual(u.tier_of('same_assembly', True, True), 'A')
        self.assertIsNone(u.tier_of('same_assembly', False, False))
        self.assertEqual(u.tier_of('cross_assembly', True, True), 'B')
        self.assertEqual(u.tier_of('cross_assembly', False, True), 'C')
        self.assertIsNone(u.tier_of('cross_assembly', False, False))

    def test_position_class(self):
        row = dict(verdict='insufficient', locus_ambiguous='2')
        self.assertEqual(u.position_class(row), 'locus_ambiguous')
        self.assertEqual(u.position_class(dict(row, locus_ambiguous='0')), 'insufficient_coverage')
        self.assertEqual(u.position_class(dict(verdict='mixed')), 'conflict')
        self.assertEqual(u.position_class(dict(verdict='consistent')), 'supported')

    def test_pair_class(self):
        self.assertEqual(u.pair_class('A', 'supported', 'supported'), 'both_positions_supported')
        self.assertEqual(u.pair_class('B', 'supported', 'locus_ambiguous'), 'one_position_supported')
        self.assertEqual(u.pair_class('A', 'supported', 'conflict'), 'conflict_at_a_position')
        self.assertEqual(u.pair_class('C', 'supported', None), 'single_side_supported')
        self.assertEqual(u.pair_class('C', None, 'insufficient_coverage'), 'single_side_not_supported')


if __name__ == '__main__':
    unittest.main()
