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

    def test_count_mode(self):
        from collections import Counter
        unit_of = {'r1': 'molecules', 'r2': 'read_records'}
        self.assertEqual(u.count_mode(Counter({('molecules', 'P_matched'): 5, ('read_records', 'off_locus'): 2}), unit_of),
                         'molecules')
        self.assertEqual(u.count_mode(Counter({('molecules', 'P_matched'): 5, ('read_records', 'A_matched'): 1}), unit_of),
                         'mixed')
        self.assertEqual(u.count_mode(Counter(), {'r1': 'assumed_molecules'}), 'assumed_molecules')

    def test_occurrence_class(self):
        self.assertEqual(u.occurrence_class(['supported', 'supported']), 'all_intervals_supported')
        self.assertEqual(u.occurrence_class(['supported', 'locus_ambiguous']), 'some_intervals_supported')
        self.assertEqual(u.occurrence_class(['supported', 'conflict']), 'conflict_at_an_interval')


class CompletenessTests(unittest.TestCase):
    positions = {'P_1': dict(sample='S1'), 'P_2': dict(sample='S1'), 'P_3': dict(sample='S2')}
    runs_used = [dict(sample='S1', run='R1', status='checked'), dict(sample='S2', run='R2', status='checked')]
    ready = {'S1': ['R1'], 'S2': ['R2']}

    def res(self, *ids):
        return {t: dict(test_id=t, sample=self.positions.get(t, dict(sample='S1'))['sample'], verdict='consistent')
                for t in ids}

    def test_complete_by_ids(self):
        c = u.completeness(self.positions, self.res('P_1', 'P_2', 'P_3'), self.runs_used, self.ready)
        self.assertTrue(c['complete'])
        self.assertEqual((c['assemblies_complete'], c['runs_checked'], c['intervals_with_result']), (2, 2, 3))

    def test_same_row_count_with_other_ids_is_not_complete(self):
        c = u.completeness(self.positions, self.res('P_1', 'P_9', 'P_3'), self.runs_used, self.ready)
        self.assertFalse(c['complete'])
        self.assertEqual(c['missing_interval_ids'], {'S1': ['P_2']})
        self.assertEqual(c['missing_unexpected_interval_ids'], {'S1': ['P_9']})

    def test_unchecked_ready_run_is_not_complete(self):
        c = u.completeness(self.positions, self.res('P_1', 'P_2', 'P_3'), self.runs_used, dict(self.ready, S2=['R2', 'R3']))
        self.assertFalse(c['complete'])
        self.assertEqual(c['missing_runs'], ['S2:R3'])
        self.assertEqual(c['assemblies_complete'], 1)

    def test_progress_snapshot_counts_checked_only(self):
        res = self.res('P_1', 'P_3')
        snap = u.progress_snapshot(u.completeness(self.positions, res, self.runs_used[:1], self.ready), self.positions, res)
        self.assertIs(snap['complete'], False)
        self.assertEqual(snap['checked_of_expected'], dict(assemblies='0/2', runs='1/2', intervals='2/3'))
        self.assertEqual(snap['intervals_not_checked'], 1)
        self.assertEqual(snap['interval_class_among_checked'], {'supported': 2})
        self.assertNotIn('interval_class', snap)


if __name__ == '__main__':
    unittest.main()
