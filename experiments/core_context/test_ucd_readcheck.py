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


class RotatedCompletenessTests(unittest.TestCase):
    def test_rotated_runs_expected_and_results_required(self):
        positions = {'T1': dict(sample='S1'), 'T2': dict(sample='S1')}
        rot = [dict(test_id='T2', sample='S1', pseudo_sample='S1@rot')]
        runs = [dict(sample='S1', run='R1', status='checked')]
        res = {'T1': dict(sample='S1'), 'T2': dict(sample='S1', verdict='not_tested:element_interval_crosses_contig_end')}
        c = u.completeness(positions, res, runs, {'S1': ['R1']}, rot)
        self.assertFalse(c['complete'])
        self.assertEqual(c['missing_rotated_interval_ids'], ['T2'])
        self.assertEqual(c['missing_runs'], ['S1@rot:R1'])
        res['T2'] = dict(sample='S1', verdict='consistent', coordinates='rotated:complete_in_header:offset=500')
        c = u.completeness(positions, res, runs + [dict(sample='S1@rot', run='R1', status='checked')], {'S1': ['R1']}, rot)
        self.assertTrue(c['complete'])
        self.assertEqual((c['rotated_runs_checked'], c['rotated_intervals_with_result']), (1, 1))
        c = u.completeness(positions, res, runs + [dict(sample='S1@rot', run='R1', status='checked')], {'S1': ['R1']})
        self.assertEqual(c['missing_unexpected_runs'], ['S1@rot:R1'])


class RotateTests(unittest.TestCase):
    def test_interval_across_the_start_is_moved_inside_the_contig(self):
        import random
        import tempfile
        from pathlib import Path
        from unittest import mock
        random.seed(1)
        q = ''.join(random.choice('ACGT') for _ in range(1000))
        with tempfile.TemporaryDirectory() as t:
            t = Path(t)
            fa = t / 'S1.fna'
            fa.write_text('>c1 Some bacterium chromosome, complete genome\n%s\n>c2 plasmid\nACGTACGT\n' % q)
            rd = t / 'round'
            (rd / 'check' / 'S1').mkdir(parents=True)
            u.pc.write_tsv(rd / 'positions.tsv', [dict(test_id='T1', sample='S1', contig='c1', lo=-40, hi=60),
                                                 dict(test_id='T2', sample='S1', contig='c1', lo=950, hi=1030),
                                                 dict(test_id='T3', sample='S1', contig='c1', lo=100, hi=200)],
                           ['test_id', 'sample', 'contig', 'lo', 'hi'])
            u.pc.write_tsv(rd / 'check' / 'S1' / 'positions_result.tsv',
                           [dict(test_id='T1', verdict='not_tested:element_interval_crosses_contig_end'),
                            dict(test_id='T2', verdict='not_tested:element_interval_crosses_contig_end'),
                            dict(test_id='T3', verdict='consistent')], ['test_id', 'verdict'])
            real = u.fasta_of
            with mock.patch.object(u, 'PROJ', t), \
                    mock.patch.object(u, 'fasta_of', lambda o, s: str(fa) if s == 'S1' else real(o, s)), \
                    mock.patch('sys.stdout'):
                u.cmd_rotate(mock.Mock(round_dir=rd, repo_out=None))
                rot = u.rc3.read_fasta(real(None, 'S1' + u.ROT))
            rows = {r['test_id']: r for r in u.pc.read_tsv(rd / 'rotated_positions.tsv')}
            self.assertEqual(sorted(rows), ['T1', 'T2'])
            self.assertEqual(rows['T1']['header_annotation'], 'complete_in_header')
            self.assertEqual(rot['c2'], 'ACGTACGT')
            wrap = {'T1': q[-40:] + q[:60], 'T2': q[950:] + q[:30]}
            for k, r in rows.items():
                self.assertEqual(rot['c1'][int(r['rot_lo']):int(r['rot_hi'])], wrap[k])
            self.assertEqual((rd / 'rotated_samples.txt').read_text(), 'S1@rot\n')


if __name__ == '__main__':
    unittest.main()
