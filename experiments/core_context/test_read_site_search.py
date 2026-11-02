import unittest

import read_site_search as rs


def aln(t='c', ts=0, te=1000, strand='+', qs=0, qe=1000, qlen=1000, score=1000):
    return dict(t=t, ts=ts, te=te, strand=strand, qs=qs, qe=qe, qlen=qlen, score=score)


class PlacementTests(unittest.TestCase):
    def test_inner_end_is_next_to_the_core(self):
        self.assertEqual(rs.placement([aln(ts=100, te=1100)], 'L'), ('placed', 'c', 1100, '+'))
        self.assertEqual(rs.placement([aln(ts=2000, te=3000)], 'R'), ('placed', 'c', 2000, '+'))
        self.assertEqual(rs.placement([aln(ts=100, te=1100, strand='-')], 'L'), ('placed', 'c', 100, '-'))
        self.assertEqual(rs.placement([aln(ts=2000, te=3000, strand='-')], 'R'), ('placed', 'c', 3000, '-'))

    def test_competing_placement_is_ambiguous_and_short_alignment_unplaced(self):
        two = [aln(ts=100, te=1100, score=1000), aln(t='d', ts=5000, te=6000, score=990)]
        self.assertEqual(rs.placement(two, 'L')[0], 'ambiguous')
        far = [aln(ts=100, te=1100, score=1000), aln(t='d', ts=5000, te=6000, score=900)]
        self.assertEqual(rs.placement(far, 'L')[0], 'placed')
        self.assertEqual(rs.placement([aln(qs=0, qe=500)], 'L')[0], 'unplaced')
        self.assertEqual(rs.placement([], 'R')[0], 'unplaced')


class GeometryTests(unittest.TestCase):
    def test_forward_and_reverse_give_the_sequence_between_the_flanks(self):
        self.assertEqual(rs.geometry('c', 1100, '+', 'c', 2100, '+', 1000), ('proper', 1100, 2100))
        self.assertEqual(rs.geometry('c', 2100, '-', 'c', 1100, '-', 1000), ('proper', 1100, 2100))
        self.assertEqual(rs.geometry('c', 1100, '+', 'c', 1100, '+', 1000), ('proper', 1100, 1100))  # empty

    def test_junctions(self):
        self.assertEqual(rs.geometry('c', 1100, '+', 'd', 2100, '+', 1000)[0], 'junction')
        self.assertEqual(rs.geometry('c', 1100, '+', 'c', 2100, '-', 1000)[0], 'junction')
        self.assertEqual(rs.geometry('c', 1100, '+', 'c', 50000, '+', 1000)[0], 'junction')
        self.assertEqual(rs.geometry('c', 2100, '+', 'c', 1100, '+', 1000)[0], 'junction')


class CanonicalTests(unittest.TestCase):
    def test_both_strands_of_one_position_give_the_same_key(self):
        fwd = dict(left_contig='c', left_pos=1100, left_strand='+', right_contig='c', right_pos=2100, right_strand='+')
        rev = dict(left_contig='c', left_pos=2100, left_strand='-', right_contig='c', right_pos=1100, right_strand='-')
        self.assertEqual(rs.canonical(fwd), rs.canonical(rev))
        self.assertEqual(rs.geometry(*rs.canonical(rev), 1000), ('proper', 1100, 2100))


class ClusterTests(unittest.TestCase):
    def test_observations_within_tolerance_merge(self):
        obs = [dict(c_left_pos=1000, c_right_pos=2000), dict(c_left_pos=1030, c_right_pos=2010),
               dict(c_left_pos=1500, c_right_pos=2500)]
        self.assertEqual(sorted(len(c) for c in rs.cluster(obs)), [1, 2])


class ShiftedTests(unittest.TestCase):
    def test_position_overlapping_a_known_one_is_marked_and_keeps_its_class(self):
        rows = [dict(cls='known_catalog_position', contig='c', start=1000, end=2000),
                dict(cls='not_in_assembly:other_sequence', contig='c', start=1150, end=2000),
                dict(cls='catalog_omission', contig='c', start=5000, end=6000),
                dict(cls='not_in_assembly:junction', contig='', start='', end='', left_contig='c', left_pos=1990,
                     right_contig='d', right_pos=10)]
        out = rs.mark_shifted(rows)
        self.assertEqual([r['cls'] for r in out], ['known_catalog_position', 'shifted_known_position',
                                                    'catalog_omission', 'shifted_known_position'])
        self.assertEqual(out[1]['raw_cls'], 'not_in_assembly:other_sequence')


if __name__ == '__main__':
    unittest.main()
