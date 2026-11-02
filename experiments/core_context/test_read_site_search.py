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


class ClusterTests(unittest.TestCase):
    def test_observations_within_tolerance_merge(self):
        obs = [dict(left_pos=1000, right_pos=2000), dict(left_pos=1030, right_pos=2010),
               dict(left_pos=1500, right_pos=2500)]
        self.assertEqual(sorted(len(c) for c in rs.cluster(obs)), [1, 2])


if __name__ == '__main__':
    unittest.main()
