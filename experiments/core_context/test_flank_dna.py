import random
import unittest

import core_review as cr
import flank_dna as fd
import s6_review as rv

rnd = random.Random(7)


def seq(n):
    return ''.join(rnd.choice('ACGT') for _ in range(n))


class OuterWindowTests(unittest.TestCase):
    def test_windows_are_cut_before_alignment(self):
        f = dict(up='U' * 2500 + 'i' * 500, down='j' * 500 + 'D' * 2500)
        o = fd.outer_windows(f)
        self.assertEqual((o['up'], o['down']), ('U' * 2500, 'D' * 2500))
        self.assertEqual(fd.outer_windows(dict(up='A' * 300, down='C' * 300)), dict(up='', down=''))


class CounterexampleTests(unittest.TestCase):
    def test_outer_matching_only_the_other_inner_window_is_not_shared(self):
        # review counterexample: each copy's outer window holds the other copy's inner (excluded) window, nothing else
        # is shared. v0.1 aligned the full flanks and called this shared; v2 cuts both outer windows first.
        ix, iy = seq(500), seq(500)
        x = dict(up=iy + seq(500) + ix, down=seq(1500))
        y = dict(up=ix + seq(500) + iy, down=seq(1500))
        self.assertEqual(rv.flank_check(x, y)['up'][2], 'shared')        # the v0.1 behaviour, kept for the record
        v = fd.outer_flank_check(x, y)
        self.assertEqual(v['dna_up'], 'not_shared')
        self.assertEqual(v['cov_xup_yup_plus'], 0.0)


class OrientationTests(unittest.TestCase):
    def test_same_side_forward(self):
        s = seq(2500)
        v = fd.outer_flank_check(dict(up=s + seq(500), down=seq(3000)), dict(up=s + seq(500), down=seq(3000)))
        self.assertEqual((v['dna_up'], v['dna_down']), ('shared', 'not_shared'))
        self.assertFalse(v['other_orientation_match'])

    def test_same_side_reverse_is_recorded_not_counted_as_shared(self):
        s = seq(2500)
        v = fd.outer_flank_check(dict(up=s + seq(500), down=seq(3000)),
                                 dict(up=cr.revcomp(s) + seq(500), down=seq(3000)))
        self.assertEqual(v['dna_up'], 'not_shared')
        self.assertGreaterEqual(v['cov_xup_yup_minus'], 0.99)
        self.assertTrue(v['other_orientation_match'])

    def test_cross_side_match(self):
        # x's upstream outer window is y's downstream outer window (e.g. a copy in the opposite orientation in a
        # shared segment, or two copies flanking one shared region)
        s = seq(2500)
        v = fd.outer_flank_check(dict(up=s + seq(500), down=seq(3000)), dict(up=seq(3000), down=seq(500) + s))
        self.assertEqual((v['dna_up'], v['dna_down']), ('not_shared', 'not_shared'))
        self.assertGreaterEqual(v['cov_xup_ydown_plus'], 0.99)
        self.assertGreaterEqual(v['cross_side_max'], 0.99)

    def test_short_outer_flank(self):
        v = fd.outer_flank_check(dict(up=seq(1200), down=seq(3000)), dict(up=seq(3000), down=seq(3000)))
        self.assertEqual(v['dna_up'], 'short_flank')
        self.assertEqual(v['outer_len_up'], 700)


if __name__ == '__main__':
    unittest.main()
