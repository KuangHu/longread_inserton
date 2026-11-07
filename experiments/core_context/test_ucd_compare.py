"""Tests of the U-insertion-D comparison (ucd_compare.py) on seeded synthetic sequences (needs blastn on PATH)."""
import random
import unittest

import ucd_compare as uc

R = random.Random(20261106)


def seq(n):
    return ''.join(R.choice('ACGT') for _ in range(n))


def mutate(s, rate):
    r = random.Random(len(s) * 7 + int(rate * 1000))
    return ''.join(r.choice('ACGT'.replace(c, '')) if r.random() < rate else c for c in s)


X = seq(1200)
U1, D1, U2, D2 = seq(1000), seq(1000), seq(1000), seq(1000)


def occ(i, u, x, d):
    return dict(id=i, U=u, X=x, D=d)


def flip(o):
    """The same occurrence as seen on the other strand."""
    return dict(o, U=uc.rc(o['D']), X=uc.rc(o['X']), D=uc.rc(o['U']))


def verdict(a, b):
    return uc.compare(a, b)['verdict']


class Structure(unittest.TestCase):
    def test_same_insertion_both_sides_different(self):
        r = uc.compare(occ('a', U1, X, D1), occ('b', U2, mutate(X, 0.02), D2))
        self.assertEqual(r['verdict'], 'ucd_supported')
        self.assertGreaterEqual(r['insertion_identity'], 97)
        self.assertEqual((r['upstream'], r['downstream']), ('different', 'different'))

    def test_same_insertion_same_context(self):
        r = uc.compare(occ('a', U1, X, D1), occ('b', mutate(U1, 0.05), X, mutate(D1, 0.05)))
        self.assertEqual(r['verdict'], 'same_context')

    def test_only_one_side_different(self):
        self.assertEqual(verdict(occ('a', U1, X, D1), occ('b', mutate(U1, 0.03), X, D2)), 'only_one_side_different')
        self.assertEqual(verdict(occ('a', U1, X, D1), occ('b', U2, X, mutate(D1, 0.03))), 'only_one_side_different')

    def test_reverse_complement_side_is_unified(self):
        r = uc.compare(occ('a', U1, X, D1), flip(occ('b', U2, X, D2)))
        self.assertEqual((r['verdict'], r['b_flipped']), ('ucd_supported', True))
        # the same context seen from the other strand is the same context, not a different one
        self.assertEqual(verdict(occ('a', U1, X, D1), flip(occ('b', U1, X, D1))), 'same_context')

    def test_inverted_insertion_in_one_context_is_rearranged(self):
        # B = U1 - rc(X) - D1: after orienting by the insertion, B's flanks are rc(D1) / rc(U1): cross-similar
        r = uc.compare(occ('a', U1, X, D1), occ('b', U1, uc.rc(X), D1))
        self.assertEqual(r['verdict'], 'swapped_or_rearranged_context')
        # flanks swapped between the two sides
        self.assertEqual(verdict(occ('a', U1, X, D1), occ('b', D1, X, U1)), 'swapped_or_rearranged_context')

    def test_partial_and_dissimilar_insertion(self):
        r = uc.compare(occ('a', U1, X, D1), occ('b', U2, X[:600], D2))
        self.assertEqual(r['verdict'], 'insertion_partial')
        self.assertLess(r['insertion_cov_a'], 0.9)
        self.assertEqual(verdict(occ('a', U1, X, D1), occ('b', U2, seq(1200), D2)), 'insertion_dissimilar')
        self.assertEqual(verdict(occ('a', U1, X, D1), occ('b', U2, mutate(X, 0.2), D2)), 'insertion_dissimilar')

    def test_missing_or_low_complexity_context_is_not_different(self):
        self.assertEqual(verdict(occ('a', U1, X, D1), occ('b', U2[:300], X, D2)), 'insufficient_context')
        self.assertEqual(verdict(occ('a', U1, X, D1), occ('b', 'N' * 700 + U2[:300], X, D2)), 'insufficient_context')
        self.assertEqual(verdict(occ('a', U1, X, D1), occ('b', 'AT' * 500, X, D2)), 'insufficient_context')

    def test_reused_physical_flank_in_one_assembly(self):
        # tandem copies 400 bp apart on one contig: B's upstream flank is A's insertion and downstream flank
        g = U1 + X + D1[:400] + mutate(X, 0.01) + D2
        lo_b = len(U1) + len(X) + 400
        A = dict(uc.sides_of(g, len(U1), len(U1) + len(X), '+', 1000), id='a', sample='s', contig='c')
        B = dict(uc.sides_of(g, lo_b, lo_b + len(X), '+', 1000), id='b', sample='s', contig='c')
        r = uc.compare(A, B)
        self.assertEqual(r['verdict'], 'shared_physical_context')
        self.assertIn('D_A~U_B', r['reasons'])
        # the same sequences without coordinates are judged on sequence alone
        self.assertNotEqual(verdict(dict(A, sample=''), dict(B, sample='')), 'shared_physical_context')


class Invariance(unittest.TestCase):
    CASES = [(occ('a', U1, X, D1), occ('b', U2, mutate(X, 0.02), D2)),
             (occ('a', U1, X, D1), occ('b', mutate(U1, 0.05), X, D2)),
             (occ('a', U1, X, D1), occ('b', U1, uc.rc(X), D1)),
             (occ('a', U1, X, D1), occ('b', U2, X[:600], D2))]

    def test_swap_sides(self):
        for a, b in self.CASES:
            self.assertEqual(uc.compare(a, b), uc.compare(b, a))

    def test_reverse_complement_both_or_either(self):
        for a, b in self.CASES:
            v = verdict(a, b)
            self.assertEqual(verdict(flip(a), b), v)
            self.assertEqual(verdict(a, flip(b)), v)
            self.assertEqual(verdict(flip(a), flip(b)), v)

    def test_input_order(self):
        rows = [verdict(a, b) for a, b in self.CASES]
        order = list(range(len(self.CASES)))
        random.Random(1).shuffle(order)
        self.assertEqual([verdict(*self.CASES[i]) for i in order], [rows[i] for i in order])


class Helpers(unittest.TestCase):
    def test_sides_of_minus_strand(self):
        g = U1 + X + D1
        d = uc.sides_of(g, 1000, 2200, '-', 1000)
        self.assertEqual((d['X'], d['U'], d['D']), (uc.rc(X), uc.rc(D1), uc.rc(U1)))
        self.assertEqual((d['U_iv'], d['D_iv']), ((2200, 3200), (0, 1000)))

    def test_contig_edge_truncates(self):
        d = uc.sides_of(U1[:200] + X + D1, 200, 1400, '+', 1000)
        self.assertEqual(len(d['U']), 200)

    def test_rc_iupac_and_case(self):
        self.assertEqual(uc.rc('ACGTRYKMSWBDHVN'), 'NBDHVWSKMRYACGT')
        self.assertEqual(uc.rc('acgtrykm'), 'kmryacgt')
        for s in ('ACGTRYSWKMBDHVN' * 3, 'acgtnrykmswbdhv', 'AcGtRy'):
            self.assertEqual(uc.rc(uc.rc(s)), s)
        self.assertEqual(uc.rc(X), X.translate(str.maketrans('ACGT', 'TGCA'))[::-1])  # A/C/G/T unchanged

    def test_minus_strand_extraction_with_ambiguity_codes(self):
        g = U1[:500] + 'ACRYGT' + X[:300] + 'KMNSW' + D1[:500]
        d = uc.sides_of(g, 500, 811, '-', 500)
        comp = dict(zip('ACGTRYKMSWN', 'TGCAYRMKSWN'))
        self.assertEqual(d['X'], ''.join(comp[c] for c in reversed(g[500:811])))
        self.assertNotIn('?', d['X'])

    def test_union_coverage_counts_overlaps_once(self):
        self.assertEqual(uc.union_len([(0, 100), (50, 150), (300, 310)]), 160)


if __name__ == '__main__':
    unittest.main()
