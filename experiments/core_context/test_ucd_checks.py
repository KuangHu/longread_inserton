"""Tests of Check 1 (source extraction) and Check 2 (independent U-X-D comparison) in ucd_checks.py (needs blastn)."""
import hashlib
import random
import tempfile
import unittest
from pathlib import Path

import ucd_checks as ck
import ucd_compare as uc

R = random.Random(20261107)


def seq(n):
    return ''.join(R.choice('ACGT') for _ in range(n))


def mutate(s, rate):
    r = random.Random(len(s) * 11 + int(rate * 1000))
    return ''.join(r.choice('ACGT'.replace(c, '')) if r.random() < rate else c for c in s)


X = seq(1200)
U1, D1, U2, D2 = seq(1000), seq(1000), seq(1000), seq(1000)
CFG, CC = uc.CONFIG, ck.CHECK_CONFIG


def occ(i, u, x, d, sample='', contig=''):
    n = len(u)
    return dict(id=i, U=u, X=x, D=d, sample=sample, contig=contig,
                ivs=dict(U=(0, n), X=(n, n + len(x)), D=(n + len(x), n + len(x) + len(d))))


def check(a, b):
    with tempfile.TemporaryDirectory() as t:
        return ck.check_ucd(a, b, CFG, CC, t)[0]


def fasta(tmp, recs, width=60):
    p = Path(tmp) / 'g.fna'
    p.write_text(''.join('>%s desc\n%s\n' % (n, '\n'.join(s[i:i + width] for i in range(0, len(s), width)))
                         for n, s in recs.items()))
    return p, hashlib.sha256(p.read_bytes()).hexdigest()


class Check1(unittest.TestCase):
    def setUp(self):
        self.t = tempfile.TemporaryDirectory()
        self.g = U1 + X + D1
        self.p, self.h = fasta(self.t.name, {'c1': seq(500), 'c2': self.g})

    def tearDown(self):
        self.t.cleanup()

    def test_both_strands_match_the_screen(self):
        src = ck.Source(self.p, self.h)
        for st in '+-':
            screen = uc.sides_of(uc.read_fasta(self.p)['c2'], 1000, 2200, st, 1000)
            c = ck.check1_side(src, 'c2', 1000, 2200, st, 1000, screen)
            self.assertEqual(c['status'], 'pass', c['reasons'])
        self.assertEqual(c['seqs']['X'], uc.rc(X))
        self.assertEqual(c['seqs']['U'], uc.rc(D1))
        self.assertEqual(c['ivs']['U'], (2200, 3200))

    def test_fetch_across_line_breaks(self):
        src = ck.Source(self.p, self.h)
        for lo, hi in ((0, 1), (59, 61), (119, 240), (3199, 3200)):
            self.assertEqual(src.fetch('c2', lo, hi), self.g[lo:hi])

    def test_errors_are_reported_not_repaired(self):
        src = ck.Source(self.p, self.h)
        self.assertEqual(ck.check1_side(src, 'c2', 3000, 3300, '+', 1000)['reasons'], ['interval_out_of_bounds'])
        self.assertEqual(ck.check1_side(src, 'nope', 0, 10, '+', 1000)['reasons'], ['contig_missing'])
        self.assertEqual(ck.check1_side(ck.Source(self.p, '0' * 64), 'c2', 0, 10, '+', 1000)['reasons'],
                         ['fasta_sha256_mismatch'])
        screen = uc.sides_of(self.g, 1000, 2200, '+', 1000)
        screen['X'] = 'A' + screen['X'][1:] if screen['X'][0] != 'A' else 'C' + screen['X'][1:]
        r = ck.check1_side(src, 'c2', 1000, 2200, '+', 1000, screen)
        self.assertEqual((r['status'], r['reasons']), ('fail', ['screen_sequence_differs:X@0']))

    def test_contig_edge_is_recorded_not_an_error(self):
        c = ck.check1_side(ck.Source(self.p, self.h), 'c2', 300, 1000, '+', 1000)
        self.assertEqual((c['status'], c['edge'], c['ivs']['U']), ('pass', 'U_at_contig_edge', (0, 300)))

    def test_irregular_lines_refused(self):
        p = Path(self.t.name) / 'bad.fna'
        p.write_text('>c\nACGTACGT\nACG\nACGTACGT\n')
        src = ck.Source(p, hashlib.sha256(p.read_bytes()).hexdigest())
        self.assertEqual(ck.check1_side(src, 'c', 0, 5, '+', 2)['reasons'], ['irregular_line_width'])

    def test_iupac_reverse_complement(self):
        self.assertEqual(ck.revcomp('ACGTRYKMN'), 'NKMRYACGT')


class Check2(unittest.TestCase):
    def test_ucd_supported_and_same_context(self):
        r = check(occ('a', U1, X, D1), occ('b', U2, mutate(X, 0.02), D2))
        self.assertEqual(r['ucd_check'], 'ucd_supported')
        self.assertGreaterEqual(min(r['chk_insertion_cov_a'], r['chk_insertion_cov_b']), 0.9)
        self.assertEqual(check(occ('a', U1, X, D1), occ('b', mutate(U1, 0.05), X, mutate(D1, 0.05)))['ucd_check'],
                         'same_context')
        self.assertEqual(check(occ('a', U1, X, D1), occ('b', mutate(U1, 0.03), X, D2))['ucd_check'],
                         'only_one_side_different')

    def test_reverse_complement_unified(self):
        b = occ('b', uc.rc(D2), uc.rc(X), uc.rc(U2))
        r = check(occ('a', U1, X, D1), b)
        self.assertEqual((r['ucd_check'], r['chk_b_flipped']), ('ucd_supported', True))
        self.assertEqual(check(occ('a', U1, X, D1), occ('b', uc.rc(D1), uc.rc(X), uc.rc(U1)))['ucd_check'],
                         'same_context')
        self.assertEqual(check(occ('a', U1, X, D1), occ('b', U1, uc.rc(X), D1))['ucd_check'],
                         'swapped_or_rearranged_context')

    def test_fragments_out_of_order_are_not_a_full_match(self):
        # X' = second half + first half: every base matches, but not collinearly
        r = check(occ('a', U1, X, D1), occ('b', U2, X[600:] + X[:600], D2))
        self.assertEqual(r['ucd_check'], 'rearranged_or_fragmented_insertion')
        self.assertLess(r['chk_insertion_cov_a'], 0.9)
        self.assertGreaterEqual(r['chk_insertion_union_cov_a'], 0.9)
        # an inverted middle segment
        r = check(occ('a', U1, X, D1), occ('b', U2, X[:400] + uc.rc(X[400:800]) + X[800:], D2))
        self.assertEqual(r['ucd_check'], 'rearranged_or_fragmented_insertion')

    def test_repeated_hit_not_counted_twice(self):
        # X' = the first half of X twice: union over X' is full, X is only half covered
        r = check(occ('a', U1, X, D1), occ('b', U2, X[:600] + X[:600], D2))
        self.assertEqual(r['ucd_check'], 'insertion_partial')
        self.assertLessEqual(r['chk_insertion_cov_a'], 0.51)

    def test_partial_dissimilar_insufficient(self):
        self.assertEqual(check(occ('a', U1, X, D1), occ('b', U2, X[:600], D2))['ucd_check'], 'insertion_partial')
        self.assertEqual(check(occ('a', U1, X, D1), occ('b', U2, seq(1200), D2))['ucd_check'], 'insertion_dissimilar')
        self.assertEqual(check(occ('a', U1, X, D1), occ('b', 'N' * 700 + U2[:300], X, D2))['ucd_check'],
                         'insufficient_context')

    def test_agrees_with_screen_on_simple_cases(self):
        cases = [(occ('a', U1, X, D1), occ('b', U2, mutate(X, 0.02), D2)),
                 (occ('a', U1, X, D1), occ('b', mutate(U1, 0.05), X, D2)),
                 (occ('a', U1, X, D1), occ('b', U1, uc.rc(X), D1)),
                 (occ('a', U1, X, D1), occ('b', U2, X[:600], D2))]
        for a, b in cases:
            self.assertEqual(check(a, b)['ucd_check'], uc.compare(a, b)['verdict'])

    def test_swap_invariance(self):
        a, b = occ('a', U1, X, D1), occ('b', U2, mutate(X, 0.02), D2)
        self.assertEqual(check(a, b), check(b, a))

    def test_chain_overlap_counted_once(self):
        hs = [dict(qs=0, qe=600, ss=0, se=600, nident=600, length=600),
              dict(qs=590, qe=1200, ss=590, se=1200, nident=610, length=610)]
        ch = ck.best_chain(hs, 'plus', 20)
        self.assertEqual(len(ch), 2)
        self.assertEqual(ck.cover(ch, 1200, 1200)['cov_a'], 1.0)
        self.assertEqual(len(ck.best_chain(hs, 'plus', 5)), 1)


if __name__ == '__main__':
    unittest.main()
