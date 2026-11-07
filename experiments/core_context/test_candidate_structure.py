"""Tests of the full-insertion check rules (candidate_structure.py): every ORF kept, gene count as a range, target."""
import unittest

import candidate_structure as cs

GENE = 'ATG' + 'GCT' * 80 + 'TAA'  # one complete 81-aa ORF on the plus strand


class OrfTests(unittest.TestCase):
    def test_overlapping_orfs_are_kept_not_dropped(self):
        o = cs.orfs(GENE + 'A' * 30)
        self.assertEqual(sum(x['primary'] for x in o), 1)
        self.assertTrue(len(o) > 1)  # the overlapping frames are listed with primary False
        self.assertTrue(all(not x['primary'] for x in o if not x['complete']))

    def test_partial_orf_at_the_end_is_listed(self):
        o = cs.orfs('C' + 'GCT' * 70)  # no start, no stop: partial
        self.assertTrue(any(not x['complete'] for x in o))


class TargetTests(unittest.TestCase):
    def test_gene_range(self):
        self.assertEqual(cs.target(0, 1, 0, 1300), 'met')
        self.assertEqual(cs.target(0, 1, 1, 1300), 'met')  # 1-2 genes either way
        self.assertEqual(cs.target(0, 1, 2, 1300), 'unresolved_unexplained_orfs')  # could be 3 genes
        self.assertEqual(cs.target(0, 0, 1, 1300), 'unresolved_unexplained_orfs')  # no known complete gene
        self.assertEqual(cs.target(0, 3, 0, 2000), 'not_met')
        self.assertEqual(cs.target(0, 0, 0, 900), 'not_met')

    def test_size_strict_and_ambiguity(self):
        self.assertEqual(cs.target(0, 1, 0, 2499), 'met')
        self.assertEqual(cs.target(0, 1, 0, 2500), 'not_met')
        self.assertEqual(cs.target(4, 1, 0, 1300), 'unresolved_ambiguous_consensus')


if __name__ == '__main__':
    unittest.main()
