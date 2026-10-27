import unittest

import protein_context as pc


def contig(strands, circular=False):
    return [dict(gene_uid='g%d' % i, strand=s, circular_contig=str(circular)) for i, s in enumerate(strands)]


class ContextTests(unittest.TestCase):
    def test_minus_strand_core_swaps_sides_and_strand_sense(self):
        gs = contig('++-++++')
        labels = {g['gene_uid']: g['gene_uid'] for g in gs}
        up, down = pc.context(gs, 3, 3, False, labels)
        self.assertEqual([x[2] for x in up], ['g4', 'g5', 'g6'])      # upstream of a minus-strand core is to the right
        self.assertEqual([x[2] for x in down], ['g2', 'g1', 'g0'])
        self.assertEqual([x[1] for x in down], [True, False, False])   # g2 is on the minus strand like the core

    def test_linear_edge_is_short_and_circular_wraps(self):
        gs = contig('+++++')
        labels = {g['gene_uid']: g['gene_uid'] for g in gs}
        self.assertEqual(len(pc.context(gs, 0, 0, True, labels)[0]), 0)
        gs = contig('+++++', circular=True)
        up, down = pc.context(gs, 0, 0, True, labels)
        self.assertEqual([x[2] for x in up], ['g4', 'g3', 'g2'])
        self.assertEqual([x[2] for x in down], ['g1', 'g2', 'g3'])

    def test_wrap_stops_at_core_on_small_circle(self):
        gs = contig('+++', circular=True)
        labels = {g['gene_uid']: g['gene_uid'] for g in gs}
        up, down = pc.context(gs, 0, 0, True, labels)
        self.assertEqual(len(up), 2)
        self.assertEqual(len(down), 2)


class RelationTests(unittest.TestCase):
    def occ(self, up, down):
        return dict(up=[(x, True, x) for x in up], down=[(x, True, x) for x in down])

    def test_relations(self):
        a = self.occ(['A', 'B', 'C'], ['D', 'E', 'F'])
        self.assertEqual(pc.relation(a, self.occ(['X', 'B', 'Y'], ['Z', 'Z', 'F'])), 'similar_context')
        self.assertEqual(pc.relation(a, self.occ(['X', 'Y', 'Z'], ['P', 'Q', 'R'])), 'different_context_candidate')
        self.assertEqual(pc.relation(a, self.occ(['A', 'Y', 'Z'], ['P', 'Q', 'R'])), 'one_side_shared_unresolved')
        self.assertEqual(pc.relation(a, self.occ(['A', 'B'], ['D', 'E', 'F'])), 'edge_or_missing_context')

    def test_unlabelled_neighbours_never_shared(self):
        a = self.occ([None, None, None], ['D', 'E', 'F'])
        self.assertEqual(pc.relation(a, self.occ([None, None, None], ['D', 'E', 'F'])), 'one_side_shared_unresolved')


class PairKeyTests(unittest.TestCase):
    def test_same_structure_on_either_strand_shares_key(self):
        k1, f1 = pc.pair_key('PF_a', '+', 'PF_b', '-')
        k2, f2 = pc.pair_key('PF_b', '+', 'PF_a', '-')      # same pair read on the other strand
        self.assertEqual(k1, k2)
        self.assertNotEqual(f1, f2)
        self.assertNotEqual(pc.pair_key('PF_a', '+', 'PF_b', '+')[0], k1)   # relative orientation is part of the key


if __name__ == '__main__':
    unittest.main()
