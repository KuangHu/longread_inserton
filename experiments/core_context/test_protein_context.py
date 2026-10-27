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

    def test_unknown_labels_are_not_different(self):
        # one side entirely unlabelled: unknown, not evidence of a different context
        a = self.occ([None, None, None], ['D', 'E', 'F'])
        self.assertEqual(pc.relation(a, self.occ(['X', 'Y', 'Z'], ['P', 'Q', 'R'])), 'insufficient_context_labels')
        self.assertEqual(pc.relation(a, self.occ([None, None, None], ['D', 'E', 'F'])), 'insufficient_context_labels')
        # a single unlabelled neighbour on an unshared side also blocks the different call
        b = self.occ(['A', None, 'C'], ['D', 'E', 'F'])
        self.assertEqual(pc.relation(b, self.occ(['X', 'Y', 'Z'], ['P', 'Q', 'R'])), 'insufficient_context_labels')
        # a shared label still counts on a partly unlabelled side
        self.assertEqual(pc.relation(b, self.occ(['A', 'Y', 'Z'], ['D', 'Q', 'R'])), 'similar_context')

    def test_one_side_shared_needs_labelled_other_side(self):
        a = self.occ(['A', 'B', 'C'], ['D', 'E', None])
        self.assertEqual(pc.relation(a, self.occ(['A', 'Y', 'Z'], ['P', 'Q', 'R'])), 'insufficient_context_labels')
        a = self.occ(['A', 'B', 'C'], ['D', 'E', 'F'])
        self.assertEqual(pc.relation(a, self.occ(['A', 'Y', 'Z'], ['P', 'Q', 'R'])), 'one_side_shared_unresolved')


class CompleteCdsTests(unittest.TestCase):
    def unit(self, **kw):
        u = dict(type='CDS', partial=False, pseudo=False, crosses_origin=False)
        u.update(kw)
        return u

    def test_s4_strict_rule(self):
        self.assertTrue(pc.is_complete(self.unit(), 'match', [], True))
        self.assertTrue(pc.is_complete(self.unit(), 'match_with_transl_except', [], True))
        self.assertFalse(pc.is_complete(self.unit(), 'match', ['no_start'], True))      # translation matches, bad start
        self.assertFalse(pc.is_complete(self.unit(), 'match', ['internal_stop'], True))
        self.assertFalse(pc.is_complete(self.unit(partial=True), 'match', [], True))
        self.assertFalse(pc.is_complete(self.unit(crosses_origin=True), 'match', [], True))
        self.assertFalse(pc.is_complete(self.unit(), 'mismatch', ['translation_mismatch'], True))
        self.assertFalse(pc.is_complete(self.unit(), 'match', [], False))


class SpanTests(unittest.TestCase):
    def test_two_gene_span_includes_spacing(self):
        g, h = dict(start0='100', end0='1300'), dict(start0='1500', end0='2550')
        self.assertEqual(pc.pair_span(g, h), 2450)                      # 1200 + 200 spacing + 1050
        self.assertTrue(pc.pair_span(g, h) < pc.TARGET_MAX)
        self.assertFalse(pc.pair_span(g, dict(start0='1700', end0='2650')) < pc.TARGET_MAX)   # spacing pushes it over


class ScopeTests(unittest.TestCase):
    def occ(self, oid, span, up, down):
        return dict(occurrence_id=oid, sample='s' + oid, span=span, in_target_length=span < pc.TARGET_MAX, spacing='',
                    products='p', up=[(x, True, x) for x in up], down=[(x, True, x) for x in down])

    def test_target_scope_drops_long_instances(self):
        xs = [self.occ('1', 900, 'ABC', 'DEF'), self.occ('2', 3000, 'XYZ', 'PQR'), self.occ('3', 950, 'ABC', 'DEF')]
        groups = {('single_gene', 'PF_x'): xs}
        pairs, cores = pc.compare(groups, 'all', lambda o: True)
        self.assertEqual(cores[0]['different_context_candidate'], 2)
        pairs, cores = pc.compare(groups, 'target_lt2500', lambda o: o['in_target_length'])
        self.assertEqual(cores[0]['occurrences'], 2)
        self.assertEqual(cores[0]['different_context_candidate'], 0)
        self.assertEqual(cores[0]['span_range'], '900..950')
        pairs, cores = pc.compare({('single_gene', 'PF_y'): xs[1:2] + [self.occ('4', 2600, 'ABC', 'DEF')]},
                                  'target_lt2500', lambda o: o['in_target_length'])
        self.assertEqual(cores, [])


class PairKeyTests(unittest.TestCase):
    def test_same_structure_on_either_strand_shares_key(self):
        k1, f1 = pc.pair_key('PF_a', '+', 'PF_b', '-')
        k2, f2 = pc.pair_key('PF_b', '+', 'PF_a', '-')      # same pair read on the other strand
        self.assertEqual(k1, k2)
        self.assertNotEqual(f1, f2)
        self.assertNotEqual(pc.pair_key('PF_a', '+', 'PF_b', '+')[0], k1)   # relative orientation is part of the key


if __name__ == '__main__':
    unittest.main()
