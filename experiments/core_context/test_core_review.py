import unittest

import core_review as cr


def contig(strands, circular=False):
    return [dict(gene_uid='g%d' % i, strand=s, circular_contig=str(circular)) for i, s in enumerate(strands)]


class WindowTests(unittest.TestCase):
    def test_circular_slice_wraps(self):
        self.assertEqual(cr.cslice('ABCDEFGH', -2, 3), 'GHABC')
        self.assertEqual(cr.cslice('ABCDEFGH', 6, 10), 'GHAB')

    def test_window_orientation_and_core_offsets(self):
        seq = 'AAAACCGGTTTT'
        sub, c0, c1 = cr.window(seq, 4, 8, True, 2, False)
        self.assertEqual((sub, sub[c0:c1]), ('AACCGGTT', 'CCGG'))
        sub, c0, c1 = cr.window(seq, 4, 8, False, 2, False)
        self.assertEqual(sub[c0:c1], cr.revcomp('CCGG'))
        self.assertEqual(sub, cr.revcomp('AACCGGTT'))
        sub, c0, c1 = cr.window(seq, 1, 3, True, 3, False)      # linear: clipped at the contig start
        self.assertEqual((sub, c0, c1), ('AAAACC', 1, 3))
        sub, c0, c1 = cr.window('AAAACCGG', 1, 3, True, 3, True)   # circular, 8 bp: not longer than core + 2w, clipped
        self.assertEqual((sub, c0), ('AAAACC', 1))
        sub, c0, c1 = cr.window('ACGT' * 10, 1, 3, True, 3, True)
        self.assertEqual((len(sub), c0, c1), (8, 3, 5))


class ChainTests(unittest.TestCase):
    def test_extension_beyond_core_and_chaining(self):
        # A: core at [100, 400) of 700 bp; B: core at [300, 600) of 900 bp; shared from A 60..460 in two HSPs
        hsps = [(61, 250, 261, 450), (261, 460, 461, 660)]
        pieces = cr.chain(hsps, 100, 400, 300, 600)
        self.assertEqual(pieces, [(60, 250, 260, 450), (260, 460, 460, 660)])
        e = cr.extensions(pieces, 100, 400, 700, 300, 600, 900)
        self.assertEqual((e['a_up'], e['a_down'], e['b_up'], e['b_down']), (40, 60, 40, 60))
        self.assertEqual((e['up_cover'], e['down_cover']), ('covered', 'covered'))
        self.assertFalse(e['up_reaches_end'] or e['down_reaches_end'])

    def test_off_diagonal_far_hsp_not_chained_and_end_flag(self):
        hsps = [(101, 400, 301, 600), (600, 690, 50, 140)]      # second HSP: elsewhere in B (e.g. a nested element)
        pieces = cr.chain(hsps, 100, 400, 300, 600)
        self.assertEqual(pieces, [(100, 400, 300, 600)])
        e = cr.extensions(pieces, 100, 400, 700, 300, 600, 900)
        self.assertEqual((e['a_up'], e['a_down']), (0, 0))      # outermost piece stops exactly at both ends: measured 0
        self.assertEqual(e['up_cover'], 'covered_at_chain_end')
        e = cr.extensions([(2, 400, 202, 600)], 100, 400, 700, 300, 600, 900)
        self.assertTrue(e['up_reaches_end'])
        self.assertIsNone(cr.chain([(500, 690, 50, 240)], 100, 400, 300, 600))

    def test_local_hit_inside_core_is_unknown_not_zero(self):
        # review counterexample 1: 1,000 bp core, a 100 bp local alignment in its middle
        pieces = cr.chain([(501, 600, 501, 600)], 100, 1100, 100, 1100)
        e = cr.extensions(pieces, 100, 1100, 1300, 100, 1100, 1300)
        self.assertEqual((e['a_up'], e['a_down'], e['b_up'], e['b_down']), (None, None, None, None))
        self.assertEqual((e['up_cover'], e['down_cover']), ('not_reached', 'not_reached'))
        self.assertFalse(e['up_reaches_end'] or e['down_reaches_end'])
        self.assertEqual(cr.side_boundary([e['a_up']], False), (None, 'end_not_covered', 1))
        # covered in A but not in B is still unknown
        e = cr.extensions([(50, 1150, 120, 1150)], 100, 1100, 1300, 100, 1100, 1300)
        self.assertIsNone(e['a_up'])
        self.assertEqual(e['a_down'], 50)

    def test_core_end_in_chain_gap_is_unknown(self):
        # review counterexample 2: two chained pieces, the core start (200) lies in the 20 bp unaligned gap between them
        pieces = cr.chain([(101, 190, 101, 190), (211, 600, 211, 600)], 200, 500, 200, 500)
        self.assertEqual(len(pieces), 2)
        e = cr.extensions(pieces, 200, 500, 800, 200, 500, 800)
        self.assertEqual(e['up_cover'], 'end_in_chain_gap')
        self.assertIsNone(e['a_up'])
        self.assertEqual((e['down_cover'], e['a_down']), ('covered', 100))
        # the core start exactly where an inner piece begins also borders the gap: unknown, not 100 bp
        pieces = cr.chain([(81, 180, 81, 180), (201, 600, 201, 600)], 200, 500, 200, 500)
        e = cr.extensions(pieces, 200, 500, 800, 200, 500, 800)
        self.assertEqual((e['up_cover'], e['a_up']), ('end_in_chain_gap', None))
        # an end inside one piece with a gap further out is measured, and the gap is counted
        pieces = cr.chain([(61, 120, 61, 120), (141, 600, 141, 600)], 200, 500, 200, 500)
        e = cr.extensions(pieces, 200, 500, 800, 200, 500, 800)
        self.assertEqual((e['up_cover'], e['a_up'], e['up_gaps_beyond_end']), ('covered', 140, 1))

    def test_side_boundary(self):
        self.assertEqual(cr.side_boundary([40, 45, 50], False), (45, 'resolved', 0))
        self.assertEqual(cr.side_boundary([40, 400], False)[1], 'inconsistent_between_partners')
        self.assertEqual(cr.side_boundary([40], True)[1], 'reaches_window_or_contig_end')
        self.assertEqual(cr.side_boundary([], False), (None, 'end_not_covered', 0))
        # agreeing known partners do not resolve a side while another partner does not cover the end
        self.assertEqual(cr.side_boundary([40, 41, None], False), (40, 'end_not_covered_by_some_partners', 1))


class CopyPatternTests(unittest.TestCase):
    def test_patterns_are_descriptive(self):
        self.assertEqual(cr.copy_pattern([2, 2, 2]), 'same_count_in_all_assemblies')
        self.assertEqual(cr.copy_pattern([2, 0, 2]), 'same_count_where_present')
        self.assertEqual(cr.copy_pattern([4, 2, 0]), 'count_varies')


class EndHintTests(unittest.TestCase):
    def test_inverted_repeat_and_direct_repeat(self):
        ir = 'GGTCTGACGCATTACGGATC'
        element = ir + 'A' * 50 + cr.revcomp(ir)
        seq = 'TTTTT' + 'CATGG' + element + 'CATGG' + 'AAAAA'
        tir, tsd = cr.end_hints(seq, 10, 10 + len(element))
        self.assertEqual(tir, 1.0)
        self.assertEqual(tsd, 5)


class AniTests(unittest.TestCase):
    def test_bins_are_descriptive(self):
        self.assertEqual([cr.ani_bin(x) for x in ('', 0.4, 1.0, 3.0, 12.5)],
                         ['no_cross_pair_with_ani', 'le_1', 'le_1', '1_to_5', 'gt_5'])


class ScanTests(unittest.TestCase):
    def test_mobile_neighbours_skipped_for_host_context(self):
        gs = contig('+++++++++')
        for g in gs:
            g['product'] = 'IS3 family transposase' if g['gene_uid'] in ('g5', 'g6') else 'host'
        labels = {g['gene_uid']: g['gene_uid'] for g in gs}
        mobile = lambda g: 'keyword' if 'transposase' in g['product'] else ''
        up, down = cr.scan(gs, 4, 4, True, labels, mobile, max_scan=10)
        self.assertEqual([x[2] for x in down], ['g5', 'g6', 'g7', 'g8'])
        self.assertEqual([x[2] for x in cr.host(down)], ['g7', 'g8'])        # only two host genes before the contig end
        self.assertEqual([x[2] for x in cr.host(up)], ['g3', 'g2', 'g1'])
        _, down = cr.scan(gs, 4, 4, True, labels, mobile, max_scan=2)        # scan cap reached before any host gene
        self.assertEqual(cr.host(down), [])

    def test_near_different_host_similar_is_nested(self):
        a = dict(up=[('X1', True, 'a1'), ('H1', True, 'a2'), ('H2', True, 'a3')], down=[('D1', True, 'a4')] * 3)
        b = dict(up=[('Y1', True, 'b1'), ('Y2', True, 'b2'), ('Y3', True, 'b3')], down=[('E1', True, 'b4')] * 3)
        self.assertEqual(cr.pc.relation(a, b), 'different_context_candidate')
        ha = dict(up=[('H1', True, 'a2'), ('H2', True, 'a3'), ('H3', True, 'a5')], down=[('D1', True, 'a4')] * 3)
        hb = dict(up=[('H1', True, 'b5'), ('H9', True, 'b6'), ('H8', True, 'b7')], down=[('D1', True, 'b8')] * 3)
        self.assertEqual(cr.pc.relation(ha, hb), 'similar_context')


if __name__ == '__main__':
    unittest.main()
