import unittest

import scale_review as sr


def win(labels, tag):
    return [(l, True, '%s%d' % (tag, i)) for i, l in enumerate(labels)]


def occ(oid, sample, up, down, start, contig='c', uids=None):
    u = uids or oid
    return dict(occurrence_id=oid, sample=sample, contig=contig, start0=str(start), end0=str(start + 900), strand='+',
                host=dict(up=win(up, u + 'u'), down=win(down, u + 'd')), near_mobile_keyword='0',
                skipped_mobile_up='0', skipped_mobile_down='0')


class ReviewCoreTests(unittest.TestCase):
    def test_same_assembly_ucd_pair_with_supported_positions(self):
        # S1 holds two copies at different positions; S2 and S3 repeat each position (site groups), so both S1 copies
        # are supported positions with host context different on both sides
        ys = [occ('a1', 'S1', 'ABCDE', 'FGHIJ', 1000), occ('a2', 'S1', 'PQRST', 'UVWXY', 900000),
              occ('b1', 'S2', 'ABCDE', 'FGHIJ', 5000), occ('b2', 'S2', 'PQRST', 'UVWXY', 700000),
              occ('c1', 'S3', 'ABCDE', 'FGHIJ', 3000)]
        occ_rows, groups, sa, cross, rel, summary = sr.review_core(ys)
        self.assertEqual(summary['site_groups'], 2)
        ucd = [r for r in sa if r['ucd_structural_pair']]
        self.assertEqual({(r['occurrence_a'], r['occurrence_b']) for r in ucd}, {('a1', 'a2'), ('b1', 'b2')})
        self.assertEqual({r['s6_grouping'] for r in ucd}, {'supported_positions'})
        self.assertTrue(all(r['ucd_pair_candidate_strict_v0'] for r in ucd))
        self.assertEqual(summary['ucd_assemblies'], 2)
        self.assertEqual({r['relation'] for r in rel}, {'different'})

    def test_interlocking_windows_are_not_ucd(self):
        # two nearby copies: the genes between them are x's downstream and y's upstream (same gene ids)
        x = occ('x', 'S1', 'ABCDE', 'KLMNO', 1000)
        y = occ('y', 'S1', 'VWXYZ', 'PQRST', 4000)
        y['host']['up'] = [('K2', True, 'xd0'), ('L2', True, 'xd1')] + y['host']['up'][:3]
        occ_rows, groups, sa, cross, rel, summary = sr.review_core([x, y])
        (r,) = sa
        self.assertEqual(r['context_class'], 'context_different_both_sides')
        self.assertEqual(r['host_genes_shared_physically'], 2)
        self.assertFalse(r['ucd_structural_pair'])
        self.assertFalse(r['ucd_pair_candidate_strict_v0'])

    def test_unlabelled_neighbours_stay_unresolved_not_dropped(self):
        a = occ('a', 'S1', 'ABCDE', 'FGHIJ', 1000)
        b = occ('b', 'S2', 'PQRST', 'UVWXY', 1000)
        b['host']['up'][0] = (None, True, 'bu0')
        occ_rows, groups, sa, cross, rel, summary = sr.review_core([a, b])
        self.assertEqual({r['s6_status'] for r in occ_rows}, {'unresolved'})
        self.assertEqual(summary['unresolved'], 2)


    def test_other_occurrences_change_s6_grouping_but_not_the_pair_structure(self):
        # the pair (a1, a2): sequences, family and local contexts fixed. Adding c (S3) and f (S4) chains a1's site group
        # into a component that is not a clique (c ~ a1, f ~ c, f vs a1 conflicting), so a1 and b1 become unresolved.
        # That changes the S6 grouping field and the strict v0 flag only, never the local U-core-D judgement.
        base = [occ('a1', 'S1', 'ABCDE', 'FGHIJ', 1000), occ('a2', 'S1', 'PQRST', 'UVWXY', 900000),
                occ('b1', 'S2', 'ABCDE', 'FGHIJ', 5000), occ('b2', 'S2', 'PQRST', 'UVWXY', 700000)]
        extra = [occ('c', 'S3', 'DEKLM', 'FGHIJ', 3000), occ('f', 'S4', 'KLMVZ', 'FGHIJ', 3000)]
        got = {}
        for name, ys in (('before', base), ('after', base + extra)):
            occ_rows, groups, sa, cross, rel, summary = sr.review_core(ys)
            got[name] = ({r['occurrence_id']: (r['s6_status'], r['unresolved_reason']) for r in occ_rows},
                         {(r['occurrence_a'], r['occurrence_b']): r for r in sa})
        (st0, p0), (st1, p1) = got['before'], got['after']
        self.assertEqual(st0['a1'], ('site_group', ''))
        self.assertEqual(st1['a1'], ('unresolved', 'conflict_in_component'))
        for k in (('a1', 'a2'), ('b1', 'b2')):
            for f in ('context_class', 'up_support', 'down_support', 'host_genes_shared_physically', 'ucd_structural_pair'):
                self.assertEqual(p0[k][f], p1[k][f], (k, f))
            self.assertTrue(p1[k]['ucd_structural_pair'])
            self.assertEqual((p0[k]['s6_grouping'], p1[k]['s6_grouping']), ('supported_positions', 'grouping_unresolved'))
            self.assertEqual((p0[k]['ucd_pair_candidate_strict_v0'], p1[k]['ucd_pair_candidate_strict_v0']), (True, False))
        self.assertEqual(p1[('a1', 'a2')]['s6_status_a'], 'unresolved:conflict_in_component')

    def test_cross_assembly_single_copy_pairs_are_kept(self):
        # one copy in each of S1, S2, S3, with a different host context on both sides: no same-assembly pair exists, the
        # cross-assembly pairs carry the candidate. S4 holds two copies, so its pairs are left to the same-assembly scope.
        ys = [occ('a', 'S1', 'ABCDE', 'FGHIJ', 1000), occ('b', 'S2', 'PQRST', 'UVWXY', 1000),
              occ('c', 'S3', 'ABCDE', 'FGHIJ', 1000), occ('d1', 'S4', 'KLMNO', 'klmno', 1000),
              occ('d2', 'S4', 'abcde', 'fghij', 500000)]
        occ_rows, groups, sa, cross, rel, summary = sr.review_core(ys)
        self.assertEqual({(r['occurrence_a'], r['occurrence_b']) for r in cross}, {('a', 'b'), ('b', 'c')})
        self.assertEqual({r['s6_grouping'] for r in cross}, {'supported_positions'})
        self.assertEqual(summary['cross_single_copy_pairs'], 2)
        self.assertGreater(summary['cross_different_position_pairs'], 2)
        self.assertEqual([(r['occurrence_a'], r['occurrence_b']) for r in sa if r['ucd_structural_pair']], [('d1', 'd2')])


if __name__ == '__main__':
    unittest.main()
