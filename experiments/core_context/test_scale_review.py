import random
import unittest

import core_review as cr

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

    def test_cross_assembly_pairs_are_kept_with_copy_class(self):
        # one copy in each of S1, S2, S3 and two in S4, every pair with a different host context on both sides except
        # a ~ c: all cross-assembly pairs are kept, classed by copy number, never filtered by it
        ys = [occ('a', 'S1', 'ABCDE', 'FGHIJ', 1000), occ('b', 'S2', 'PQRST', 'UVWXY', 1000),
              occ('c', 'S3', 'ABCDE', 'FGHIJ', 1000), occ('d1', 'S4', 'KLMNO', 'klmno', 1000),
              occ('d2', 'S4', 'abcde', 'fghij', 500000)]
        occ_rows, groups, sa, cross, rel, summary = sr.review_core(ys)
        cls = {(r['occurrence_a'], r['occurrence_b']): r['copy_class'] for r in cross}
        self.assertEqual(cls, {('a', 'b'): 'single_x_single', ('b', 'c'): 'single_x_single',
                               ('a', 'd1'): 'single_x_multi', ('a', 'd2'): 'single_x_multi',
                               ('b', 'd1'): 'single_x_multi', ('b', 'd2'): 'single_x_multi',
                               ('c', 'd1'): 'single_x_multi', ('c', 'd2'): 'single_x_multi'})
        self.assertEqual((summary['cross_different_position_pairs'], summary['cross_single_x_single'],
                          summary['cross_single_x_multi'], summary['cross_multi_x_multi']), (8, 2, 6, 0))
        self.assertEqual([(r['occurrence_a'], r['occurrence_b']) for r in sa if r['ucd_structural_pair']], [('d1', 'd2')])


    def test_two_similar_copies_in_one_assembly_differ_from_a_single_copy_elsewhere(self):
        # assembly A: two copies with similar contexts (no same-assembly U-core-D pair); assembly B: one copy with a
        # different context. Both A-B pairs must be kept as cross-assembly pairs (single_x_multi).
        ys = [occ('a1', 'A', 'ABCDE', 'FGHIJ', 1000), occ('a2', 'A', 'ABCDE', 'FGHIJ', 800000),
              occ('b', 'B', 'PQRST', 'UVWXY', 1000)]
        occ_rows, groups, sa, cross, rel, summary = sr.review_core(ys)
        self.assertEqual([r['ucd_structural_pair'] for r in sa], [False])
        self.assertEqual(sa[0]['context_class'], 'context_repeated_both_sides')
        self.assertEqual({(r['occurrence_a'], r['occurrence_b'], r['copy_class'], r['copies_a'], r['copies_b'])
                          for r in cross}, {('a1', 'b', 'single_x_multi', 2, 1), ('a2', 'b', 'single_x_multi', 2, 1)})

    def test_multi_copy_on_both_sides(self):
        ys = [occ('a1', 'A', 'ABCDE', 'FGHIJ', 1000), occ('a2', 'A', 'ABCDE', 'FGHIJ', 800000),
              occ('b1', 'B', 'PQRST', 'UVWXY', 1000), occ('b2', 'B', 'PQRST', 'UVWXY', 600000)]
        occ_rows, groups, sa, cross, rel, summary = sr.review_core(ys)
        self.assertEqual(len(cross), 4)
        self.assertEqual({r['copy_class'] for r in cross}, {'multi_x_multi'})
        self.assertFalse(any(r['ucd_structural_pair'] for r in sa))


class BoundaryReuseTests(unittest.TestCase):
    def test_reaches_end_recovered_from_cached_extensions(self):
        rnd = random.Random(4)
        for _ in range(500):
            alen, blen = rnd.randint(3000, 7000), rnd.randint(3000, 7000)
            a0, b0 = rnd.randint(0, 3000), rnd.randint(0, 3000)
            a1, b1 = min(alen, a0 + rnd.randint(300, 1500)), min(blen, b0 + rnd.randint(300, 1500))
            qs, ss = rnd.randint(0, a0), rnd.randint(0, b0)
            qe, se_ = rnd.randint(a1, alen), rnd.randint(b1, blen)
            pieces = [(qs, qe, ss, se_)] if rnd.random() < 0.7 else [(qs, a1, ss, b1)]
            ext = cr.extensions(pieces, a0, a1, alen, b0, b1, blen)
            written = {k: '' if ext[k] is None else str(ext[k]) for k in ('a_up', 'b_up', 'a_down', 'b_down')}
            row = dict(identical_protein='False', protein_identity_min='', protein_coverage_min='',
                       identical_core_dna='False', core_dna_identity='', core_dna_coverage_min='',
                       up_cover=ext['up_cover'], down_cover=ext['down_cover'], **written)
            v = sr.from_cache(row)
            A, B = dict(c0=a0, c1=a1, ext='N' * alen), dict(c0=b0, c1=b1, ext='N' * blen)
            self.assertEqual(sr.reaches_end(v, A, B), (ext['up_reaches_end'], ext['down_reaches_end']))
            self.assertEqual([v[k] for k in ('a_up', 'b_up', 'a_down', 'b_down')],
                             [ext[k] for k in ('a_up', 'b_up', 'a_down', 'b_down')])

    def test_core_match_class(self):
        self.assertEqual(sr.core_match_class(True, 100.0, 1.0), 'identical_core_dna')
        self.assertEqual(sr.core_match_class(False, 99.0, 0.95), 'homolog_full_coverage')
        self.assertEqual(sr.core_match_class(False, 99.0, 0.5), 'partial_core_match')
        self.assertEqual(sr.core_match_class(False, '', 0), 'no_core_dna_hit')


if __name__ == '__main__':
    unittest.main()
