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
    def test_same_assembly_ucd_pair_needs_supported_positions_and_no_physical_overlap(self):
        # S1 holds two copies at different positions; S2 and S3 repeat each position (site groups), so both S1 copies
        # are supported positions with host context different on both sides
        ys = [occ('a1', 'S1', 'ABCDE', 'FGHIJ', 1000), occ('a2', 'S1', 'PQRST', 'UVWXY', 900000),
              occ('b1', 'S2', 'ABCDE', 'FGHIJ', 5000), occ('b2', 'S2', 'PQRST', 'UVWXY', 700000),
              occ('c1', 'S3', 'ABCDE', 'FGHIJ', 3000)]
        occ_rows, groups, sa, rel, summary = sr.review_core(ys)
        self.assertEqual(summary['site_groups'], 2)
        ucd = [r for r in sa if r['ucd_pair_candidate']]
        self.assertEqual({(r['occurrence_a'], r['occurrence_b']) for r in ucd}, {('a1', 'a2'), ('b1', 'b2')})
        self.assertEqual(summary['ucd_assemblies'], 2)
        self.assertEqual({r['relation'] for r in rel}, {'different'})

    def test_interlocking_windows_are_not_ucd(self):
        # two nearby copies: the genes between them are x's downstream and y's upstream (same gene ids)
        x = occ('x', 'S1', 'ABCDE', 'KLMNO', 1000)
        y = occ('y', 'S1', 'VWXYZ', 'PQRST', 4000)
        y['host']['up'] = [('K2', True, 'xd0'), ('L2', True, 'xd1')] + y['host']['up'][:3]
        occ_rows, groups, sa, rel, summary = sr.review_core([x, y])
        (r,) = sa
        self.assertEqual(r['context_class'], 'context_different_both_sides')
        self.assertEqual(r['host_genes_shared_physically'], 2)
        self.assertFalse(r['ucd_pair_candidate'])

    def test_unlabelled_neighbours_stay_unresolved_not_dropped(self):
        a = occ('a', 'S1', 'ABCDE', 'FGHIJ', 1000)
        b = occ('b', 'S2', 'PQRST', 'UVWXY', 1000)
        b['host']['up'][0] = (None, True, 'bu0')
        occ_rows, groups, sa, rel, summary = sr.review_core([a, b])
        self.assertEqual({r['s6_status'] for r in occ_rows}, {'unresolved'})
        self.assertEqual(summary['unresolved'], 2)


if __name__ == '__main__':
    unittest.main()
