import unittest

import read_site_search as rs


def aln(t='c', ts=0, te=1000, strand='+', qs=0, qe=1000, qlen=1000, score=1000):
    return dict(t=t, ts=ts, te=te, strand=strand, qs=qs, qe=qe, qlen=qlen, score=score)


FAM = {'PF_a': {'PH_1'}, 'PF_b': {'PH_1'}, 'PF_c': {'PH_2'}}


class PlacementTests(unittest.TestCase):
    def test_inner_end_is_next_to_the_core(self):
        self.assertEqual(rs.placement([aln(ts=100, te=1100)], 'L'), ('placed', 'c', 1100, '+'))
        self.assertEqual(rs.placement([aln(ts=2000, te=3000)], 'R'), ('placed', 'c', 2000, '+'))
        self.assertEqual(rs.placement([aln(ts=100, te=1100, strand='-')], 'L'), ('placed', 'c', 100, '-'))
        self.assertEqual(rs.placement([aln(ts=2000, te=3000, strand='-')], 'R'), ('placed', 'c', 3000, '-'))

    def test_competing_placement_is_ambiguous_and_short_alignment_unplaced(self):
        two = [aln(ts=100, te=1100, score=1000), aln(t='d', ts=5000, te=6000, score=990)]
        self.assertEqual(rs.placement(two, 'L')[0], 'ambiguous')
        far = [aln(ts=100, te=1100, score=1000), aln(t='d', ts=5000, te=6000, score=900)]
        self.assertEqual(rs.placement(far, 'L')[0], 'placed')
        self.assertEqual(rs.placement([aln(qs=0, qe=500)], 'L')[0], 'unplaced')
        self.assertEqual(rs.placement([], 'R')[0], 'unplaced')

    def test_inner_end_must_be_aligned(self):
        # counterexample to v0.1: 85% of the flank aligned, but not its inner end -> the inner coordinate is not seen
        self.assertEqual(rs.placement([aln(qs=0, qe=850, ts=100, te=950)], 'L')[0], 'inner_end_unaligned')
        self.assertEqual(rs.placement([aln(qs=150, qe=1000, ts=2000, te=2850)], 'R')[0], 'inner_end_unaligned')
        self.assertEqual(rs.placement([aln(qs=0, qe=850, ts=100, te=950, strand='-')], 'L')[0], 'inner_end_unaligned')
        # the outer end unaligned is fine; the inner end within INNER_CLIP is fine
        self.assertEqual(rs.placement([aln(qs=150, qe=1000, ts=250, te=1100)], 'L'), ('placed', 'c', 1100, '+'))
        self.assertEqual(rs.placement([aln(qs=0, qe=850, ts=2000, te=2850)], 'R'), ('placed', 'c', 2000, '+'))
        self.assertEqual(rs.placement([aln(qs=0, qe=990, ts=110, te=1100)], 'L')[0], 'placed')
        self.assertEqual(rs.placement([aln(qs=0, qe=979, ts=121, te=1100)], 'L')[0], 'inner_end_unaligned')


class GeometryTests(unittest.TestCase):
    def test_forward_and_reverse_give_the_sequence_between_the_flanks(self):
        self.assertEqual(rs.geometry('c', 1100, '+', 'c', 2100, '+', 1000), ('proper', 1100, 2100))
        self.assertEqual(rs.geometry('c', 2100, '-', 'c', 1100, '-', 1000), ('proper', 1100, 2100))
        self.assertEqual(rs.geometry('c', 1100, '+', 'c', 1100, '+', 1000), ('proper', 1100, 1100))  # empty

    def test_junctions(self):
        self.assertEqual(rs.geometry('c', 1100, '+', 'd', 2100, '+', 1000)[0], 'junction')
        self.assertEqual(rs.geometry('c', 1100, '+', 'c', 2100, '-', 1000)[0], 'junction')
        self.assertEqual(rs.geometry('c', 1100, '+', 'c', 50000, '+', 1000)[0], 'junction')
        self.assertEqual(rs.geometry('c', 2100, '+', 'c', 1100, '+', 1000)[0], 'junction')


class CanonicalTests(unittest.TestCase):
    def test_both_strands_of_one_position_give_the_same_key(self):
        fwd = dict(left_contig='c', left_pos=1100, left_strand='+', right_contig='c', right_pos=2100, right_strand='+')
        rev = dict(left_contig='c', left_pos=2100, left_strand='-', right_contig='c', right_pos=1100, right_strand='-')
        self.assertEqual(rs.canonical(fwd), rs.canonical(rev))
        self.assertEqual(rs.geometry(*rs.canonical(rev), 1000), ('proper', 1100, 2100))


class ClusterTests(unittest.TestCase):
    def test_observations_within_tolerance_merge(self):
        obs = [dict(c_left_pos=1000, c_right_pos=2000), dict(c_left_pos=1030, c_right_pos=2010),
               dict(c_left_pos=1500, c_right_pos=2500)]
        self.assertEqual(sorted(len(c) for c in rs.cluster(obs)), [1, 2])


def row(cls, ck='PF_a', cs='+', contig='c', start=1000, end=2000, **kw):
    return dict(dict(cls=cls, core_key=ck, core_strand=cs, contig=contig, start=start, end=end), **kw)


class ShiftedTests(unittest.TestCase):
    def test_position_overlapping_a_known_one_is_marked_and_keeps_its_class(self):
        rows = [row('known_catalog_position'),
                row('not_in_assembly:other_sequence', start=1150),
                row('catalog_omission', start=5000, end=6000),
                row('not_in_assembly:junction', contig='', start='', end='', left_contig='c', left_pos=1990,
                    right_contig='d', right_pos=10)]
        out = rs.mark_shifted(rows, FAM)
        self.assertEqual([r['cls'] for r in out], ['known_catalog_position', 'shifted_known_position',
                                                    'catalog_omission', 'shifted_known_position'])
        self.assertEqual(out[1]['raw_cls'], 'not_in_assembly:other_sequence')

    def test_core_relation_is_checked(self):
        # counterexamples to v0.1: an unrelated core's known position, or a homolog in the opposite orientation,
        # does not make an overlapping position shifted; a homologous core in a consistent orientation does
        rows = [row('known_catalog_position', ck='PF_c'),
                row('not_in_assembly:other_sequence', ck='PF_a', start=1150),
                row('known_catalog_position', ck='PF_b', contig='d'),
                row('catalog_omission', ck='PF_a', contig='d', start=1100, cs='-'),
                row('catalog_omission', ck='PF_a', contig='d', start=1100, end=1900)]
        out = rs.mark_shifted(rows, FAM)
        self.assertEqual([r['cls'] for r in out], ['known_catalog_position', 'not_in_assembly:other_sequence',
                                                    'known_catalog_position', 'catalog_omission',
                                                    'shifted_known_position'])


def occ(oid, ck='PF_a', start=1000, end=2000, strand='+', contig='c'):
    return dict(occurrence_id=oid, core_key=ck, contig=contig, start0=start, end0=end, strand=strand)


class CoreRelationTests(unittest.TestCase):
    def test_relation_and_sign(self):
        self.assertEqual(rs.core_relation('PF_a', 'PF_a', FAM), ('same', '+'))
        self.assertEqual(rs.core_relation('PF_a', 'PF_b', FAM), ('homolog', '+'))  # same homology group
        self.assertEqual(rs.core_relation('PF_a', 'PF_c', FAM), ('unrelated', ''))
        self.assertEqual(rs.core_relation('PF_a+,PF_c+', 'PF_c-,PF_x+', FAM), ('homolog', '-'))
        self.assertEqual(rs.core_relation('PF_a+,PF_c+', 'PF_a+,PF_c-', FAM), ('homolog', ''))  # pairs disagree
        self.assertEqual(rs.core_relation('NEG_recA', 'PF_a', FAM), ('unrelated', ''))


class MatchCatalogTests(unittest.TestCase):
    def test_same_core_orientation_and_coverage(self):
        # boundary offset of 60 bp (> TOL; v0.1 containment failed) but reciprocal coverage holds -> matched
        self.assertEqual(rs.match_catalog('c', 1000, 2000, '+', 'PF_a', [occ('o1', start=940, end=2000)], FAM)[:2],
                         ('known_catalog_position', 'o1'))
        # counterexample: same core, opposite orientation -> not known
        self.assertEqual(rs.match_catalog('c', 1000, 2000, '-', 'PF_a', [occ('o1')], FAM)[0],
                         'catalog_mismatch:orientation')

    def test_overlap_alone_is_not_a_match(self):
        # an unrelated core wholly inside the interval (v0.1 called this known_catalog_position_other_core)
        cls, ids, notes = rs.match_catalog('c', 1000, 2000, '+', 'PF_a', [occ('o2', ck='PF_c', start=1100, end=1900)], FAM)
        self.assertEqual((cls, ids, notes), (None, '', 'o2:unrelated'))
        # a homolog covering only half the interval (a one-gene occurrence inside a two-gene core): coverage fails
        self.assertEqual(rs.match_catalog('c', 1000, 3000, '+', 'PF_a+,PF_c+', [occ('o3', ck='PF_b', end=2000)],
                                          FAM)[0], 'catalog_mismatch:coverage')
        # the same core's occurrence overlapping by 30%: not a match either
        self.assertEqual(rs.match_catalog('c', 1000, 2000, '+', 'PF_a', [occ('o4', start=1700, end=2700)], FAM)[0],
                         'catalog_mismatch:coverage')
        # same core preferred over a matching homolog
        self.assertEqual(rs.match_catalog('c', 1000, 2000, '+', 'PF_a', [occ('o5', ck='PF_b'), occ('o6')], FAM)[:2],
                         ('known_catalog_position', 'o6'))
        self.assertEqual(rs.match_catalog('c', 1000, 2000, '+', 'PF_a', [occ('o5', ck='PF_b')], FAM)[:2],
                         ('known_catalog_position_other_core', 'o5'))


def galn(ts, te, cg, t='c'):
    return dict(t=t, ts=ts, te=te, cg=cg, strand='+', qs=0, qe=te - ts)


class EmptySiteTests(unittest.TestCase):
    def test_local_cigar(self):
        clean = galn(3000, 7000, '4000M')
        self.assertEqual(rs.empty_site_reads([clean], 5000, 5008, 1200), 'empty')
        # counterexample to v0.1 (length difference < 200 bp overall): a 1.2 kb insertion at the site and a 1.2 kb
        # deletion 1.5 kb away cancel in length, but the read carries an insertion at the site
        ins = galn(3000, 8200, '2000M1200I1500M1200D500M')
        self.assertEqual(rs.empty_site_reads([ins], 2000 + 3000, 2008 + 3000, 1200), 'indel')
        self.assertEqual(rs.empty_site_reads([galn(4500, 7000, '2500M')], 5000, 5008, 1200), '')  # does not span

    def test_molecules(self):
        mol = {'r1': 'm1', 'r1d': 'm1', 'r2': 'm2', 'r3': 'm3', 'r4': 'm4'}
        got = rs.empty_site_counts({'r3'}, {'r1': 'empty', 'r1d': 'empty', 'r2': 'indel', 'r3': 'empty',
                                            'r4': 'empty'}, mol)
        # r1 and its duplex record r1d are one molecule; r3 also carries the core there -> both_forms, not empty
        self.assertEqual(got, dict(empty_site_molecules=2, empty_site_read_records=3, empty_site_indel_molecules=1,
                                   both_forms_molecules=1))


if __name__ == '__main__':
    unittest.main()
