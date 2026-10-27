import unittest

import core_review as cr
import site_groups as sg


def win(labels, same=True):
    return [(l, same, 'u%d' % i) for i, l in enumerate(labels)]


def occ(sample, up, down, start=0):
    return dict(sample=sample, contig='c', start0=str(start), end0=str(start + 100), host=dict(up=win(up), down=win(down)))


class SideRelationTests(unittest.TestCase):
    def test_ordered_oriented_matches(self):
        self.assertEqual(sg.side_relation(win('ABCDE'), win('ABXDE'))[:2], ('same', 4))
        self.assertEqual(sg.side_relation(win('ABCDE'), win('PQRST'))[0], 'different')
        # shared labels in the opposite order give at most one ordered match: not same, and not different either
        self.assertEqual(sg.side_relation(win('ABPQR'), win('BAXYZ'))[:2], ('unknown', 1))
        # orientation must agree
        y = win('ABXYZ'); y[1] = ('B', False, 'u1')
        self.assertEqual(sg.side_relation(win('ABPQR'), y)[:2], ('unknown', 1))

    def test_missing_labels_are_unknown_not_different(self):
        self.assertEqual(sg.side_relation(win('ABCD') + [(None, True, 'u4')], win('PQRST'))[0], 'unknown')
        self.assertEqual(sg.side_relation(win('ABCD'), win('PQRS'))[0], 'unknown')      # short window (contig edge)

    def test_repeated_flank_label_never_matches(self):
        rel = sg.side_relation(win('AAPQR'), win('AAXYZ'))
        self.assertEqual(rel[0], 'unknown')
        self.assertTrue(rel[2])


class ReverseInstanceTests(unittest.TestCase):
    def test_minus_strand_copy_matches_plus_strand_copy(self):
        # A: L2 L1 [core+] R1 R2 ... on the + strand; B: the same locus written in reverse complement
        names = ['L5', 'L4', 'L3', 'L2', 'L1', 'core', 'R1', 'R2', 'R3', 'R4', 'R5']
        strands = ['+', '-', '+', '+', '-', '+', '+', '-', '-', '+', '+']
        ga = [dict(gene_uid='a%d' % i, strand=s, circular_contig='False') for i, s in enumerate(strands)]
        flip = {'+': '-', '-': '+'}
        gb = [dict(gene_uid='b%d' % i, strand=flip[s], circular_contig='False') for i, s in enumerate(strands[::-1])]
        labels = {g['gene_uid']: names[i] for i, g in enumerate(ga)}
        labels.update({g['gene_uid']: names[::-1][i] for i, g in enumerate(gb)})
        no = lambda g: ''
        ua, da = cr.scan(ga, 5, 5, True, labels, no, max_scan=10, keep=5)
        ub, db = cr.scan(gb, 5, 5, False, labels, no, max_scan=10, keep=5)
        self.assertEqual([x[0] for x in ub], [x[0] for x in ua])
        self.assertEqual([x[1] for x in ub], [x[1] for x in ua])
        self.assertEqual(sg.side_relation([x[:3] for x in ua], [x[:3] for x in ub])[:2], ('same', 5))
        self.assertEqual(sg.side_relation([x[:3] for x in da], [x[:3] for x in db])[:2], ('same', 5))


class GroupTests(unittest.TestCase):
    def run_group(self, os_):
        ids = sorted(os_)
        verdicts = {frozenset((x, y)): sg.pair_verdict(os_[x], os_[y]) for i, x in enumerate(ids) for y in ids[i + 1:]}
        return sg.group(ids, {x: os_[x]['sample'] for x in ids}, verdicts), verdicts

    def test_same_position_across_samples(self):
        os_ = {'a': occ('S1', 'ABCDE', 'FGHIJ'), 'b': occ('S2', 'ABCDE', 'FGHIJ'), 'c': occ('S3', 'ABXDE', 'FGHYJ')}
        (res, groups, amb), _ = self.run_group(os_)
        self.assertEqual(groups, [['a', 'b', 'c']])
        self.assertEqual({r[0] for r in res.values()}, {'site_group'})

    def test_one_strain_two_positions(self):
        os_ = {'a': occ('S1', 'ABCDE', 'FGHIJ'), 'b': occ('S1', 'PQRST', 'UVWXY', start=5000)}
        (res, groups, amb), v = self.run_group(os_)
        self.assertEqual(v[frozenset('ab')]['verdict'], 'different_position_same_assembly')
        self.assertEqual(groups, [])
        self.assertEqual({r[0] for r in res.values()}, {'single_observation_site'})

    def test_same_assembly_repeated_flank_is_flagged_and_not_merged(self):
        os_ = {'a': occ('S1', 'ABCDE', 'FGHIJ'), 'b': occ('S1', 'ABCDE', 'FGHIJ', start=9000)}
        (res, groups, amb), v = self.run_group(os_)
        self.assertTrue(v[frozenset('ab')]['repeated_flank'])
        self.assertEqual(groups, [])

    def test_two_copies_matching_one_position_are_ambiguous(self):
        # b1 and b2 sit in one assembly at different places, both matching a's flanks (segmental duplicate)
        os_ = {'a': occ('S1', 'ABCDE', 'FGHIJ'), 'b1': occ('S2', 'ABCDE', 'FGHIJ'),
               'b2': occ('S2', 'ABCDE', 'FGHIJ', start=9000)}
        (res, groups, amb), _ = self.run_group(os_)
        self.assertEqual(groups, [])
        self.assertEqual(res['a'], ('unresolved', 'ambiguous_multiple_correspondence'))

    def test_no_transitive_merge(self):
        # A ~ B (shares the inner genes), B ~ C (shares the outer genes), A vs C: no shared label -> conflict
        os_ = {'A': occ('S1', 'ABPQR', 'FGPQR'), 'B': occ('S2', 'ABCDE', 'FGHIJ'), 'C': occ('S3', 'XYZDE', 'KLMIJ')}
        (res, groups, amb), v = self.run_group(os_)
        self.assertEqual(v[frozenset('AB')]['verdict'], 'same_position_candidate')
        self.assertEqual(v[frozenset('BC')]['verdict'], 'same_position_candidate')
        self.assertNotEqual(v[frozenset('AC')]['verdict'], 'same_position_candidate')
        self.assertEqual(groups, [])
        self.assertEqual({r for r in res.values()}, {('unresolved', 'conflict_in_component')})

    def test_unresolved_is_not_a_site(self):
        os_ = {'a': occ('S1', 'ABCDE', 'FGHIJ'), 'b': occ('S2', 'XYZ', 'PQRST')}      # short up context on b: unknown, not different
        (res, groups, amb), _ = self.run_group(os_)
        self.assertEqual(res['a'], ('unresolved', 'unresolved_pairs'))
        self.assertEqual(res['b'], ('unresolved', 'unresolved_pairs'))

    def test_lis(self):
        self.assertEqual(sg.lis([0, 2, 1, 3, 4]), 4)
        self.assertEqual(sg.lis([]), 0)


if __name__ == '__main__':
    unittest.main()
