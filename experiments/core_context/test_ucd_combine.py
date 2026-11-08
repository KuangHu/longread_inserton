"""Tests of the strict length rule and helpers of ucd_combine / ucd_batch (no blastn needed)."""
import gzip
import tempfile
import unittest
from pathlib import Path

import ucd_batch as ub
import ucd_combine as cb


def row(la, lb, iva=None, ivb=None):
    return {'a_x_interval': iva or '100-%d' % (100 + la), 'b_x_interval': ivb or '0-%d' % lb,
            'x_len_a': str(la), 'x_len_b': str(lb)}


class Strict(unittest.TestCase):
    def test_both_below_2500(self):
        self.assertEqual(cb.strict_len(row(2499, 1)), ([], []))

    def test_exactly_2500_fails(self):
        self.assertEqual(cb.strict_len(row(2500, 900)), (['a'], []))
        self.assertEqual(cb.strict_len(row(2600, 2500)), (['a', 'b'], []))

    def test_length_from_interval_must_match_sequence(self):
        self.assertEqual(cb.strict_len(row(900, 900, iva='0-901')), ([], ['a']))
        self.assertEqual(cb.strict_len({'a_x_interval': '', 'x_len_a': '', 'b_x_interval': '0-5', 'x_len_b': '5'}),
                         ([], ['a']))


class Helpers(unittest.TestCase):
    def test_ids_of_table_and_list(self):
        with tempfile.TemporaryDirectory() as t:
            p = Path(t) / 'x.tsv.gz'
            with gzip.open(p, 'wt') as f:
                f.write('structure\tpair_id\nsingle\tM1\nsingle\tM2\n')
            self.assertEqual(ub.ids_of(p), {'M1', 'M2'})
            q = Path(t) / 'ids.txt'
            q.write_text('pair_id\nM3\nM4\n')
            self.assertEqual(ub.ids_of(q), {'M3', 'M4'})

    def test_candidate_interval_reason_is_explicit(self):
        r = {'a_input_error': '', 'a_boundary_up': 'end_not_covered', 'a_boundary_down': 'resolved', 'a_elem_lo': '',
             'a_start0': '10', 'a_end0': '20'}
        self.assertEqual(ub.interval(r, 'a', 'cds_core'), (10, 20))
        self.assertEqual(ub.interval(r, 'a', 'candidate_insertion'),
                         'candidate_interval_not_available(boundary_up=end_not_covered,pair_extension_not_determined)')


class Fingerprint(unittest.TestCase):
    R = {'a_sample': 'S1', 'a_assembly': 'G1', 'a_contig': 'c', 'a_strand': '+', 'a_input_error': '',
         'a_start0': '10', 'a_end0': '20', 'a_boundary_up': 'resolved', 'a_boundary_down': 'resolved', 'a_elem_lo': '5',
         'a_elem_hi': '25', 'b_sample': 'S2', 'b_assembly': 'G2', 'b_contig': 'd', 'b_strand': '-', 'b_input_error': '',
         'b_start0': '1', 'b_end0': '9', 'b_boundary_up': 'resolved', 'b_boundary_down': 'end_not_covered',
         'b_elem_lo': '', 'b_elem_hi': ''}
    MAN = {'S1': {'fna_sha256': 'x1'}, 'S2': {'fna_sha256': 'x2'}}
    CODE = {'ucd_compare.py': 'c1', 'ucd_checks.py': 'k1'}

    def fp(self, r=None, man=None, code=None, cfg=None, kind='cds_core'):
        return ub.fingerprint(r or self.R, kind, man or self.MAN, cfg or {'context': 1000}, {'m': 20},
                              code or self.CODE)

    def test_same_input_same_fingerprint(self):
        self.assertEqual(self.fp(), self.fp(dict(self.R)))

    def test_any_input_change_changes_it(self):
        base = self.fp()
        self.assertNotEqual(base, self.fp(dict(self.R, a_strand='-')))
        self.assertNotEqual(base, self.fp(dict(self.R, a_end0='21')))
        self.assertNotEqual(base, self.fp(man={'S1': {'fna_sha256': 'y'}, 'S2': {'fna_sha256': 'x2'}}))
        self.assertNotEqual(base, self.fp(code={'ucd_compare.py': 'c2', 'ucd_checks.py': 'k1'}))
        self.assertNotEqual(base, self.fp(cfg={'context': 500}))

    def test_candidate_depends_on_boundary_status(self):
        c = self.fp(kind='candidate_insertion')
        self.assertNotEqual(c, self.fp(dict(self.R, b_boundary_down='resolved', b_elem_lo='0', b_elem_hi='12'),
                                       kind='candidate_insertion'))
        self.assertEqual(self.fp(), self.fp(dict(self.R, b_boundary_down='resolved')))  # cds_core unaffected


if __name__ == '__main__':
    unittest.main()
