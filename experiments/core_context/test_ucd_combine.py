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


if __name__ == '__main__':
    unittest.main()
