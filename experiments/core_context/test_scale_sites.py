"""Tests of the physical-site layer (scale_sites.py): several seeds at one place are one site."""
import unittest

import scale_sites as ss


def pos(ck, cls, lp, rp, mol=5, lc='c1', rc='c1', boundary='no_extension_detected', contig='c1', start='', end='',
        shifted=''):
    return dict(core_key=ck, cls=cls, left_contig=lc, left_pos=str(lp), right_contig=rc, right_pos=str(rp),
                molecules=str(mol), boundary=boundary, insert_len='1200', insert_size='lt_2500', catalog='',
                contig=contig, start=str(start), end=str(end), shifted_from=shifted)


class SiteTests(unittest.TestCase):
    def test_homologous_seeds_at_one_place_are_one_site(self):
        g = ss.group_sites([pos('PF_a', 'not_in_assembly:core_absent', 1000, 1010),
                            pos('PF_b', 'not_in_assembly:core_absent', 1030, 1000),
                            pos('PF_a', 'not_in_assembly:core_absent', 9000, 9010)])
        self.assertEqual(sorted(len(x) for x in g), [1, 2])

    def test_swapped_ends_and_tolerance(self):
        g = ss.group_sites([pos('PF_a', 'known_catalog_position', 1000, 2200),
                            pos('PF_b', 'known_catalog_position_other_core', 2210, 1020)])
        self.assertEqual(len(g), 1)
        g = ss.group_sites([pos('PF_a', 'known_catalog_position', 1000, 2200),
                            pos('PF_b', 'known_catalog_position', 1000, 2200 + 51)])
        self.assertEqual(len(g), 2)  # one end beyond TOL: not the same site

    def test_shifted_joins_its_known_position_and_category_order(self):
        rows = [pos('PF_a', 'known_catalog_position', 1000, 2200, contig='c1', start=1000, end=2200),
                pos('PF_a', 'shifted_known_position', 1400, 2200, shifted='c1:1000-2200:PF_a'),
                pos('PF_b', 'not_in_assembly:junction', 1000, 2200, mol=4)]
        g = ss.group_sites(rows)
        self.assertEqual(len(g), 1)
        r = ss.site_row('s:S0', 's', g[0])
        self.assertEqual((r['category'], r['cores'], 'cores_disagree' in r['flags']), ('known', 2, True))

    def test_below_threshold_flag(self):
        r = ss.site_row('s:S0', 's', [pos('PF_a', 'not_in_assembly:core_absent', 1, 5, mol=2)])
        self.assertIn('below_molecule_threshold', r['flags'])


if __name__ == '__main__':
    unittest.main()
