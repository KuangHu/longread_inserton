import unittest

import scale_candidates as sc


class EvidenceStateTests(unittest.TestCase):
    def test_layers_are_kept_apart(self):
        self.assertEqual(sc.evidence_state('identical_core_dna', 'lt2500', 'lt2500'), 'match_both_lt2500')
        self.assertEqual(sc.evidence_state('homolog_full_coverage', 'lt2500', 'lt2500'), 'match_both_lt2500')
        # one copy's span unknown: pending, even if the other is resolved short
        self.assertEqual(sc.evidence_state('identical_core_dna', 'lt2500', 'unknown'), 'match_boundary_unknown_pending')
        self.assertEqual(sc.evidence_state('identical_core_dna', 'ge2500', 'unknown'), 'match_boundary_unknown_pending')
        self.assertEqual(sc.evidence_state('identical_core_dna', 'lt2500', 'ge2500'), 'match_span_ge2500')
        # a local core match stays pending whatever the spans
        self.assertEqual(sc.evidence_state('partial_core_match', 'lt2500', 'lt2500'), 'partial_core_match_pending')
        self.assertEqual(sc.evidence_state('no_core_dna_hit', 'lt2500', 'lt2500'), 'no_core_dna_hit')


if __name__ == '__main__':
    unittest.main()
