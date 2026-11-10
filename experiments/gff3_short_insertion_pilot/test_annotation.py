"""Coordinate and annotation challenge cases, without network access."""
import unittest
from annotate_candidates import parse_gff, audit, attributes


def gff(kind, start, end, attrs, strand='+'):
    return f'c\ttest\t{kind}\t{start}\t{end}\t.\t{strand}\t0\t{attrs}\n'


def run(lines, start=100, end=300):
    features, _, _ = parse_gff(lines)
    return audit({'seed_id': 'test', 'filled_contig': 'c', 'start0': str(start),
                  'end0': str(end), 'length': str(end-start)}, features, {'p': 'M'})[0]


class AnnotationTests(unittest.TestCase):
    def test_exact_bounds_and_adjacent_gene(self):
        result = run([gff('gene', 101, 300, 'ID=g'),
                      gff('CDS', 101, 300, 'ID=d;Parent=g;protein_id=p'),
                      gff('gene', 301, 400, 'ID=adjacent')])
        self.assertEqual(result['annotation_status'], 'eligible_1_2_CDS')

    def test_one_base_cut_is_review(self):
        result = run([gff('gene', 100, 300, 'ID=g'),
                      gff('CDS', 100, 300, 'ID=d;Parent=g;protein_id=p')])
        self.assertEqual(result['complete_CDS_count'], 0)
        self.assertIn('boundary_intersects_feature', result['review_reasons'])

    def test_multipart_minus_strand_counts_once(self):
        result = run([gff('gene', 101, 300, 'ID=g', '-'),
                      gff('CDS', 101, 150, 'ID=d;Parent=g;protein_id=p', '-'),
                      gff('CDS', 201, 300, 'ID=d;Parent=g;protein_id=p', '-')])
        self.assertEqual(result['complete_CDS_count'], 1)
        self.assertEqual(result['annotation_status'], 'eligible_1_2_CDS')

    def test_parent_pseudogene_and_partial(self):
        result = run([gff('pseudogene', 101, 300, 'ID=g;pseudo=true;partial=true'),
                      gff('CDS', 101, 300, 'ID=d;Parent=g;protein_id=p')])
        self.assertIn('pseudogene', result['review_reasons'])
        self.assertIn('partial', result['review_reasons'])

    def test_missing_protein_and_rna(self):
        result = run([gff('gene', 101, 200, 'ID=g'),
                      gff('CDS', 101, 200, 'ID=d;Parent=g;protein_id=missing'),
                      gff('tRNA', 220, 290, 'ID=r')])
        self.assertIn('protein_sequence_missing', result['review_reasons'])
        self.assertIn('RNA_overlap', result['review_reasons'])

    def test_encoded_delimiters(self):
        self.assertEqual(attributes('ID=g;product=alpha%2Cbeta%3Bgamma')['product'], ['alpha,beta;gamma'])


if __name__ == '__main__':
    unittest.main()
