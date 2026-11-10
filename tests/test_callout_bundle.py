import copy
import csv
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('bundle_check', ROOT/'scripts/check_callout_bundle.py')
checker = importlib.util.module_from_spec(spec)
spec.loader.exec_module(checker)


class BundleTests(unittest.TestCase):
    def setUp(self):
        self.s = ROOT/'experiments/read_callout/results/PaAR2038'
        c = self.s/'callout'
        with (self.s/'evidence/member_read_support.tsv').open() as f:
            evidence = list(csv.DictReader(f, delimiter='\t'))
        self.args = [checker.load(c/'summary.json'), checker.load(c/'validate.json'),
                     checker.load(c/'all_families_annotated.json'),
                     [l.split('\t') for l in (c/'candidates.bed').read_text().splitlines()],
                     evidence,
                     checker.load(self.s/'evidence/summary.json')]

    def test_all_four_archives_reconcile_without_becoming_release_approval(self):
        for sample in ('PaAR2038','DdD110','AmHOT1A3','RoMoAcy1'):
            s = self.s.parent/sample
            r = checker.audit(s, s/'reads_verify.json', snapshot=True)
            self.assertEqual(r['status'], 'snapshot_tables_consistent')
            self.assertFalse(r['scale_release_approved'])
            self.assertEqual(r['audited_read_check'], 'not_performed')
            self.assertTrue(all(f['long_v1_full_ge2_complete_unambiguous_triplet_records'] is None for f in r['family_rows']))

    def test_wrong_coordinate_version_rejected(self):
        self.args[4][0]['start0'] = str(int(self.args[4][0]['start0']) + 1)
        with self.assertRaisesRegex(ValueError, 'coordinates differ'):
            checker.check_tables(*self.args)

    def test_missing_read_row_rejected(self):
        self.args[4].pop()
        with self.assertRaisesRegex(ValueError, 'Missing/extra'):
            checker.check_tables(*self.args)

    def test_duplicate_read_row_rejected(self):
        self.args[4].append(copy.deepcopy(self.args[4][0]))
        with self.assertRaisesRegex(ValueError, 'Duplicate evidence'):
            checker.check_tables(*self.args)

    def test_tnp_subset_cannot_pass_all_family_contract(self):
        self.args[0]['subset'] = 'tnp-supported'
        with self.assertRaisesRegex(ValueError, 'subset all'):
            checker.check_tables(*self.args)

    def test_wrong_download_sample_rejected(self):
        data = checker.load(self.s/'reads_verify.json')
        data['run']['sample_id'] = 'another_isolate'
        with tempfile.TemporaryDirectory() as d:
            path = Path(d)/'verify.json'
            path.write_text(json.dumps(data))
            with self.assertRaisesRegex(ValueError, 'different sample'):
                checker.audit(self.s, path, snapshot=True)

    def test_v1_support_does_not_resolve_censored_boundaries_or_ns(self):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        e = [copy.deepcopy(self.args[4][0])]
        m = e[0]['member']
        b = {m: {'flanks_complete':True}}
        triplets = {m:SeqRecord(Seq('ACGT'))}
        e[0]['boundary_source'] = 'homology_extent_censored_at_window'
        self.assertEqual(checker.summarize(e,b,triplets)[0]['long_v1_full_ge2_complete_unambiguous_triplet_records'],0)
        e[0]['boundary_source'] = 'provisional_homology_median'
        triplets[m].seq = Seq('ACNT')
        self.assertEqual(checker.summarize(e,b,triplets)[0]['long_v1_full_ge2_complete_unambiguous_triplet_records'],0)
        triplets[m].seq = Seq('ACGT')
        self.assertEqual(checker.summarize(e,b,triplets)[0]['long_v1_full_ge2_complete_unambiguous_triplet_records'],1)

    def test_ar2038_is5_record_regression(self):
        rs = [r for r in self.args[4] if r['family']=='PaAR2038_RS_R1']
        self.assertEqual(len(rs),11)
        self.assertEqual((min(int(r['length']) for r in rs), max(int(r['length']) for r in rs)),(1171,1173))
        self.assertEqual((min(int(r['full']) for r in rs), max(int(r['full']) for r in rs)),(145,181))
        self.assertEqual(len({(r['contig'],r['start0'],r['end0']) for r in rs}),11)

    def test_live_zero_hit_run_passes_but_changed_input_fails(self):
        from multisite_insertion.pipeline import run
        with tempfile.TemporaryDirectory() as d:
            p = Path(d)
            (p/'native').mkdir()
            (p/'genome.fa').write_text('>c\nACGT\n')
            for f in ('families.fasta','ranges.tsv','families.gff','library.fasta'):
                (p/'native'/f).write_text('')
            (p/'manifest.json').write_text(json.dumps({'schema_version':1,'samples':[
                {'id':'s','fasta':'genome.fa','native':{'RS':'native','PILER':'native'}}]}))
            run(p/'manifest.json', p/'callout', methods=['PILER','RS'], subset='all')
            (p/'evidence').mkdir()
            (p/'evidence/member_read_support.tsv').write_text('member\tfamily\n')
            (p/'evidence/summary.json').write_text(json.dumps({'sample':'s','members':0,'member_classes':{}}))
            (p/'verify.json').write_text(json.dumps({'status':'verified','run':{'sample_id':'s'},'checks':
                {'bytes_match':True,'md5_match':True,'base_count_match':True}}))
            r = checker.audit(p, p/'verify.json')
            self.assertEqual(r['status'],'live_exports_and_v1_tables_consistent')
            self.assertEqual(r['long_records'],0)
            self.assertFalse(r['scale_release_approved'])
            (p/'callout/inputs/s/genome.fasta').write_text('>c\nAAAA\n')
            with self.assertRaisesRegex(ValueError,'Staged input changed'):
                checker.audit(p,p/'verify.json')


if __name__ == '__main__':
    unittest.main()
