import gzip
import hashlib
import random
import shutil
import tempfile
import unittest
from argparse import Namespace
from pathlib import Path

import protein_context as pc
import scale_entry as se

MMSEQS = '/global/scratch/users/kh36969/longread_callout/env/s4/bin/mmseqs'
AA = 'ACDEFGHIKLMNPQRSTVWY'


def h(s):
    return hashlib.sha256(s.encode()).hexdigest()


class FamilyMappingTests(unittest.TestCase):
    """Batch 2 extends the batch-1 mapping: mapped proteins keep their ids, a close variant joins its family by the
    same member thresholds, an unrelated protein gets a new family."""

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        rnd = random.Random(11)
        prot = lambda n: 'M' + ''.join(rnd.choice(AA) for _ in range(n))
        self.p = [prot(300) for _ in range(4)]
        v = list(self.p[0])
        for i in rnd.sample(range(1, 300), 12):          # 96% identity, full length
            v[i] = rnd.choice([x for x in AA if x != v[i]])
        self.variant, self.unrelated = ''.join(v), prot(280)
        (self.tmp / 'genes').mkdir()
        for smp, seqs in (('S1', self.p), ('S2', [self.p[1], self.variant, self.unrelated])):
            with gzip.open(self.tmp / 'genes' / (smp + '.faa.gz'), 'wt') as f:
                for s in seqs:
                    f.write('>%s\n%s\n' % (h(s), s))
        pc.write_tsv(self.tmp / 'manifest.tsv', [dict(batch=1, sample='S1'), dict(batch=2, sample='S2')], ['batch', 'sample'])

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def test_incremental_mapping(self):
        a = lambda b: Namespace(manifest=self.tmp / 'manifest.tsv', work=self.tmp, batch=b, mmseqs=MMSEQS, threads=2, pool=False)
        se.cmd_families(a(1))
        b1 = dict(se.load_catalog(self.tmp))
        self.assertEqual(len(b1), 4)
        se.cmd_families(a(2))
        cat = se.load_catalog(self.tmp)
        for k, r in b1.items():
            self.assertEqual((cat[k]['family'], cat[k]['homology_group']), (r['family'], r['homology_group']))
        v, u = cat[h(self.variant)], cat[h(self.unrelated)]
        self.assertEqual(v['family'], cat[h(self.p[0])]['family'])
        self.assertEqual(v['assigned'], 'family:search,homology_group:search')
        self.assertEqual(u['assigned'], 'family:cluster,homology_group:cluster')
        self.assertNotIn(u['family'], {r['family'] for r in b1.values()})
        self.assertEqual(len(cat), 6)        # p[1] in batch 2 was already mapped and is not added twice


    def test_pool_mapping_is_one_clustering_of_every_assembly(self):
        for smp, seqs in (('S1', self.p), ('S2', [self.p[1], self.variant, self.unrelated])):
            with open(self.tmp / (smp + '.faa'), 'w') as f:
                for i, s in enumerate(seqs):
                    f.write('>%s_%d\n%s\n' % (smp, i, s))
        pc.write_tsv(self.tmp / 'pool.tsv', [dict(batch=1, sample='S1', faa=self.tmp / 'S1.faa'),
                                             dict(batch=2, sample='S2', faa=self.tmp / 'S2.faa')], ['batch', 'sample', 'faa'])
        w = self.tmp / 'pw'
        se.cmd_pool_proteins(Namespace(manifest=self.tmp / 'pool.tsv', work=w))
        se.cmd_families(Namespace(manifest=self.tmp / 'pool.tsv', work=w, batch=1, mmseqs=MMSEQS, threads=2, pool=True))
        cat = se.load_catalog(w)
        self.assertEqual(len(cat), 6)
        self.assertEqual(cat[h(self.variant)]['family'], cat[h(self.p[0])]['family'])
        self.assertEqual({r['assigned'] for r in cat.values()}, {'cluster'})
        with self.assertRaises(ValueError):          # the pool mapping is built once, never extended in place
            se.cmd_families(Namespace(manifest=self.tmp / 'pool.tsv', work=w, batch=1, mmseqs=MMSEQS, threads=2, pool=True))


class ShardTests(unittest.TestCase):
    def test_shard_is_stable_and_in_range(self):
        self.assertEqual(se.shard_of('single_gene', 'PF_abc'), se.shard_of('single_gene', 'PF_abc'))
        self.assertTrue(all(0 <= se.shard_of('two_gene', 'PF_%d' % i) < se.SHARDS for i in range(200)))


if __name__ == '__main__':
    unittest.main()
