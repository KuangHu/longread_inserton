import copy
import random
import unittest

from core_context import extract, identity, oriented_segments, relation, sha
from Bio.Seq import Seq


def fixture(same_flanks=False, reverse=False):
    rng = random.Random(43)
    dna = lambda n: ''.join(rng.choices('ACGT', k=n))
    core = dna(800)
    first = [dna(500), core, dna(500)]
    second = first if same_flanks else [dna(500), core, dna(500)]
    if reverse:
        second = [str(Seq(s).reverse_complement()) for s in reversed(second)]
    graph = dict(paths={}, nodes={}, blocks={})
    genomes, mapping = {}, {}
    for i, pieces in enumerate([first, second]):
        name = 'sample' + str(i); seq = ''.join(pieces)
        genomes[name] = seq
        mapping[name] = dict(sample=name, contig='chr', length=len(seq), sha256=sha(seq))
        ids = [i*3+j for j in range(3)]
        graph['paths'][str(i)] = dict(id=i, name=name, tot_len=len(seq), circular=False, nodes=ids)
        start = 0
        for j, (nid, part) in enumerate(zip(ids, pieces)):
            bid = 100 if j == 1 else (j if same_flanks else i*10+j)
            graph['blocks'][str(bid)] = dict(id=bid)
            graph['nodes'][str(nid)] = dict(id=nid, block_id=bid, path_id=i,
                strand='-' if reverse and i == 1 else '+', position=[start, start+len(part)])
            start += len(part)
    return graph, genomes, mapping


class Tests(unittest.TestCase):
    def groups(self, **kwargs):
        return extract(*fixture(**kwargs), 500, 50, 2500)[0]

    def test_both_filled_no_empty_site(self):
        # Both assemblies contain the core. There are no gap calls or gap inputs.
        xs = self.groups()['PG_100']
        self.assertEqual(len(xs), 2)
        self.assertEqual(len({x['sample'] for x in xs}), 2)
        a, b = xs
        self.assertEqual(relation(identity(a['core'], b['core']), identity(a['upstream'], b['upstream']),
            identity(a['downstream'], b['downstream']), a['anchors'], b['anchors'], True),
            'different_adjacency_low_flank_similarity_candidate')
        self.assertTrue(all(not x['gap_required'] for x in xs))

    def test_same_context_not_new_site(self):
        a, b = self.groups(same_flanks=True)['PG_100']
        self.assertEqual(relation(1, 1, 1, a['anchors'], b['anchors'], True),
                         'similar_flanks_candidate_same_context')

    def test_reverse_strand(self):
        a, b = self.groups(reverse=True)['PG_100']
        self.assertEqual(a['core'], b['core'])
        self.assertEqual(len(b['upstream']), 500)
        self.assertEqual(b['anchors'], [('12', '+'), ('10', '+')])

    def test_low_identity_alone_not_different(self):
        self.assertEqual(relation(1, .5, .5, [('A','+'), ('B','+')], [('A','+'), ('B','+')], True),
                         'context_unresolved')

    def test_edge_not_negative_context(self):
        self.assertEqual(relation(1, None, .5, [None, None], [None, None], False), 'edge_or_missing_context')

    def test_core_mismatch_retained_as_diagnostic(self):
        self.assertEqual(relation(.7, .5, .5, [None, None], [None, None], True), 'core_similarity_below_threshold')

    def test_multiple_copies_same_sample_retained(self):
        g, f, m = fixture()
        m['sample1']['sample'] = 'sample0'
        xs = extract(g, f, m, 500, 50, 2500)[0]['PG_100']
        self.assertEqual(len(xs), 2)
        self.assertEqual(len({x['sample'] for x in xs}), 1)

    def test_strict_length_limit(self):
        g, f, m = fixture()
        self.assertNotIn('PG_100', extract(g, f, m, 500, 50, 800)[0])

    def test_input_hash_mismatch(self):
        g, f, m = fixture(); f['sample0'] = 'N' + f['sample0'][1:]
        with self.assertRaises(ValueError): extract(g, f, m, 500, 50, 2500)

    def test_circular_deferred(self):
        g, f, m = fixture(); g['paths']['0']['circular'] = True
        groups, excluded = extract(g, f, m, 500, 50, 2500)
        self.assertNotIn('PG_100', groups)
        self.assertEqual(len(excluded), 3)


if __name__ == '__main__':
    unittest.main()
