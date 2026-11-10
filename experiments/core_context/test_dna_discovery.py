"""Unit tests of the DNA entry's coordinate, merge, design and extension logic (no external tools except blastn)."""
import random
import unittest

import core_review as cr
import dna_discovery as dd


def rnd(n, seed):
    r = random.Random(seed)
    return ''.join(r.choice('ACGT') for _ in range(n))


def occ(sample, contig, lo, hi, strand, contig_seq, W=dd.PARAMS['window']):
    wl, wh = max(0, lo - W), min(len(contig_seq), hi + W)
    return dict(library_id='L0_0', sample=sample, contig=contig, contig_len=str(len(contig_seq)), lo=str(lo),
                hi=str(hi), strand=strand, identity='1.0', win_lo=str(wl), win_hi=str(wh), window=contig_seq[wl:wh])


class Coordinates(unittest.TestCase):
    def test_oriented_minus_strand(self):
        g = rnd(20000, 1)
        o = occ('s', 'c', 9000, 10000, '-', g)
        w, c0, c1 = dd.oriented(o)
        self.assertEqual(w[c0:c1], cr.revcomp(g[9000:10000]))
        o2 = occ('s', 'c', 9000, 10000, '+', g)
        w2, a0, a1 = dd.oriented(o2)
        self.assertEqual(w2[a0:a1], g[9000:10000])

    def test_genomic_extension_both_strands(self):
        o = dict(lo='1000', hi='2000', strand='+')
        self.assertEqual(dd.genomic(o, 30, 70), (970, 2070))
        o = dict(lo='1000', hi='2000', strand='-')  # library upstream is genome downstream
        self.assertEqual(dd.genomic(o, 30, 70), (930, 2030))


class Merge(unittest.TestCase):
    def test_merge_loci(self):
        self.assertEqual(dd.merge_loci([(0, 100), (5, 100), (500, 600)], 0.8), [(0, 100), (500, 600)])
        self.assertEqual(dd.merge_loci([(0, 100), (50, 150)], 0.8), [(0, 100), (50, 150)])

    def test_merge_occurrences_keeps_strands_apart_and_window(self):
        g = rnd(30000, 2)
        a = occ('s', 'c', 10000, 11000, '+', g)
        b = occ('s', 'c', 10500, 11800, '+', g)
        c = occ('s', 'c', 10200, 11000, '-', g)
        m = dd.merge_occurrences([a, b, c])
        self.assertEqual(len(m), 2)
        p = [x for x in m if x['strand'] == '+'][0]
        self.assertEqual((p['lo'], p['hi']), ('10000', '11800'))
        self.assertEqual(p['window'], g[int(p['win_lo']):int(p['win_hi'])])


class Design(unittest.TestCase):
    def test_all_pairs_and_cover_design(self):
        g = rnd(1000, 3)
        os_ = [dict(occ('s%d' % (i % 30), 'c', 0, 10, '+', g)) for i in range(300)]
        for i, o in enumerate(os_):
            o['lo'] = str(i)
        pairs, total, how = dd.design(os_[:100], dd.PARAMS)
        self.assertEqual((len(pairs), total, how), (4950, 4950, 'all_pairs'))
        pairs, total, how = dd.design(os_, dd.PARAMS)
        self.assertEqual(how, 'occurrence_cover_design')
        self.assertEqual(total, 300 * 299 // 2)
        seen = {i for p in pairs for i in p}
        self.assertEqual(seen, set(range(300)))  # every occurrence is in a formed pair
        self.assertEqual(pairs, dd.design(os_, dd.PARAMS)[0])  # deterministic


class Extension(unittest.TestCase):
    """Same X (1,200 bp) in two different contexts: the extension finds X's ends from a seed inside it, on either
    strand, and a pair whose shared block reaches the window end stays unresolved."""

    def setUp(self):
        self.X = rnd(1200, 10)
        self.gA = rnd(8000, 11) + self.X + rnd(8000, 12)
        self.gB = rnd(9000, 13) + cr.revcomp(self.X) + rnd(7000, 14)

    def pair(self, A, B):
        for o in (A, B):
            o['_w'] = dd.oriented(o)
        hits = cr.blast('blastn', {'0': A['_w'][0], '1': B['_w'][0]}, 'blastn').get(('0', '1'), [])
        return dd.extend_pair(A, B, hits, dd.PARAMS)

    def test_resolved_minus_strand(self):
        A = occ('a', 'c', 8300, 8900, '+', self.gA)            # seed inside X, plus strand
        B = occ('b', 'c', 9000 + 1200 - 900, 9000 + 1200 - 300, '-', self.gB)  # same seed, X reversed in B
        e = self.pair(A, B)
        self.assertEqual(e['status'], 'resolved')
        self.assertEqual(dd.genomic(A, e['a_up'], e['a_down']), (8000, 9200))
        self.assertEqual(dd.genomic(B, e['b_up'], e['b_down']), (9000, 10200))

    def test_same_context_unresolved(self):
        A = occ('a', 'c', 8300, 8900, '+', self.gA)
        B = occ('b', 'c', 8300, 8900, '+', self.gA)  # identical locus copy: shared to the window ends
        e = self.pair(A, B)
        self.assertEqual(e['status'], 'unresolved')
        self.assertEqual(e['up'], 'shared_to_window_or_contig_end')


class AlignV02(unittest.TestCase):
    """v0.2: canonical alignment (first window forward, second in relative orientation, one subject per call) and
    the transform into each source's frame give the same extension as a direct search in that frame."""

    def setUp(self):
        self.X = rnd(1200, 20)
        self.gA = rnd(8000, 21) + self.X + rnd(8000, 22)
        self.gB = rnd(9000, 23) + cr.revcomp(self.X) + rnd(7000, 24)

    def occ(self, sample, g, lo, hi, strand):
        o = dict(sample=sample, contig='c', lo=str(lo), hi=str(hi), strand=strand, identity='1.0', clen=str(len(g)),
                 contig_len=str(len(g)))
        o['wlo'], o['whi'] = dd.grid_window(lo, hi, len(g), dd.PARAMS)
        return o

    def direct(self, A, B, gA, gB):
        recs = {}
        for k, O, g in (('0', A, gA), ('1', B, gB)):
            s = g[int(O['wlo']):int(O['whi'])]
            recs[k] = s if O['strand'] == '+' else cr.revcomp(s)
            c0, c1 = dd.seed_in_window(O)
            O['_w'] = (recs[k], c0, c1)
        hs = cr.blast('blastn', recs, 'blastn').get(('0', '1'), [])
        return dd.extend_pair(A, B, hs, dd.PARAMS)

    def canonical(self, A, B, gA, gB):
        key, first, second, rel = dd.align_key(A, B)
        seqs = {first: (gA if first[0] == A['sample'] else gB)[first[2]:first[3]],
                second: (gA if second[0] == A['sample'] else gB)[second[2]:second[3]]}
        subj = seqs[second] if rel == '+' else cr.revcomp(seqs[second])
        hs = dd._blast_group((subj, [(key, seqs[first])]))[key]
        wa = (A['sample'], A['contig'], int(A['wlo']), int(A['whi']))
        a_first = wa == first
        flip = (A['strand'] != '+') if a_first else (A['strand'] != rel)
        hsps = dd.to_source(hs, first[3] - first[2], second[3] - second[2], a_first, flip)
        for O in (A, B):
            c0, c1 = dd.seed_in_window(O)
            O['_w'] = (None, c0, c1)
            O['_len'] = int(O['whi']) - int(O['wlo'])
        return dd.extend_pair(A, B, hsps, dd.PARAMS)

    def test_all_orientations_and_orders(self):
        n = 0
        for sa in '+-':
            sb = '-' if sa == '+' else '+'  # X is reversed in B: library orientation differs between the copies
            A = self.occ('a', self.gA, 8300, 8900, sa)
            B = self.occ('b', self.gB, 9300, 9900, sb)
            for X, Y, gx, gy in ((A, B, self.gA, self.gB), (B, A, self.gB, self.gA)):
                d = self.direct(dict(X), dict(Y), gx, gy)
                c = self.canonical(dict(X), dict(Y), gx, gy)
                keys = ('status', 'a_up', 'a_down', 'b_up', 'b_down')
                self.assertEqual({k: d.get(k) for k in keys}, {k: c.get(k) for k in keys})
                self.assertEqual(c['status'], 'resolved')
                ext = {'a': (X, c['a_up'], c['a_down']), 'b': (Y, c['b_up'], c['b_down'])}
                for s, (O, up, down) in ext.items():
                    x = dd.genomic(O, up, down)
                    self.assertEqual(x, (8000, 9200) if O['sample'] == 'a' else (9000, 10200))
                n += 1
        self.assertEqual(n, 4)

    def test_key_symmetry(self):
        A = self.occ('a', self.gA, 8300, 8900, '+')
        B = self.occ('b', self.gB, 9300, 9900, '-')
        A2, B2 = dict(A, strand='-'), dict(B, strand='+')  # same pair seen in the other global orientation
        self.assertEqual(dd.align_key(A, B)[0], dd.align_key(B, A)[0])
        self.assertEqual(dd.align_key(A, B)[0], dd.align_key(A2, B2)[0])

    def test_pair_chunks_cover_every_pair_once(self):
        g = rnd(1000, 3)
        os_ = [dict(occ('s%d' % (i % 7), 'c', i, i + 10, '+', g)) for i in range(250)]
        for i, o in enumerate(os_):
            o['lo'] = str(i)
        seen = []
        for c, how, it, size in dd.pair_chunks(len(os_), os_, dd.PARAMS):
            got = list(it)
            self.assertEqual(len(got), size)
            seen += got
        self.assertEqual(len(seen), 250 * 249 // 2)
        self.assertEqual(len(set(seen)), len(seen))


if __name__ == '__main__':
    unittest.main()
