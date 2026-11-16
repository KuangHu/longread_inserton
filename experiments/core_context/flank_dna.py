"""Auxiliary flank DNA re-check v2 for S6 review pairs. Reported only: it never changes an S6 grouping, a host-context
class or an admission. Evidence stays assembly_only; production_ready=false.

Fixes the v0.1 check in s6_review.flank_check, which aligned the whole 3 kb flanks and measured coverage only on the
query's outer part: an outer window could count as shared while it matched only the other copy's excluded inner
window (the 500 bp next to the core end, where element ends may extend past the CDS).

v2: both copies' outer windows (INNER..W bp beyond each core end, in core orientation) are cut first; only the outer
windows are aligned (blastn, both strands; its own blast call, the shared core_review.blast used by the boundary
analysis is not touched). For each of the four window combinations (x up / x down vs y up / y down) and each
strand, coverage = min(fraction of x's window covered, fraction of y's window covered) by the HSPs between those two
windows. Reported per pair:
  same-side forward   x up vs y up, x down vs y down, plus strand: the repeated-context signal
  same-side reverse   same windows, minus strand
  cross-side          x up vs y down and x down vs y up, each strand
Class per side from same-side forward coverage: shared >= 0.5, not_shared <= 0.1, partial otherwise, short_flank
if either outer window is < MIN_OUTER bp. Copies on one contig closer than the flank length can have physically
overlapping windows (one copy's downstream window is the other's upstream region): copy_gap_bp and
host_genes_shared_physically (the same gene in both S6 host windows) flag this, since a match or a 'same' host side
there is the same DNA, not a duplicate. 'not_shared' means no outer-flank alignment was detected by this check
(blastn, e <= 1e-5, these windows only); it is not proof that the copies share no DNA.
"""
import argparse
import json
import subprocess
import tempfile
from collections import Counter
from pathlib import Path

import core_review as cr
import protein_context as pc
import s6_review as rv

VERSION = 'flank_dna_v2'
W, INNER, MIN_OUTER = rv.W, rv.INNER, rv.MIN_OUTER
SIDES = ('up', 'down')


def outer_windows(f):
    """{'up': seq, 'down': seq} in core orientation -> the outer windows (up flanks end at the core)."""
    return dict(up=f['up'][:max(0, len(f['up']) - INNER)], down=f['down'][INNER:])


def align(qs, ss):
    """blastn of named query windows against named subject windows, both strands. Returns
    {(query, subject, strand): [(q0, q1, s0, s1)]} with 0-based half-open intervals on each sequence."""
    qs, ss = {k: v for k, v in qs.items() if v}, {k: v for k, v in ss.items() if v}
    out = {}
    if not qs or not ss:
        return out
    with tempfile.NamedTemporaryFile('w', suffix='.fa') as fq, tempfile.NamedTemporaryFile('w', suffix='.fa') as fs:
        for f, recs in ((fq, qs), (fs, ss)):
            for k, v in recs.items():
                f.write('>%s\n%s\n' % (k, v))
            f.flush()
        res = subprocess.run(['blastn', '-query', fq.name, '-subject', fs.name, '-task', 'blastn', '-strand', 'both',
                              '-dust', 'no', '-evalue', '1e-5', '-outfmt', '6 qseqid sseqid qstart qend sstart send'],
                             check=True, capture_output=True, text=True).stdout
    for line in res.splitlines():
        q, s, a, b, c, d = line.split('\t')
        a, b, c, d = int(a), int(b), int(c), int(d)
        strand = 'plus' if c <= d else 'minus'
        out.setdefault((q, s, strand), []).append((a - 1, b, min(c, d) - 1, max(c, d)))
    return out


def outer_flank_check(fx, fy):
    """fx, fy: full flanks {'up', 'down'} of two copies in core orientation. Returns a flat dict of coverages,
    lengths and per-side classes (see module doc)."""
    ox, oy = outer_windows(fx), outer_windows(fy)
    hits = align({'x' + s: ox[s] for s in SIDES}, {'y' + s: oy[s] for s in SIDES})
    out = {}
    for sx in SIDES:
        for sy in SIDES:
            for strand in ('plus', 'minus'):
                hs = hits.get(('x' + sx, 'y' + sy, strand), [])
                lx, ly = len(ox[sx]), len(oy[sy])
                cov = 0.0
                if lx and ly:
                    cov = min(rv.covered([(h[0], h[1]) for h in hs], 0, lx) / lx,
                              rv.covered([(h[2], h[3]) for h in hs], 0, ly) / ly)
                out['cov_x%s_y%s_%s' % (sx, sy, strand)] = round(cov, 3)
    for s in SIDES:
        length = min(len(ox[s]), len(oy[s]))
        out['outer_len_' + s] = length
        out['dna_' + s] = rv.flank_class(out['cov_x%s_y%s_plus' % (s, s)], length)
    out['same_side_reverse_max'] = max(out['cov_xup_yup_minus'], out['cov_xdown_ydown_minus'])
    out['cross_side_max'] = max(v for k, v in out.items() if k.startswith('cov_') and
                                ('xup_ydown' in k or 'xdown_yup' in k))
    out['other_orientation_match'] = max(out['same_side_reverse_max'], out['cross_side_max']) >= 0.5
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--subset-dir', type=Path, required=True)
    p.add_argument('--work', type=Path, required=True, help='protein work dir (genes.tsv.gz)')
    a = p.parse_args()
    rdir = a.subset_dir / 'protein_entry' / 's6' / 'review'
    occs = {o['occurrence_id']: o for o in pc.read_tsv(a.subset_dir / 'protein_entry' / 's6' / 's6_occurrences.tsv.gz')}
    circular = {(g['sample'], g['contig']): g['circular_contig'] == 'True' for g in pc.read_tsv(a.work / 'genes.tsv.gz')}
    fasta = {r['sample']: r['fasta'] for r in pc.read_tsv(a.subset_dir.parent / (a.subset_dir.name + '.samples.tsv'))}
    for smp in sorted({o['sample'] for o in occs.values()}):
        seqs = pc.read_fna(fasta[smp])
        for o in occs.values():
            if o['sample'] == smp:
                sub, c0, c1 = cr.window(seqs[o['contig']], int(o['start0']), int(o['end0']), o['strand'] == '+', W,
                                        circular[(smp, o['contig'])])
                o['flank'] = dict(up=sub[:c0], down=sub[c1:])

    # every same-assembly pair and every case pair of the committed S6 review, with the v0.1 classes next to v2
    todo = [('same_assembly', r) for r in pc.read_tsv(rdir / 's6_same_assembly_pairs.tsv.gz')]
    todo += [('case:' + r['case'], r) for r in pc.read_tsv(rdir / 'review_cases.tsv')]
    rows, done = [], {}
    for src, r in todo:
        x, y = r['occurrence_a'], r['occurrence_b']
        if (x, y) not in done:
            done[(x, y)] = outer_flank_check(occs[x]['flank'], occs[y]['flank'])
        v = done[(x, y)]
        old = (r['dna_up'], r['dna_down'])
        ox, oy = occs[x], occs[y]
        gap = (max(int(ox['start0']), int(oy['start0'])) - min(int(ox['end0']), int(oy['end0']))
               if (ox['sample'], ox['contig']) == (oy['sample'], oy['contig']) else '')
        hg = lambda o: set(filter(None, (o['host_up_genes'] + ',' + o['host_down_genes']).split(',')))
        rows.append(dict(source=src, copy_gap_bp=gap, windows_can_overlap=gap != '' and gap < 2 * W,
                         host_genes_shared_physically=len(hg(ox) & hg(oy)), structure=r['structure'], core_key=r['core_key'], occurrence_a=x, occurrence_b=y,
                         sample_a=occs[x]['sample'], sample_b=occs[y]['sample'],
                         context_class=r['context_class'], v01_dna_up=old[0], v01_dna_down=old[1],
                         v01_cov_up=r['dna_outer_cov_up'], v01_cov_down=r['dna_outer_cov_down'], **v,
                         changed='%s/%s' % old != '%s/%s' % (v['dna_up'], v['dna_down'])))
    cols = list(rows[0])
    pc.write_tsv(rdir / 'flank_dna_v2.tsv.gz', rows, cols)

    summary = dict(version=VERSION, subset=a.subset_dir.name, evidence_source='assembly_only', production_ready=False,
                   use='auxiliary report only; never changes S6 groups, context classes or admission',
                   parameters=dict(flank_bp=W, inner_bp_excluded=INNER, min_outer_bp=MIN_OUTER, blastn='-task blastn '
                                   '-strand both -dust no -evalue 1e-5', shared='>= 0.5', not_shared='<= 0.1 (no outer '
                                   'alignment detected by this check; not proof of no shared DNA)'), structures={})
    for st in ('single_gene', 'two_gene'):
        sa = [r for r in rows if r['source'] == 'same_assembly' and r['structure'] == st]
        summary['structures'][st] = dict(
            same_assembly_pairs=len(sa), changed=sum(r['changed'] for r in sa),
            transitions=dict(Counter('%s/%s -> %s/%s' % (r['v01_dna_up'], r['v01_dna_down'], r['dna_up'], r['dna_down'])
                                     for r in sa if r['changed'])),
            v2_by_context={k: dict(Counter('%s/%s' % (r['dna_up'], r['dna_down']) for r in sa if r['context_class'] == k))
                           for k in sorted({r['context_class'] for r in sa})},
            other_orientation_match=dict(Counter(r['context_class'] for r in sa if r['other_orientation_match'])),
            other_orientation_match_windows_can_overlap=sum(r['other_orientation_match'] and r['windows_can_overlap'] for r in sa),
            windows_can_overlap=dict(Counter(r['context_class'] for r in sa if r['windows_can_overlap'])),
            host_genes_shared_physically=dict(Counter(r['context_class'] for r in sa if r['host_genes_shared_physically'])))
    cs = [r for r in rows if r['source'].startswith('case:')]
    summary['case_pairs'] = len(cs)
    summary['case_pairs_changed'] = [dict(case=r['source'][5:], a=r['occurrence_a'], b=r['occurrence_b'],
                                          v01='%s/%s' % (r['v01_dna_up'], r['v01_dna_down']),
                                          v2='%s/%s' % (r['dna_up'], r['dna_down'])) for r in cs if r['changed']]
    (rdir / 'flank_dna_v2_summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
