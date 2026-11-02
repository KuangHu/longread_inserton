"""Strict P/A evidence, per molecule, for candidate structures of the read-site search (not_in_assembly:core_absent).

A candidate is a merged position whose two flanks are adjacent in the sample's own assembly while reads carry an
insertion between them: not_in_assembly:core_absent positions with >= MAIN_COUNT molecules of the current search
(v0.3: host anchors found beyond extra sequence), plus the v0.1 ones (development regression cases, from the v0.1
backup tables). lp / rp are the junctions (inner ends of the host anchors; rp < lp when the flanks overlap). The
records of a candidate are the current search observations of its core (flank MAIN_FLANK) with both flanks placed at
lp / rp (canonical, within TOL). The inserted sequence is the read sequence between the two junctions: extra left +
core hit + extra right (insert = read[qs - left extra, qe + right extra)), taken from the best-identity record, in
assembly orientation, so P is built from the reads themselves. Each candidate is
tested with the read check's own rules (s3_readcheck v0.3 via ucd_readcheck, thresholds unchanged) against two local
alleles built from the same assembly:
  P = asm[lp - FLANK, lp) + insert + asm[rp, rp + FLANK)   the structure the core-carrying reads show
  A = asm[lp - FLANK, rp + FLANK)                           the assembly as it is (empty site)
Empty-site molecules are also counted with the v0.2 search rule (read_site_search.empty_site_reads / _counts: genome
alignment over the site +/- EMPTY_SPAN, local CIGAR, molecules) for comparison with v0.1's count. Reads: every record with a whole-assembly alignment
(read check genome PAF) within FLANK of the site, plus the core-carrying records; aligned to that candidate's two
alleles only (align_each_interval), counted by rc3 rules (span + ANCHOR, CIGAR limit, allele margin, genome-wide anchor
competition). Records are merged into molecules by FASTQ-header grouping across all runs of the sample.

Per molecule, strict_call is the caller's call or, for a molecule without one, the step where its records leave the
count (trace_intervals steps). Nothing is dropped:
  evidence P / A           P / A (matched when the read is at least as long as the longer tested span)
  evidence conflict        records of the molecule disagree (P vs A vs ambiguous, within or across runs), or the
                           search saw the core in this molecule at this position while the strict call is A
  evidence unknown         ambiguous, cigar_fail, off_locus, locus_ambiguous, partial:*, no_allele_hit
Both forms present is the expected outcome for a structure in part of the DNA; no empty site is required and no
verdict is turned into an insertion event or mobility. Figures: one SVG per candidate, one row per molecule, its best
alignment drawn on the P allele (U | core | D); an A alignment is projected onto P with the core skipped (dashed).

Usage: python candidate_evidence.py [--threads N]   (SLURM; never downloads: ready-manifest runs only)
"""
import argparse
import gzip
import json
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

import protein_context as pc
import read_site_search as rs
import s3_readcheck as rc3
import ucd_readcheck as u

OUT = rs.REPO / 'read_site_search_pilot' / 'candidates'
CLS = 'not_in_assembly:core_absent'
EVIDENCE = {'P_matched': 'P', 'P': 'P', 'A_matched': 'A', 'A': 'A', 'ambiguous_molecule_disagrees': 'conflict',
            'molecule_disagrees_across_runs': 'conflict'}
ORDER = ['P_matched', 'P', 'A_matched', 'A', 'conflict', 'ambiguous', 'cigar_fail', 'locus_ambiguous', 'off_locus',
         'partial:clipped', 'partial:ends', 'partial:short_read', 'no_allele_hit']


def candidates(work, out):
    """v0.2 main-setting core_absent positions with >= MAIN_COUNT molecules, plus the v0.1 ones (from the backup of
    the v0.1 per-sample tables), merged by (sample, core, contig, flank ends within TOL)."""
    got = []

    def add(r, version):
        k = (r['sample'], r['core_key'], r['contig'], int(r['left_pos']), int(r['right_pos']))
        for c in got:
            if c['key'][:3] == k[:3] and abs(c['key'][3] - k[3]) <= rs.TOL and abs(c['key'][4] - k[4]) <= rs.TOL:
                c['found_in'].add(version)
                c.setdefault(version, r)
                return
        got.append(dict(key=k, found_in={version}, **{version: r}))
    for r in pc.read_tsv(out / 'positions.tsv.gz'):
        if int(r['flank']) == rs.MAIN_FLANK and r['cls'] == CLS and int(r['molecules']) >= rs.MAIN_COUNT:
            add(r, 'current')
    for s in rs.pilot_samples(out):
        if not (work / 'v0.1_backup' / s).exists():  # a batch without a v0.1 run (validation)
            continue
        for r in pc.read_tsv(work / 'v0.1_backup' / s / 'positions.tsv.gz'):
            if int(r['flank']) == rs.MAIN_FLANK and r['cls'] == CLS and int(r['molecules']) >= rs.MAIN_COUNT:
                add(dict(r, sample=s), 'v0.1')
    for c in got:
        c['cid'] = '%s:%s:%d-%d:%s' % (c['key'][0], c['key'][2], c['key'][3], c['key'][4], c['key'][1])
    return got


def core_records(work, c):
    """Current search observations (MAIN_FLANK) of this core with both flanks placed whose canonical flank ends are
    those of the candidate within TOL."""
    s, ck, contig, lp, rp = c['key']
    out = []
    for o in pc.read_tsv(work / 'samples' / s / 'observations.tsv.gz'):
        if int(o['flank']) != rs.MAIN_FLANK or o['core_key'] != ck or o['status'] != 'placed':
            continue
        for k in ('left_pos', 'right_pos', 'qs', 'qe', 'left_extra', 'right_extra'):
            o[k] = int(o[k])
        lc, a, ls, rc_, b, rs_ = rs.canonical(o)
        if lc == rc_ == contig and ls == rs_ == '+' and abs(a - lp) <= rs.TOL and abs(b - rp) <= rs.TOL:
            out.append(o)
    return out


def test_of(cid, contig, lp, rp, core, seq):
    """core: the inserted sequence (here core + unaligned flank ends)."""
    """rc3 test row and the two alleles (see module doc)."""
    F = rc3.FLANK
    assert lp >= F and rp + F <= len(seq), (cid, 'site too close to a contig end')
    u_, d_ = seq[lp - F:lp], seq[rp:rp + F]
    p, a = u_ + core + d_, seq[lp - F:rp + F]
    t = dict(test_id=cid, locus_id=cid)
    for role in 'CE':
        t.update({'%s_site' % role: '%s:%d-%d' % (contig, lp, rp),
                  '%sP_j' % role: '%d,%d,0,0' % (F, F + len(core)), '%sP_len' % role: len(p),
                  '%sA_j' % role: '%d,%d,0,0' % (F, F + rp - lp), '%sA_len' % role: len(a)})
    return t, [('%s|C|P' % cid, p), ('%s|C|A' % cid, a)]


def read_steps(apaf, galns, tests, geo):
    """Per (read, test): the caller's read-level call (rc3.read_calls) or the trace step, plus the read's best
    alignment per allele (for the figure). Same best-alignment rule as rc3.classify."""
    best, detail, qlen = defaultdict(dict), defaultdict(dict), {}
    for b in rc3.parse_paf(apaf):
        tid, role, al = b['t'].split('|')
        lo, hi, _ = geo[tid][role]['iv'][al]
        st = (2 if rc3.cigar_ok(b, lo, hi, geo[tid][role]['limit']) else 1) if b['ts'] <= lo and b['te'] >= hi else 0
        qlen[b['q']] = b['qlen']
        cur = best[b['q']].get((tid, role, al))
        if cur is None or (st, b['score']) > cur[:2]:
            best[b['q']][(tid, role, al)] = (st, b['score'], (b['strand'], rc3.q_at(b, lo), rc3.q_at(b, hi - 1))
                                             if st == 2 else None)
            detail[(b['q'], tid)][al] = b
    out = {}
    for q, hits in best.items():
        calls = rc3.read_calls(hits, galns.get(q, []), geo, qlen[q])
        for tid in {k[0] for k in hits}:
            call = calls.get((tid, 'C'))
            d = detail[(q, tid)]
            if not call:
                al = max(d, key=lambda x: d[x]['te'] - d[x]['ts'])
                b, (ilo, ihi, _) = d[al], geo[tid]['C']['iv'][al]
                miss_l, miss_r = max(0, b['ts'] - ilo), max(0, ihi - b['te'])
                left_q = b['qs'] if b['strand'] == '+' else b['qlen'] - b['qe']
                right_q = b['qlen'] - b['qe'] if b['strand'] == '+' else b['qs']
                call = ('partial:clipped' if (miss_l and left_q > 100) or (miss_r and right_q > 100) else
                        'partial:short_read' if qlen[q] < ihi - ilo else 'partial:ends')
            out[(q, tid)] = (call, d)
    return out


def molecule_call(calls):
    """Molecule-level call from its records' calls (rc3.classify rule, extended across runs)."""
    base = {c.split('_matched')[0] for c in calls}
    if len(base & {'P', 'A', 'ambiguous'}) > 1:
        return 'ambiguous_molecule_disagrees'
    for x in ('P', 'A', 'ambiguous', 'locus_ambiguous', 'off_locus', 'cigar_fail', 'partial:clipped', 'partial:ends',
              'partial:short_read'):
        if x in base:
            return x + '_matched' if x in ('P', 'A') and x + '_matched' in calls else x
    return 'no_allele_hit'


def run_sample(s, cands, work, threads):
    od = work / 'candidates' / s
    od.mkdir(parents=True, exist_ok=True)
    runs = u.run_info()[s]
    u.reads_ready.require([r['run'] for r in runs])
    seq = rc3.read_fasta(rs.manifest()[s]['fasta'])
    cores = {c['cid']: core_records(work, c) for c in cands}
    # core sequence: the best-identity core-carrying record (ties: read name), in assembly orientation
    pick = {cid: max(os_, key=lambda o: (float(o['identity']), o['read'])) for cid, os_ in cores.items() if os_}
    want = {o['read'] for os_ in cores.values() for o in os_}
    win = defaultdict(list)
    for c in cands:
        lo, hi = min(c['key'][3], c['key'][4]), max(c['key'][3], c['key'][4])
        win[c['key'][2]].append((lo - rc3.FLANK, hi + rc3.FLANK, c['cid']))
    near, gl = defaultdict(set), defaultdict(list)  # cid -> reads near; read -> genome PAF lines
    for paf in sorted((u.PROJ / 'core_context' / 'readmap' / s).glob('*.paf.gz')):
        with gzip.open(str(paf), 'rt') as h:
            for line in h:
                f = line.split('\t', 9)
                for lo, hi, cid in win.get(f[5], ()):
                    if int(f[7]) < hi and int(f[8]) > lo:
                        near[cid].add(f[0])
    for cid, os_ in cores.items():
        near[cid] |= {o['read'] for o in os_}
    names = set().union(*near.values())
    for paf in sorted((u.PROJ / 'core_context' / 'readmap' / s).glob('*.paf.gz')):
        with gzip.open(str(paf), 'rt') as h:
            for line in h:
                if line.split('\t', 1)[0] in names:
                    gl[line.split('\t', 1)[0]].append(line)
    recs, seq_of, run_of, order = [], {}, {}, defaultdict(list)
    for r in runs:
        with gzip.open(r['local_file'], 'rt') as h:
            keep = False
            for i, line in enumerate(h):
                if i % 4 == 0:
                    name = line[1:].split(None, 1)[0]
                    keep = name in names
                    if keep:
                        recs.append((name, rc3.record_ids(line[1:].strip(), r['grouping'])))
                        run_of[name] = r['run']
                        order[r['run']].append(name)
                elif i % 4 == 1 and keep:
                    seq_of[name] = line
    mol = rc3.molecules(recs)
    tests, fa, info = [], [], {}
    for c in cands:
        cid = c['cid']
        if cid not in pick:
            info[cid] = dict(status='no_placed_record_in_current_search')
            continue
        o = pick[cid]
        core = seq_of[o['read']].strip()[o['qs'] - o['left_extra']:o['qe'] + o['right_extra']]
        core = core if o['left_strand'] == '+' else rc3.rc(core)
        t, alleles = test_of(cid, c['key'][2], c['key'][3], c['key'][4], core, seq[c['key'][2]])
        tests.append(t)
        fa += alleles
        # extra sequence in assembly orientation (a record read from the other strand has its sides swapped)
        cl = [(x['left_extra'], x['right_extra']) if x['left_strand'] == '+' else
              (x['right_extra'], x['left_extra']) for x in cores[cid]]
        info[cid] = dict(status='tested', core_from_read=o['read'], core_len=o['qe'] - o['qs'], insert_len=len(core),
                         core_identity=o['identity'], P_len=len(alleles[0][1]), A_len=len(alleles[1][1]),
                         left_extra_median=rs.median([x[0] for x in cl]), right_extra_median=rs.median([x[1] for x in cl]),
                         left_extra_range='%d-%d' % (min(x[0] for x in cl), max(x[0] for x in cl)),
                         right_extra_range='%d-%d' % (min(x[1] for x in cl), max(x[1] for x in cl)),
                         search_statuses=json.dumps(dict(Counter('%s|%s' % (x['left_status'], x['right_status'])
                                                                 for x in cores[cid]))))
    (od / 'alleles.fa').write_text(''.join('>%s\n%s\n' % x for x in fa))
    geo = {t['test_id']: rc3.test_geometry(t) for t in tests}
    per_rec, audit, galn_all = {}, {}, defaultdict(list)
    for r in runs:
        if not order[r['run']]:
            continue
        apaf, gpaf = od / ('alleles_%s.paf' % r['run']), od / ('genome_%s.paf' % r['run'])
        gpaf.write_text(''.join(l for n in order[r['run']] for l in gl[n]))
        reads_of = {t['test_id']: {n for n in near[t['test_id']] if run_of.get(n) == r['run']} for t in tests}
        u.align_each_interval(tests, fa, reads_of, order[r['run']], seq_of, rc3.PRESET[r['read_type']], od / 'tmp',
                              apaf, threads)
        galns = defaultdict(list)
        for g in rc3.parse_paf(gpaf):
            galns[g['q']].append(g)
            galn_all[g['q']].append(g)
        per_rec.update(read_steps(apaf, galns, tests, geo))
        counts, _ = rc3.classify(tests, [(apaf, gpaf, mol, r['run'])])
        audit[r['run']] = {t['test_id']: dict(counts.get((t['test_id'], 'C'), Counter())) for t in tests}
    for c in cands:  # empty-site molecules by the v0.2 search rule
        cid = c['cid']
        if info[cid]['status'] != 'tested':
            continue
        contig, lo, hi = c['key'][2], min(c['key'][3:5]), max(c['key'][3:5])
        calls = {q: rs.empty_site_reads([g for g in galn_all[q] if g['t'] == contig], lo, hi, info[cid]['insert_len'])
                 for q in near[cid]}
        info[cid].update({'v02rule_' + k: v for k, v in rs.empty_site_counts(
            {o['read'] for o in cores[cid]}, {q: x for q, x in calls.items() if x}, mol).items()})
    return info, near, per_rec, mol, cores, audit, run_of, geo, {n: len(x.strip()) for n, x in seq_of.items()}


def svg(cid, info, rows, geo, w=900, row=5):
    """One row per molecule on the P allele: U | core | D. Solid: best P alignment; dashed with a gap over the core:
    A alignment projected onto P. Colour by evidence."""
    F, L = rc3.FLANK, info['insert_len']
    n = info['P_len']
    sx = lambda x: 150 + max(0, min(n, x)) * (w - 160) / n
    col = {'P': '#2a7ab0', 'A': '#c2452d', 'conflict': '#7b4fa0', 'unknown': '#aaaaaa'}
    lo_p, hi_p = geo['C']['iv']['P'][:2]
    h = 90 + row * len(rows)
    out = ['<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" font-family="sans-serif" font-size="11">' % (w, h),
           '<rect width="100%" height="100%" fill="white"/>',
           '<text x="10" y="14">%s</text>' % cid,
           '<text x="10" y="28">rows = molecules; solid = alignment to P (U+core+D); dashed = alignment to A drawn with '
           'the core skipped; dotted lines = span needed (core +/- %d bp)</text>' % rc3.ANCHOR,
           '<rect x="%.1f" y="36" width="%.1f" height="10" fill="#cccccc"/>' % (sx(0), sx(F) - sx(0)),
           '<rect x="%.1f" y="36" width="%.1f" height="10" fill="#333333"/>' % (sx(F), sx(F + L) - sx(F)),
           '<rect x="%.1f" y="36" width="%.1f" height="10" fill="#cccccc"/>' % (sx(F + L), sx(n) - sx(F + L)),
           '<text x="%.1f" y="60">U (assembly, %d bp)</text>' % (sx(0), F),
           '<text x="%.1f" y="60">insert %d bp from reads (core hit %d bp + extra sequence)</text>' % (
               sx(F), L, info['core_len']),  # extra sequence = insert - core hit
           '<text x="%.1f" y="60">D (assembly)</text>' % (sx(F + L) + 4)]
    for x in (lo_p, hi_p):
        out.append('<line x1="%.1f" x2="%.1f" y1="34" y2="%d" stroke="#888" stroke-dasharray="1,3"/>' % (sx(x), sx(x), h))
    out.append('<text x="10" y="%d" fill="%s">P</text><text x="25" y="%d" fill="%s">A</text>'
               '<text x="40" y="%d" fill="%s">conflict</text><text x="90" y="%d" fill="%s">unknown</text>' % (
                   78, col['P'], 78, col['A'], 78, col['conflict'], 78, col['unknown']))
    prev = None
    for i, r in enumerate(rows):
        y = 88 + i * row
        c = col[r['evidence']]
        if r['strict_call'] != prev:
            out.append('<text x="10" y="%d" font-size="9">%s</text>' % (y + 3, r['strict_call']))
            prev = r['strict_call']
        al, b = r['_aln']
        if b is None:
            continue
        if al == 'P':
            out.append('<line x1="%.1f" x2="%.1f" y1="%d" y2="%d" stroke="%s" stroke-width="3"/>' % (
                sx(b['ts']), sx(b['te']), y, y, c))
        else:  # A coordinates -> P: positions after the A junction shift by core length - A middle
            j1, j2 = F, F + (info['A_len'] - 2 * F)
            m = lambda a: a if a <= j1 else a - (j2 - j1) + L
            a0, a1 = m(b['ts']), m(b['te'])
            if a0 < F < a1:
                out.append('<line x1="%.1f" x2="%.1f" y1="%d" y2="%d" stroke="%s" stroke-width="3" stroke-dasharray="5,2"/>' % (
                    sx(a0), sx(F), y, y, c))
                out.append('<line x1="%.1f" x2="%.1f" y1="%d" y2="%d" stroke="%s" stroke-width="0.6"/>' % (
                    sx(F), sx(F + L), y, y, c))
                out.append('<line x1="%.1f" x2="%.1f" y1="%d" y2="%d" stroke="%s" stroke-width="3" stroke-dasharray="5,2"/>' % (
                    sx(F + L), sx(a1), y, y, c))
            else:
                out.append('<line x1="%.1f" x2="%.1f" y1="%d" y2="%d" stroke="%s" stroke-width="3" stroke-dasharray="5,2"/>' % (
                    sx(a0), sx(a1), y, y, c))
    out.append('</svg>')
    return '\n'.join(out)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--work', type=Path, default=rs.WORK)
    p.add_argument('--out', type=Path, default=rs.REPO / 'read_site_search_pilot')
    p.add_argument('--threads', type=int, default=16)
    a = p.parse_args()
    od = a.out / 'candidates'
    od.mkdir(parents=True, exist_ok=True)
    cands = candidates(a.work, a.out)
    by_s = defaultdict(list)
    for c in cands:
        by_s[c['key'][0]].append(c)
    mrows, summ = [], []
    for s, cs in sorted(by_s.items()):
        info, near, per_rec, mol, cores, audit, run_of, geo, qlen = run_sample(s, cs, a.work, a.threads)
        for c in cs:
            cid = c['cid']
            v1, v2 = c.get('v0.1', {}), c.get('current', {})
            base = dict(candidate=cid, sample=s, core_key=c['key'][1], contig=c['key'][2], left_pos=c['key'][3],
                        right_pos=c['key'][4], found_in='+'.join(sorted(c['found_in'])),
                        v01_molecules=v1.get('molecules', ''), v01_empty_site_molecules=v1.get('empty_site_molecules', ''),
                        molecules=v2.get('molecules', ''), boundary=v2.get('boundary', ''),
                        insert_size=v2.get('insert_size', ''), **{k: v2.get(k, '') for k in rs.EMPTY_COLS})
            if info[cid]['status'] != 'tested':
                summ.append(dict(base, status=info[cid]['status']))
                continue
            core_mols = {mol.get(o['read'], o['read']) for o in cores[cid]}
            by_mol = defaultdict(list)
            for n in near[cid]:
                by_mol[mol.get(n, n)].append(n)
            rows = []
            for m, ns in by_mol.items():
                calls = {n: per_rec.get((n, cid), ('no_allele_hit', {}))[0] for n in ns}
                call = molecule_call(set(calls.values()))
                if call != 'ambiguous_molecule_disagrees' and len({run_of.get(n) for n in ns}) > 1 and len(
                        {c_.split('_matched')[0] for c_ in calls.values()} & {'P', 'A', 'ambiguous'}) > 1:
                    call = 'molecule_disagrees_across_runs'
                ev = EVIDENCE.get(call, 'unknown')
                note = ''
                if m in core_mols and ev == 'A':
                    ev, note = 'conflict', 'search_core_here_vs_strict_A'
                elif m not in core_mols and ev == 'P':
                    note = 'strict_P_not_seen_by_search'
                # figure: the record with the most informative call; its best P (P-type) or A alignment
                n0 = min(ns, key=lambda n: (ORDER.index(calls[n]) if calls[n] in ORDER else 99, n))
                d = per_rec.get((n0, cid), ('', {}))[1]
                al = 'A' if call.startswith('A') or 'P' not in d else 'P'
                rows.append(dict(base, molecule=m, records=len(ns), runs=';'.join(sorted({run_of.get(n, '') for n in ns})),
                                 read_lengths=';'.join(str(qlen.get(n, '')) for n in sorted(ns)),
                                 record_calls=';'.join('%s:%s' % (n, calls[n]) for n in sorted(ns)),
                                 strict_call=call, evidence=ev, search_core_here='yes' if m in core_mols else 'no',
                                 note=note, _aln=(al, d.get(al))))
            rows.sort(key=lambda r: (ORDER.index(r['strict_call']) if r['strict_call'] in ORDER else
                                     ORDER.index('conflict') if r['evidence'] == 'conflict' else 99,
                                     r['_aln'][1]['ts'] if r['_aln'][1] else 0, r['molecule']))
            (od / ('%s.svg' % cid.replace(':', '_'))).write_text(svg(cid, info[cid], rows, geo[cid]))
            mrows += rows
            n = Counter(r['strict_call'] for r in rows)
            e = Counter(r['evidence'] for r in rows)
            summ.append(dict(base, status='tested', **{k: v for k, v in info[cid].items() if k != 'status'},
                             molecules_near=len(rows), search_core_molecules=len(core_mols),
                             P_matched=n['P_matched'], P=n['P'], A_matched=n['A_matched'], A=n['A'],
                             conflict=e['conflict'], unknown=e['unknown'],
                             strict_P_not_seen_by_search=sum(r['note'] == 'strict_P_not_seen_by_search' for r in rows),
                             search_core_vs_strict_A=sum(r['note'] == 'search_core_here_vs_strict_A' for r in rows),
                             unknown_by_step=json.dumps(dict(Counter(r['strict_call'] for r in rows
                                                                     if r['evidence'] == 'unknown'))),
                             assembly_allele_verdict=rc3.verdict(n['A_matched'], n['P_matched']),
                             per_run_caller_counts=json.dumps({run: x[cid] for run, x in audit.items()})))
    mcols = ['candidate', 'molecule', 'records', 'runs', 'read_lengths', 'strict_call', 'evidence', 'search_core_here',
             'note', 'record_calls']
    pc.write_tsv(od / 'evidence_molecules.tsv', mrows, mcols)
    scols = ['candidate', 'sample', 'core_key', 'contig', 'left_pos', 'right_pos', 'found_in', 'status', 'v01_molecules',
             'v01_empty_site_molecules', 'molecules', 'boundary', 'insert_size'] + rs.EMPTY_COLS + [
             'core_from_read', 'core_len', 'insert_len', 'left_extra_median', 'left_extra_range', 'right_extra_median',
             'right_extra_range', 'search_statuses', 'core_identity', 'P_len', 'A_len'] + [
             'v02rule_' + k for k in rs.EMPTY_COLS] + ['molecules_near', 'search_core_molecules',
             'P_matched', 'P', 'A_matched', 'A', 'conflict', 'unknown', 'unknown_by_step', 'strict_P_not_seen_by_search',
             'search_core_vs_strict_A', 'assembly_allele_verdict', 'per_run_caller_counts']
    pc.write_tsv(od / 'evidence_summary.tsv', summ, scols)
    for r in summ:
        print({k: r.get(k, '') for k in scols if k not in ('per_run_caller_counts',)})


if __name__ == '__main__':
    main()
