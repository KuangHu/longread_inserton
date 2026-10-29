"""Long-read check of U-core-D positions: does the sample's own long-read set support the structure its assembly shows
at each position of a candidate pair? Evidence added: long_read_position_check; assembly_only stays the source of
the pairs; production_ready=false. Passing both positions does not show a mobility mechanism or independent
insertion events.

pairs     Formal validation list (from scale_candidates): every match_both_lt2500 pair with at least one side in an
          assembly with verified local reads, one row per pair, fixed: both occurrences, coordinates, element interval
          per side (coding span plus that pair's shared extensions), core match class / identity / coverage, both
          sides' boundary status in the pair's scope, runs. Tier A = same assembly (one read set covers both
          positions), B = cross assembly, both sides with reads, C = cross assembly, one side with reads (single-side
          support at most). control_class labels conserved paralogs (product hint translation_or_ribosomal); they stay
          in the list. Repeated-context same-assembly pairs (context_repeated_both_sides, not U-core-D) of the
          same assemblies go to a separate control list, tested on the coding span.
select    Round selection: given assemblies, the list's pairs whose read side is among them, plus their controls.
genome    One assembly: SHA256 of each run rechecked against verify.json, then all reads of every run mapped once to
          the whole assembly (s3 MM_GENOME: secondaries kept); PAF kept for later rounds.
check     One assembly, one round: reads with a genome alignment near a tested position are taken from the full
          genome PAF (not from any earlier gap-pipeline subset) and mapped to two local alleles per position built from
          the same assembly: P = U + element + D (as assembled), A = U + D (element removed; no empty site needed).
          Counting is s3_readcheck v0.3 unchanged (span + ANCHOR, CIGAR limit, allele score margin, genome-wide anchor
          competition, molecule grouping, length-matched counts, role verdict).
summarize Per position and per pair evidence; representative SVG figures.

Position classes: supported (consistent), conflict (contradicted or mixed), locus_ambiguous (insufficient, with
molecules whose anchors also place elsewhere), insufficient_coverage (insufficient otherwise), not_tested (no reads or
element interval at a contig end).
"""
import argparse
import gzip
import hashlib
import json
import subprocess
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

import protein_context as pc
import scale_entry as se  # noqa: F401  (puts fna_pool on sys.path)
import s3_readcheck as rc3
from scale_candidates import POOL, READS

VERSION = 'ucd_readcheck_v0.1'
PROJ = Path('/global/scratch/users/kh36969/longread_callout')
CONTROL_HINTS = {'translation_or_ribosomal': 'conserved_paralog_control'}
TIERS = {'A': 'same_assembly', 'B': 'cross_assembly_both_sides_reads', 'C': 'cross_assembly_one_side_reads'}
PAIR_COLS = ['pair_id', 'tier', 'structure', 'core_key', 'product', 'product_hint', 'control_class', 'pair_scope',
             'copy_class', 's6_grouping', 'core_match_class', 'core_dna_identity', 'core_dna_coverage_min',
             'protein_identity_min']
SIDE_COLS = ['occurrence', 'sample', 'contig', 'start0', 'end0', 'strand', 's6_status', 'ext_up', 'ext_down', 'elem_lo',
             'elem_hi', 'elem_len', 'boundary_up', 'boundary_down', 'estimated_span', 'runs']


def run_info(runs_tsv=POOL / 'runs.tsv', reads=READS):
    """{sample: [run rows]} usable runs with a verified local FASTQ whose MD5 equals the archive (scale_candidates)."""
    out = defaultdict(list)
    for r in pc.read_tsv(runs_tsv):
        v = reads / r['run_accession'] / 'verify.json'
        if r['usable'] != 'True' or not v.exists():
            continue
        j = json.loads(v.read_text())
        if j.get('status') == 'verified' and j.get('md5') == r['fastq_md5']:
            out[r['sample_accession']].append(dict(run=r['run_accession'], read_type=r['read_type'],
                                                   model=r['instrument_model'], local_file=j['local_file'],
                                                   bytes=j['bytes'], sha256=j['sha256'],
                                                   grouping=j['summary']['molecule_grouping']))
    return out


def element(o, up, down):
    """Genomic [lo, hi) of coding span plus the shared extensions (given in core orientation)."""
    s, e = int(o['start0']), int(o['end0'])
    return (s - up, e + down) if o['strand'] == '+' else (s - down, e + up)


def tier_of(scope, ra, rb):
    if scope == 'same_assembly':
        return 'A' if ra else None
    return 'B' if ra and rb else 'C' if ra or rb else None


def cmd_pairs(a):
    tag = 'b' + a.batches.replace(',', '-')
    out = a.repo_out
    runs = run_info()
    qual = {}
    for r in pc.read_tsv(out / ('match_both_lt2500_pairs_%s.tsv.gz' % tag)):
        qual[(r['structure'], r['core_key'], r['occurrence_a'], r['occurrence_b'])] = r
    cores = {(r['structure'], r['core_key']): r for r in pc.read_tsv(out / ('core_candidates_%s.tsv.gz' % tag))}
    bnd = {}
    for scope in ('same_assembly', 'cross_assembly'):
        for r in pc.read_tsv(out / ('boundary_occurrences_%s_%s.tsv.gz' % (scope, tag))):
            bnd[(scope, r['structure'], r['occurrence_id'])] = r
    samples_with_reads = set(runs)
    pair_rows, occ_need, controls = {}, set(), []
    for k in range(64):
        for r in pc.read_tsv(a.work / 'review' / 'boundary' / tag / ('shard_%03d' % k) / 'pairs.tsv.gz'):
            key = (r['structure'], r['core_key'], r['occurrence_a'], r['occurrence_b'])
            if key in qual and tier_of(r['pair_scope'], r['sample_a'] in samples_with_reads,
                                       r['sample_b'] in samples_with_reads):
                pair_rows[key] = r
                occ_need.update(key[2:])
        for r in pc.read_tsv(a.work / 'review' / 's6' / tag / ('shard_%03d' % k) / 'same_assembly_pairs.tsv.gz'):
            if r['context_class'] == 'context_repeated_both_sides' and r['sample'] in samples_with_reads:
                controls.append(r)
                occ_need.update((r['occurrence_a'], r['occurrence_b']))
    occ = {}
    for k in range(64):
        for r in pc.read_tsv(a.work / 'review' / 's6' / tag / ('shard_%03d' % k) / 'occurrences.tsv.gz'):
            if r['occurrence_id'] in occ_need:
                occ[r['occurrence_id']] = r
    rows = []
    for key, r in pair_rows.items():
        c = cores[key[:2]]
        tier = tier_of(r['pair_scope'], r['sample_a'] in samples_with_reads, r['sample_b'] in samples_with_reads)
        row = dict(pair_id='VP_' + hashlib.sha256('|'.join(key).encode()).hexdigest()[:12], tier=tier,
                   structure=key[0], core_key=key[1], product=c['product'], product_hint=c['product_hint'],
                   control_class=CONTROL_HINTS.get(c['product_hint'], ''),
                   **{f: r[f] for f in PAIR_COLS[7:] if f in r})
        for s in 'ab':
            o = occ[r['occurrence_' + s]]
            up, down = int(r[s + '_up']), int(r[s + '_down'])
            lo, hi = element(o, up, down)
            b = bnd[(r['pair_scope'], key[0], o['occurrence_id'])]
            row.update({'%s_%s' % (s, f): v for f, v in dict(
                occurrence=o['occurrence_id'], sample=o['sample'], contig=o['contig'], start0=o['start0'], end0=o['end0'],
                strand=o['strand'], s6_status=o['s6_status'] + (':' + o['unresolved_reason'] if o['unresolved_reason'] else ''),
                ext_up=up, ext_down=down, elem_lo=lo, elem_hi=hi, elem_len=hi - lo, boundary_up=b['up_status'],
                boundary_down=b['down_status'], estimated_span=b['estimated_span'],
                runs=';'.join('%s:%s' % (x['run'], x['read_type']) for x in runs.get(o['sample'], []))).items()})
        rows.append(row)
    rank = {'identical_core_dna': 0, 'homolog_full_coverage': 1}
    rows.sort(key=lambda x: (x['tier'], x['structure'], rank[x['core_match_class']], -float(x['core_dna_identity']),
                             -float(x['core_dna_coverage_min']), x['pair_id']))
    cols = PAIR_COLS + ['%s_%s' % (s, f) for s in 'ab' for f in SIDE_COLS]
    pc.write_tsv(out / ('validation_pairs_%s.tsv.gz' % tag), rows, cols)
    crow = []
    for r in controls:
        x = dict(pair_id='VC_' + hashlib.sha256('|'.join((r['structure'], r['core_key'], r['occurrence_a'],
                                                          r['occurrence_b'])).encode()).hexdigest()[:12],
                 tier='control', structure=r['structure'], core_key=r['core_key'],
                 control_class='repeated_context_control', pair_scope='same_assembly', s6_grouping=r['s6_grouping'])
        for s in 'ab':
            o = occ[r['occurrence_' + s]]
            x.update({'%s_%s' % (s, f): v for f, v in dict(
                occurrence=o['occurrence_id'], sample=o['sample'], contig=o['contig'], start0=o['start0'],
                end0=o['end0'], strand=o['strand'], s6_status=o['s6_status'], ext_up=0, ext_down=0,
                elem_lo=o['start0'], elem_hi=o['end0'], elem_len=int(o['end0']) - int(o['start0']),
                runs=';'.join('%s:%s' % (y['run'], y['read_type']) for y in runs.get(o['sample'], []))).items()})
        crow.append(x)
    crow.sort(key=lambda x: x['pair_id'])
    pc.write_tsv(out / ('validation_controls_%s.tsv.gz' % tag), crow, cols)
    summ = dict(version=VERSION, evidence_source='assembly_only', production_ready=False, batches=a.batches,
                assemblies_with_verified_local_reads=len(runs),
                pairs={st: dict(Counter(r['tier'] for r in rows if r['structure'] == st)) for st in ('single_gene', 'two_gene')},
                pairs_conserved_paralog_control=dict(Counter(r['tier'] for r in rows if r['control_class'])),
                pair_assemblies={t: len({r[s + '_sample'] for r in rows if r['tier'] == t for s in 'ab'
                                         if r[s + '_sample'] in runs}) for t in 'ABC'},
                repeated_context_controls=dict(Counter(r['structure'] for r in crow)))
    (out / ('validation_pairs_%s.json' % tag)).write_text(json.dumps(summ, indent=1) + '\n')
    print(json.dumps(summ, indent=1))


def sides(row, runs):
    """(side, sample, contig, lo, hi) of the positions of a pair that are in assemblies with reads."""
    return [(s, row[s + '_sample'], row[s + '_contig'], int(row[s + '_elem_lo']), int(row[s + '_elem_hi']))
            for s in 'ab' if row[s + '_sample'] in runs]


def test_id(sample, contig, lo, hi):
    return 'P_' + hashlib.sha256(('%s|%s|%d|%d' % (sample, contig, lo, hi)).encode()).hexdigest()[:12]


def cmd_select(a):
    tag = 'b' + a.batches.replace(',', '-')
    want = set(a.assemblies.split(','))
    rd = a.round_dir
    rd.mkdir(parents=True, exist_ok=True)
    rows = pc.read_tsv(a.repo_out / ('validation_pairs_%s.tsv.gz' % tag))
    ctrl = pc.read_tsv(a.repo_out / ('validation_controls_%s.tsv.gz' % tag))
    pick = [r for r in rows if any(r[s + '_sample'] in want for s in 'ab')]
    pick += [r for r in ctrl if r['a_sample'] in want]
    cols = list(rows[0])
    pc.write_tsv(rd / 'pairs.tsv', pick, cols)
    pos = {}
    for r in pick:
        for s in 'ab':
            if r[s + '_sample'] in want:
                lo, hi = int(r[s + '_elem_lo']), int(r[s + '_elem_hi'])
                t = test_id(r[s + '_sample'], r[s + '_contig'], lo, hi)
                p = pos.setdefault(t, dict(test_id=t, sample=r[s + '_sample'], contig=r[s + '_contig'], lo=lo, hi=hi,
                                           occurrence=r[s + '_occurrence'], pairs=set()))
                p['pairs'].add(r['pair_id'])
    for p in pos.values():
        p['pairs'] = ';'.join(sorted(p['pairs']))
    pc.write_tsv(rd / 'positions.tsv', sorted(pos.values(), key=lambda x: (x['sample'], x['contig'], x['lo'])),
                 ['test_id', 'sample', 'contig', 'lo', 'hi', 'occurrence', 'pairs'])
    (rd / 'samples.txt').write_text('\n'.join(sorted(want)) + '\n')
    summ = dict(assemblies=len(want), pairs=dict(Counter('%s:%s' % (r['tier'], r['structure']) for r in pick)),
                positions=len(pos), pairs_with_control_class=dict(Counter(r['control_class'] for r in pick if r['control_class'])))
    (rd / 'selection.json').write_text(json.dumps(summ, indent=1) + '\n')
    print(json.dumps(summ, indent=1))


def fasta_of(repo_out, sample):
    return next(r['fasta'] for r in pc.read_tsv(repo_out / 'manifest.tsv') if r['sample'] == sample)


def cmd_genome(a):
    """Map every read of every verified run of one assembly to the whole assembly, once."""
    runs = run_info()[a.sample]
    fa = fasta_of(a.repo_out, a.sample)
    od = PROJ / 'core_context' / 'readmap' / a.sample
    od.mkdir(parents=True, exist_ok=True)
    done = []
    for r in runs:
        paf = od / ('%s.paf.gz' % r['run'])
        meta = od / ('%s.json' % r['run'])
        if meta.exists() and paf.exists():
            done.append(json.loads(meta.read_text()))
            continue
        t0 = time.time()
        sha = pc.digest(r['local_file'])
        if sha != r['sha256']:
            done.append(dict(run=r['run'], status='sha256_mismatch'))
            continue
        cmd = 'minimap2 -x %s %s -t %d %s %s | gzip -1 > %s.tmp' % (
            rc3.PRESET[r['read_type']], ' '.join(rc3.MM_GENOME), a.threads, fa, r['local_file'], paf)
        subprocess.run(['bash', '-o', 'pipefail', '-c', cmd], check=True, stderr=subprocess.DEVNULL)
        Path(str(paf) + '.tmp').rename(paf)
        m = dict(run=r['run'], status='mapped', read_type=r['read_type'], preset=rc3.PRESET[r['read_type']],
                 minimap2_args=rc3.MM_GENOME, assembly=fa, sha256=sha, seconds=round(time.time() - t0))
        meta.write_text(json.dumps(m, indent=1) + '\n')
        done.append(m)
    (od / 'genome.json').write_text(json.dumps(dict(sample=a.sample, runs=done), indent=1) + '\n')
    print(json.dumps(done, indent=1))


def make_tests(positions, seqs):
    """s3 test rows (role C = this assembly; E mirrors C and has no alleles) and the allele FASTA records."""
    tests, fa, skipped = [], [], []
    for p in positions:
        seq = seqs[p['contig']]
        lo, hi = int(p['lo']), int(p['hi'])
        if lo < 0 or hi > len(seq):
            skipped.append(dict(p, reason='element_interval_crosses_contig_end'))
            continue
        u, el, d = seq[max(0, lo - rc3.FLANK):lo], seq[lo:hi], seq[hi:hi + rc3.FLANK]
        t = dict(test_id=p['test_id'], locus_id=p['test_id'])
        for role in 'CE':
            t.update({'%s_site' % role: '%s:%d-%d' % (p['contig'], lo, hi),
                      '%sP_j' % role: '%d,%d,0,0' % (len(u), len(u) + len(el)), '%sP_len' % role: len(u + el + d),
                      '%sA_j' % role: '%d,%d,0,0' % (len(u), len(u)), '%sA_len' % role: len(u + d)})
        tests.append(t)
        fa.append(('%s|C|P' % p['test_id'], u + el + d))
        fa.append(('%s|C|A' % p['test_id'], u + d))
    return tests, fa, skipped


def cmd_check(a):
    rd = a.round_dir
    positions = [p for p in pc.read_tsv(rd / 'positions.tsv') if p['sample'] == a.sample]
    seqs = rc3.read_fasta(fasta_of(a.repo_out, a.sample))
    tests, fa, skipped = make_tests(positions, seqs)
    od = rd / 'check' / a.sample
    od.mkdir(parents=True, exist_ok=True)
    afa = od / 'alleles.fa'
    afa.write_text(''.join('>%s\n%s\n' % x for x in fa))
    near = defaultdict(list)
    for p in positions:
        near[p['contig']].append((int(p['lo']) - rc3.FLANK, int(p['hi']) + rc3.FLANK))
    gmeta = json.loads((PROJ / 'core_context' / 'readmap' / a.sample / 'genome.json').read_text())
    runs = {r['run']: r for r in run_info()[a.sample]}
    items, used, plot = [], [], defaultdict(list)
    for m in gmeta['runs']:
        if m['status'] != 'mapped':
            used.append(dict(run=m['run'], status=m['status']))
            continue
        r = runs[m['run']]
        gpaf = PROJ / 'core_context' / 'readmap' / a.sample / ('%s.paf.gz' % r['run'])
        names, lines, total = set(), [], set()
        with gzip.open(str(gpaf), 'rt') as h:
            for line in h:
                f = line.split('\t', 9)
                total.add(f[0])
                if any(int(f[7]) < hi and int(f[8]) > lo for lo, hi in near.get(f[5], ())):
                    names.add(f[0])
        sub_g = od / ('genome_%s.paf' % r['run'])
        with gzip.open(str(gpaf), 'rt') as h, open(str(sub_g), 'w') as o:
            for line in h:
                if line.split('\t', 1)[0] in names:
                    o.write(line)
        sub_fq = od / ('subset_%s.fq' % r['run'])
        recs, keep = [], False
        with gzip.open(r['local_file'], 'rt') as h, open(str(sub_fq), 'w') as o:
            for i, line in enumerate(h):
                if i % 4 == 0:
                    name = line[1:].split(None, 1)[0]
                    keep = name in names
                    if keep:
                        recs.append((name, rc3.record_ids(line[1:].strip(), r['grouping'])))
                if keep:
                    o.write(line)
        mol = rc3.molecules(recs)
        apaf = od / ('alleles_%s.paf' % r['run'])
        with open(str(apaf), 'w') as h:
            subprocess.run(['minimap2', '-x', rc3.PRESET[r['read_type']]] + rc3.MM + ['-t', str(a.threads), str(afa),
                                                                                       str(sub_fq)],
                           stdout=h, stderr=subprocess.DEVNULL, check=True)
        sub_fq.unlink()
        items.append((apaf, sub_g, mol, r['run']))
        used.append(dict(run=r['run'], status='checked', read_type=r['read_type'], grouping=r['grouping'],
                         reads_mapped_to_assembly=len(total), reads_near_positions=len(names),
                         molecules_near_positions=len(set(mol.values())),
                         count_unit='read_records' if r['grouping'] == 'molecule_grouping_unavailable' else 'molecules'))
    records = []
    counts, multi = rc3.classify(tests, items, records) if items else ({}, Counter())
    call_of = defaultdict(dict)
    for tid, role, run, m, call, reads in records:
        for x in reads.split(';'):
            call_of[(run, x.rsplit(':', 1)[0])][tid] = call
    # figure data: genome alignments of the reads that got a call at a position, on that position's contig
    pos_of = {p['test_id']: p for p in positions}
    for apaf, sub_g, mol, run in items:
        for g in rc3.parse_paf(sub_g):
            for tid, call in call_of.get((run, g['q']), {}).items():
                p = pos_of[tid]
                if g['t'] == p['contig'] and g['te'] > int(p['lo']) - rc3.FLANK and g['ts'] < int(p['hi']) + rc3.FLANK:
                    plot[tid].append((call, g['ts'], g['te'], g['strand'], g['tp'], g['q']))
    rows = []
    for t in tests:
        c = counts.get((t['test_id'], 'C'), Counter())
        v = rc3.verdict(c['P_matched'], c['A_matched']) if items else 'no_reads'
        rows.append(dict(test_id=t['test_id'], sample=a.sample, verdict=v, **{k: c[k] for k in (
            'P', 'A', 'P_matched', 'A_matched', 'ambiguous', 'locus_ambiguous', 'off_locus', 'cigar_fail')}))
    for p in skipped:
        rows.append(dict(test_id=p['test_id'], sample=a.sample, verdict='not_tested:' + p['reason']))
    pc.write_tsv(od / 'positions_result.tsv', rows, ['test_id', 'sample', 'verdict', 'P', 'A', 'P_matched', 'A_matched',
                                                     'ambiguous', 'locus_ambiguous', 'off_locus', 'cigar_fail'])
    (od / 'plot.json').write_text(json.dumps(plot) + '\n')
    (od / 'run.json').write_text(json.dumps(dict(version=VERSION, caller=rc3.VERSION, sample=a.sample, runs=used,
                                                 molecules_supporting_loci=dict(multi)), indent=1) + '\n')
    for x in od.glob('*.paf'):
        subprocess.run(['gzip', '-f', str(x)], check=True)


def position_class(r):
    v = r['verdict']
    if v == 'consistent':
        return 'supported'
    if v in ('contradicted', 'mixed'):
        return 'conflict'
    if v == 'insufficient':
        return 'locus_ambiguous' if int(r['locus_ambiguous']) > 0 else 'insufficient_coverage'
    return 'not_tested'


def pair_class(tier, a, b):
    """a, b: position classes (None = side without reads)."""
    got = [x for x in (a, b) if x is not None]
    if 'conflict' in got:
        return 'conflict_at_a_position'
    n = sum(x == 'supported' for x in got)
    if tier == 'C' or len(got) == 1:
        return 'single_side_supported' if n else 'single_side_not_supported'
    return {2: 'both_positions_supported', 1: 'one_position_supported', 0: 'no_position_supported'}[n]


SVG_COL = {'P': '#2a7ab0', 'A': '#c2452d', 'ambiguous': '#999999', 'locus_ambiguous': '#d49a1f', 'off_locus': '#bbbbbb',
           'cigar_fail': '#7b4fa0'}


def svg_position(p, segs, title, w=760, row=6):
    """One position: U / element / D bar and one line per genome alignment of reads that got a call there."""
    lo, hi = int(p['lo']), int(p['hi'])
    x0, x1 = lo - rc3.FLANK, hi + rc3.FLANK
    sx = lambda x: 10 + (min(max(x, x0), x1) - x0) * (w - 20) / (x1 - x0)
    segs = sorted(segs, key=lambda s: (s[0].split('_')[0], s[1]))[:120]
    h = 60 + row * len(segs)
    out = ['<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" font-family="sans-serif" font-size="11">' % (w, h),
           '<rect width="100%%" height="100%%" fill="white"/><text x="10" y="14">%s</text>' % title,
           '<rect x="%.1f" y="24" width="%.1f" height="10" fill="#cccccc"/>' % (sx(x0), sx(lo) - sx(x0)),
           '<rect x="%.1f" y="24" width="%.1f" height="10" fill="#333333"/>' % (sx(lo), sx(hi) - sx(lo)),
           '<rect x="%.1f" y="24" width="%.1f" height="10" fill="#cccccc"/>' % (sx(hi), sx(x1) - sx(hi)),
           '<text x="10" y="48">U (assembly)</text><text x="%.1f" y="48">element %d bp</text>' % (sx(lo), hi - lo),
           '<text x="%.1f" y="48">D (assembly)</text>' % (sx(hi) + 4)]
    for i, (call, ts, te, strand, tp, q) in enumerate(segs):
        y = 56 + i * row
        out.append('<line x1="%.1f" x2="%.1f" y1="%d" y2="%d" stroke="%s" stroke-width="3"%s/>' % (
            sx(ts), sx(te), y, y, SVG_COL.get(call.split('_matched')[0], '#000'), '' if tp == 'P' else ' stroke-dasharray="4,2"'))
    out.append('</svg>')
    return '\n'.join(out)


def cmd_summarize(a):
    rd = a.round_dir
    pairs = pc.read_tsv(rd / 'pairs.tsv')
    positions = {p['test_id']: p for p in pc.read_tsv(rd / 'positions.tsv')}
    res, runs_used, plot = {}, [], {}
    for s in sorted({p['sample'] for p in positions.values()}):
        od = rd / 'check' / s
        if not (od / 'positions_result.tsv').exists():
            continue
        for r in pc.read_tsv(od / 'positions_result.tsv'):
            res[r['test_id']] = r
        runs_used += [dict(x, sample=s) for x in json.loads((od / 'run.json').read_text())['runs']]
        plot.update(json.loads((od / 'plot.json').read_text()))
    out_pos = []
    for t, p in sorted(positions.items(), key=lambda x: (x[1]['sample'], x[1]['contig'], int(x[1]['lo']))):
        r = res.get(t, dict(verdict='not_checked'))
        out_pos.append(dict(p, **{k: r.get(k, '') for k in ('verdict', 'P', 'A', 'P_matched', 'A_matched', 'ambiguous',
                                                             'locus_ambiguous', 'off_locus', 'cigar_fail')},
                            position_class=position_class(r) if r['verdict'] not in ('not_checked',) else 'not_tested'))
    pcls = {p['test_id']: p['position_class'] for p in out_pos}
    want = {p['sample'] for p in positions.values()}
    out_pairs = []
    for r in pairs:
        cl = {}
        for s in 'ab':
            if r[s + '_sample'] in want:
                cl[s] = pcls[test_id(r[s + '_sample'], r[s + '_contig'], int(r[s + '_elem_lo']), int(r[s + '_elem_hi']))]
        out_pairs.append(dict(r, a_position_class=cl.get('a', ''), b_position_class=cl.get('b', ''),
                              pair_class=pair_class(r['tier'], cl.get('a'), cl.get('b'))))
    pc.write_tsv(a.out / 'positions.tsv.gz', out_pos, list(out_pos[0]))
    pc.write_tsv(a.out / 'pairs.tsv.gz', out_pairs, list(out_pairs[0]))
    pc.write_tsv(a.out / 'runs.tsv', runs_used, sorted({k for x in runs_used for k in x}))
    groups = defaultdict(Counter)
    for r in out_pairs:
        groups['%s|%s|%s' % (r['tier'], r['structure'], r['control_class'] or 'candidate')][r['pair_class']] += 1
    summ = dict(version=VERSION, caller=rc3.VERSION, evidence_added='long_read_position_check',
                evidence_source_of_pairs='assembly_only', production_ready=False,
                note='reads are the assembly\'s own; this checks the assembled U-core-D structure, not mobility or '
                     'independent insertion events',
                positions=len(out_pos), position_class=dict(Counter(p['position_class'] for p in out_pos)),
                position_verdict=dict(Counter(p['verdict'] for p in out_pos)),
                pairs=len(out_pairs), pair_class={k: dict(v) for k, v in sorted(groups.items())},
                runs=dict(Counter(r['status'] for r in runs_used)))
    (a.out / 'summary.json').write_text(json.dumps(summ, indent=1) + '\n')
    for f in ('selection.json', 'samples.txt'):
        (a.out / f).write_text((rd / f).read_text())
    figs = a.out / 'figures'
    figs.mkdir(exist_ok=True)
    for r in out_pairs:
        if r['pair_id'] not in a.figure_pairs.split(','):
            continue
        parts = []
        for s in 'ab':
            if r[s + '_sample'] not in want:
                continue
            t = test_id(r[s + '_sample'], r[s + '_contig'], int(r[s + '_elem_lo']), int(r[s + '_elem_hi']))
            parts.append(svg_position(positions[t], plot.get(t, []), '%s %s %s:%s-%s  %s' % (
                r['pair_id'], r[s + '_sample'], r[s + '_contig'], r[s + '_elem_lo'], r[s + '_elem_hi'], pcls[t])))
        (figs / ('%s.svg' % r['pair_id'])).write_text('\n'.join(parts) + '\n')
    print(json.dumps(summ, indent=1))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('mode', choices=['pairs', 'select', 'genome', 'check', 'summarize'])
    p.add_argument('--work', type=Path, default=PROJ / 'core_context' / 'scale_comamonadaceae_v0.1')
    p.add_argument('--repo-out', type=Path, default=Path(__file__).resolve().parent / 'scale_comamonadaceae_v0.1')
    p.add_argument('--batches', default='1,2,3')
    p.add_argument('--round-dir', type=Path)
    p.add_argument('--assemblies', default='')
    p.add_argument('--sample')
    p.add_argument('--threads', type=int, default=16)
    p.add_argument('--out', type=Path)
    p.add_argument('--figure-pairs', default='')
    a = p.parse_args()
    dict(pairs=cmd_pairs, select=cmd_select, genome=cmd_genome, check=cmd_check, summarize=cmd_summarize)[a.mode](a)


if __name__ == '__main__':
    main()
