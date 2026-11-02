"""Read-based site search, pilot (design: scale_comamonadaceae_v0.1/READ_SITE_SEARCH_PILOT.md).

Given known cores, look for them in ALL reads of a sample (not only reads near catalog positions), take the read
sequence on both sides of each core hit, place both flanks on the sample's own assembly, merge observations of one
position by molecule, and compare the merged positions with the assembly catalog. No empty site is needed, and one
read hit is never one site.

freeze   Pilot samples by a fixed rule (written to <out>/samples.tsv and freeze.json before any search):
           eligible = every checked run counted in molecules, read type ONT or PacBio HiFi/CCS, no run submitted as
           an aligned BAM, depth >= 30x (verified bases / assembly length); per genus the eligible sample with the
           most candidate cores supported in read-check round 5 (ties: depth, then id); genera whose best sample has
           >= 15 such cores; then HiFi samples by the same score until at least 2 HiFi samples are in.
panel    Core panel: the distinct core DNA sequences (coding span, in core orientation) of every occurrence in the
         pool of the cores with an A/B pair supported on both sides in round 5. Catalog: every reviewed occurrence
         (any core) of the pilot samples. Negative control: per pilot sample, single-copy conserved genes (fixed
         product list) that overlap no catalog occurrence.
search   One sample: all reads of every verified run -> panel (minimap2, read preset); core hits; flanks of exactly
         L bp next to each hit (L in FLANKS); flanks -> the sample's own assembly; molecule grouping as in the read
         check; merged positions; classes; per (L, N) counts. Never downloads (ready manifest only).
combine  Pilot tables and summary over the samples.

Classes of a merged position (both flanks placed, each without a competing placement within MIN_DELTA_AS):
  known_catalog_position          a catalog occurrence of this core lies between the two flanks
  known_catalog_position_other_core   a catalog occurrence of another core lies there (homologous core)
  catalog_omission                no catalog occurrence there, but the assembly sequence between the flanks
                                  contains the core (aligned): the assembly has it, the catalog does not list it
  not_in_assembly:core_absent     the flanks are adjacent in the assembly (or the sequence between is much shorter
                                  than the core and does not contain it): reads carry the core at a junction the
                                  assembly shows without it
  not_in_assembly:other_sequence  sequence between the flanks of core-like length that does not contain the core
  not_in_assembly:junction        flanks on different contigs, opposite strands, or far apart: the reads join two
                                  assembly regions with the core between
Not positions: unplaceable_context (a flank unplaced or placed ambiguously; reported by side), one_sided (a flank
shorter than L). Every merged position is reported with its number of independent molecules; a candidate needs
>= N (pilot N = 3, sensitivity N in COUNTS).
"""
import argparse
import gzip
import json
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

import protein_context as pc
import scale_entry as se  # noqa: F401  (puts fna_pool on sys.path)
import s3_readcheck as rc3
import ucd_readcheck as u
from scale_candidates import POOL

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'scripts'))
import reads_ready  # noqa: E402

VERSION = 'read_site_search_v0.1'
WORK = u.PROJ / 'core_context' / 'read_site_search' / 'pilot_v0.1'
REPO = Path(__file__).resolve().parent / 'scale_comamonadaceae_v0.1'
ROUND5 = REPO / 'readcheck_round5'
FLANKS = (500, 1000, 2000)  # bp of read sequence taken next to the core on each side
MAIN_FLANK, MAIN_COUNT = 1000, 3
COUNTS = (1, 2, 3, 5)  # independent molecules a merged position needs
HIT_COVERAGE = 0.8  # fraction of the panel sequence a core hit must cover
PLACED_FRACTION = 0.8  # fraction of a flank its best assembly alignment must cover
TOL = rc3.LOCUS_TOL  # bp; placements of one position agree within this
MIN_DELTA_AS = rc3.MIN_DELTA_AS
MM_PANEL = ['-c', '--secondary=yes', '-N', '200', '-p', '0.5']
MM_FLANK = ['-c', '--secondary=yes', '-N', '10', '-p', '0.5']
HIFI = 'pacbio_hifi_ccs'
NEG_PRODUCTS = ('DNA-directed RNA polymerase subunit beta', "DNA-directed RNA polymerase subunit beta'",
                'DNA gyrase subunit B', 'recombinase RecA', 'DNA polymerase III subunit beta',
                '50S ribosomal protein L2', '30S ribosomal protein S3', 'signal recognition particle protein')


def manifest():
    return {r['sample']: r for r in pc.read_tsv(REPO / 'manifest.tsv')}


def cmd_freeze(a):
    a.out.mkdir(parents=True, exist_ok=True)
    if (a.out / 'samples.tsv').exists():
        sys.exit('frozen already: %s' % (a.out / 'samples.tsv'))
    man, meta = manifest(), {r['run_accession']: r for r in pc.read_tsv(POOL / 'runs.tsv')}
    runs = defaultdict(list)
    for r in pc.read_tsv(ROUND5 / 'runs.tsv'):
        if not r['sample'].endswith(u.ROT):
            runs[r['sample']].append(r)
    pos = {p['test_id']: p for p in pc.read_tsv(ROUND5 / 'positions.tsv.gz')}
    cores = defaultdict(set)
    for m in pc.read_tsv(ROUND5 / 'interval_map.tsv.gz'):
        if u.pair_role(m) == 'candidate' and pos[m['test_id']]['position_class'] == 'supported':
            cores[pos[m['test_id']]['sample']].add((m['structure'], m['core_key']))
    info = u.run_info()
    rows = []
    for s, rs in sorted(runs.items()):
        m = man[s]
        files = {x['run']: x for x in info.get(s, [])}
        bases = sum(json.loads((Path(files[r['run']]['local_file']).parent / 'verify.json').read_text())
                    ['summary']['bases'] for r in rs if r['run'] in files)
        sub = [meta[r['run']].get('submitted_names_sra', '') for r in rs]
        types = sorted({r['read_type'] for r in rs})
        why = [w for w, bad in (
            ('not_counted_in_molecules', any(r['count_unit'] != 'molecules' for r in rs)),
            ('read_type', any(t not in ('ont', HIFI) for t in types)),
            ('aligned_bam_submission', any('align' in x.lower() for x in sub)),
            ('depth_below_30x', bases < 30 * int(m['total_bp'])),
            ('run_not_ready', any(r['run'] not in files for r in rs))) if bad]
        rows.append(dict(sample=s, genus=m['genus'], organism=m['organism'], assembly=m['assembly'],
                         assembly_level=m['assembly_level'], read_type='/'.join(types),
                         runs=';'.join(r['run'] for r in rs), depth=round(bases / int(m['total_bp']), 1),
                         supported_candidate_cores=len(cores[s]), eligible=not why, not_eligible=';'.join(why)))
    el = sorted((r for r in rows if r['eligible']), key=lambda r: (-r['supported_candidate_cores'], -r['depth'],
                                                                   r['sample']))
    pick, seen = [], set()
    for r in el:
        if r['genus'] not in seen:
            seen.add(r['genus'])
            if r['supported_candidate_cores'] >= 15:
                pick.append(dict(r, reason='best_of_genus'))
    for r in el:
        if sum(HIFI in x['read_type'] for x in pick) >= 2:
            break
        if HIFI in r['read_type'] and r['sample'] not in {x['sample'] for x in pick}:
            pick.append(dict(r, reason='hifi_minimum'))
    for r in pick:
        r['sha256'] = ';'.join(x['sha256'] for x in info[r['sample']])
    cols = ['sample', 'genus', 'organism', 'assembly', 'assembly_level', 'read_type', 'runs', 'depth',
            'supported_candidate_cores', 'reason', 'sha256']
    pc.write_tsv(a.out / 'samples.tsv', pick, cols)
    pc.write_tsv(a.out / 'candidates_considered.tsv', rows, [c for c in rows[0]])
    (a.out / 'freeze.json').write_text(json.dumps(dict(
        version=VERSION, rule=__doc__.split('freeze   ')[1].split('panel    ')[0].split(), samples=[r['sample'] for r in pick],
        considered=len(rows), eligible=len(el), round5=str(ROUND5)), indent=1) + '\n')
    for r in pick:
        print(r['sample'], r['genus'], r['read_type'], r['depth'], r['supported_candidate_cores'], r['reason'])


def pilot_samples(out):
    return [r['sample'] for r in pc.read_tsv(out / 'samples.tsv')]


def cmd_panel(a):
    """Core panel FASTA (distinct sequences), member table, catalog of the pilot samples, negative controls."""
    wd = a.work / 'panel'
    wd.mkdir(parents=True, exist_ok=True)
    pilot = set(pilot_samples(a.out))
    man = manifest()
    pairs = pc.read_tsv(ROUND5 / 'pairs.tsv.gz')
    cores = {(p['structure'], p['core_key']) for p in pairs if p['tier'] in 'AB' and p['tier']
             and u.pair_role(p) == 'candidate' and p['pair_class'] == 'both_positions_supported'}
    occ_dir = u.PROJ / 'core_context' / 'scale_comamonadaceae_v0.1' / 'review' / 's6' / 'b1-2-3'
    want, catalog = defaultdict(list), []
    for k in range(64):
        for r in pc.read_tsv(occ_dir / ('shard_%03d' % k) / 'occurrences.tsv.gz'):
            if (r['structure'], r['core_key']) in cores:
                want[r['sample']].append(r)
            if r['sample'] in pilot:
                catalog.append({c: r[c] for c in ('structure', 'core_key', 'occurrence_id', 'sample', 'contig',
                                                  'start0', 'end0', 'strand')})
    member_of, members = {}, []
    for s in sorted(want):
        seqs = rc3.read_fasta(man[s]['fasta'])
        for r in want[s]:
            d = seqs[r['contig']][int(r['start0']):int(r['end0'])].upper()
            d = d if r['strand'] == '+' else rc3.rc(d)
            key = (r['structure'], r['core_key'], d)
            if key not in member_of:
                member_of[key] = 'core|%s|%s|m%d' % (r['structure'], r['core_key'], len(members))
                members.append(dict(member=member_of[key], structure=r['structure'], core_key=r['core_key'],
                                    length=len(d), occurrences=0, first_occurrence=r['occurrence_id']))
            members[int(member_of[key].rsplit('|m', 1)[1])]['occurrences'] += 1
    with open(str(wd / 'panel.fa'), 'w') as h:
        for (st, ck, d), name in member_of.items():
            h.write('>%s\n%s\n' % (name, d))
    neg = []
    occ_by = defaultdict(list)
    for r in catalog:
        occ_by[(r['sample'], r['contig'])].append((int(r['start0']), int(r['end0'])))
    with open(str(wd / 'panel.fa'), 'a') as h:
        for s in sorted(pilot):
            genes = [g for g in pc.read_tsv(u.PROJ / 'core_context' / 'scale_comamonadaceae_v0.1' / 'genes' /
                                            ('%s.tsv.gz' % s)) if g['type'] == 'CDS']
            n = Counter(g['product'] for g in genes)
            seqs = rc3.read_fasta(man[s]['fasta'])
            for g in genes:
                lo, hi = int(g['start0']), int(g['end0'])
                if (g['product'] in NEG_PRODUCTS and n[g['product']] == 1 and g['partial'] == 'False'
                        and g['pseudo'] == 'False' and g['complete_cds'] == 'True'
                        and not any(x < hi and y > lo for x, y in occ_by[(s, g['contig'])])):
                    key = 'NEG_' + g['product'].replace(' ', '_').replace("'", 'p')
                    d = seqs[g['contig']][lo:hi].upper()
                    d = d if g['strand'] == '+' else rc3.rc(d)
                    name = 'neg|single_gene|%s|%s' % (key, s)
                    h.write('>%s\n%s\n' % (name, d))
                    neg.append(dict(member=name, structure='single_gene', core_key=key, length=len(d), occurrences=1,
                                    first_occurrence='%s|%s|%d' % (s, g['contig'], lo)))
                    catalog.append(dict(structure='single_gene', core_key=key, occurrence_id='neg|%s|%s|%d' % (
                        s, g['contig'], lo), sample=s, contig=g['contig'], start0=lo, end0=hi, strand=g['strand']))
    pc.write_tsv(wd / 'members.tsv', members + neg, list(members[0]))
    pc.write_tsv(wd / 'catalog.tsv', catalog, list(catalog[0]))
    summ = dict(version=VERSION, cores=len(cores), panel_members=len(members),
                panel_occurrences=sum(m['occurrences'] for m in members),
                panel_bp=sum(m['length'] for m in members), negative_control_members=len(neg),
                negative_control_products=dict(Counter(m['core_key'] for m in neg)),
                catalog_occurrences_pilot=len(catalog), pilot_samples=sorted(pilot))
    (wd / 'panel.json').write_text(json.dumps(summ, indent=1) + '\n')
    print(json.dumps(summ, indent=1))


def paf_rows(path):
    op = gzip.open if str(path).endswith('.gz') else open
    with op(str(path), 'rt') as h:
        for line in h:
            f = line.rstrip('\n').split('\t')
            tags = dict(x.split(':', 2)[::2] for x in f[12:])
            yield dict(q=f[0], qlen=int(f[1]), qs=int(f[2]), qe=int(f[3]), strand=f[4], t=f[5], tlen=int(f[6]),
                       ts=int(f[7]), te=int(f[8]), match=int(f[9]), block=int(f[10]), score=int(tags.get('AS', 0)),
                       tp=tags.get('tp', 'P'))


def core_hits(paf, coverage):
    """Per read: non-overlapping best hits per core key covering >= coverage of a panel sequence."""
    by = defaultdict(list)
    for b in paf_rows(paf):
        if (b['te'] - b['ts']) / b['tlen'] >= coverage:
            kind, st, ck, _ = b['t'].split('|')
            by[(b['q'], st, ck)].append(dict(b, kind=kind))
    hits = []
    for (q, st, ck), bs in by.items():
        keep = []
        for b in sorted(bs, key=lambda x: (-x['score'], x['t'], x['qs'])):
            if all(min(b['qe'], k['qe']) - max(b['qs'], k['qs']) <= 0.5 * min(b['qe'] - b['qs'], k['qe'] - k['qs'])
                   for k in keep):
                keep.append(b)
        for b in keep:
            hits.append(dict(read=q, qlen=b['qlen'], structure=st, core_key=ck, kind=b['kind'], member=b['t'],
                             qs=b['qs'], qe=b['qe'], strand=b['strand'], coverage=round((b['te'] - b['ts']) / b['tlen'], 3),
                             identity=round(b['match'] / b['block'], 3), score=b['score']))
    return hits


def placement(alns, side):
    """Best placement of one flank: (status, contig, inner position, strand). The inner position is the assembly
    coordinate of the flank end next to the core (left flank: its last base; right flank: its first base)."""
    if not alns:
        return 'unplaced', '', -1, ''
    best = max(alns, key=lambda g: g['score'])
    if (best['qe'] - best['qs']) < PLACED_FRACTION * best['qlen']:
        return 'unplaced', '', -1, ''

    def inner(g):
        return (g['te'] if g['strand'] == '+' else g['ts']) if side == 'L' else (g['ts'] if g['strand'] == '+' else g['te'])
    for g in alns:
        if g is not best and g['score'] >= best['score'] - MIN_DELTA_AS and (
                g['t'] != best['t'] or g['strand'] != best['strand'] or abs(inner(g) - inner(best)) > TOL):
            return 'ambiguous', best['t'], inner(best), best['strand']
    return 'placed', best['t'], inner(best), best['strand']


def canonical(o):
    """Flank placements in assembly orientation: an observation whose left flank (read-forward) maps to the minus
    strand is the same position read from the other strand, so the sides are swapped and the strands flipped.
    Reads from both strands of one position then merge."""
    if o['left_strand'] != '-':
        return o['left_contig'], o['left_pos'], o['left_strand'], o['right_contig'], o['right_pos'], o['right_strand']
    flip = {'+': '-', '-': '+'}
    return o['right_contig'], o['right_pos'], flip[o['right_strand']], o['left_contig'], o['left_pos'], '+'


def geometry(lc, lp, ls, rc, rp, rs, core_len):
    """('proper', assembly start, end of the sequence between the flanks) or ('junction', ...)."""
    if lc != rc or ls != rs:
        return 'junction', -1, -1
    a, b = (lp, rp) if ls == '+' else (rp, lp)
    if b - a < -TOL or b - a > core_len + 10000:
        return 'junction', -1, -1
    return 'proper', a, b


def cmd_search(a):
    s = a.sample
    od = a.work / 'samples' / s
    od.mkdir(parents=True, exist_ok=True)
    runs = u.run_info()[s]
    reads_ready.require([r['run'] for r in runs])
    fasta = manifest()[s]['fasta']
    panel = a.work / 'panel' / 'panel.fa'
    hits, recs, hit_reads = [], [], set()
    for r in runs:
        paf = od / ('panel_%s.paf.gz' % r['run'])
        if not paf.exists():
            cmd = 'minimap2 -x %s %s -t %d %s %s | gzip -1 > %s.tmp' % (
                rc3.PRESET[r['read_type']], ' '.join(MM_PANEL), a.threads, panel, r['local_file'], paf)
            subprocess.run(['bash', '-o', 'pipefail', '-c', cmd], check=True, stderr=subprocess.DEVNULL)
            Path(str(paf) + '.tmp').rename(paf)
        for h in core_hits(paf, HIT_COVERAGE):
            hits.append(dict(h, run=r['run'], hit_id='%s|%s|%d' % (h['read'], h['core_key'], h['qs'])))
            hit_reads.add(h['read'])
    # flanks: exact L bp next to each hit, in read-forward orientation; molecule ids from the FASTQ headers
    seq_of, preset = {}, {}
    for r in runs:
        with gzip.open(r['local_file'], 'rt') as h:
            name = None
            for i, line in enumerate(h):
                if i % 4 == 0:
                    name = line[1:].split(None, 1)[0]
                    if name in hit_reads:
                        recs.append((name, rc3.record_ids(line[1:].strip(), r['grouping'])))
                        preset[name] = rc3.PRESET[r['read_type']]
                elif i % 4 == 1 and name in hit_reads:
                    seq_of[name] = line.strip()
    mol = rc3.molecules(recs)
    fl = od / 'flanks.fa'
    with open(str(fl), 'w') as h:
        for x in hits:
            sq = seq_of[x['read']]
            for L in FLANKS:
                left, right = sq[max(0, x['qs'] - L):x['qs']], sq[x['qe']:x['qe'] + L]
                if len(left) == L and len(right) == L:
                    h.write('>%s|%d|L\n%s\n>%s|%d|R\n%s\n' % (x['hit_id'], L, left, x['hit_id'], L, right))
    presets = sorted(set(preset.values())) or ['map-ont']
    fpaf = od / 'flanks.paf'
    with open(str(fpaf), 'w') as h:  # one preset per sample in the pilot (single-run samples)
        subprocess.run(['minimap2', '-x', presets[0]] + MM_FLANK + ['-t', str(a.threads), fasta, str(fl)],
                       stdout=h, stderr=subprocess.DEVNULL, check=True)
    alns = defaultdict(list)
    for g in paf_rows(fpaf):
        alns[g['q']].append(g)
    seqs = rc3.read_fasta(fasta)
    catalog = [c for c in pc.read_tsv(a.work / 'panel' / 'catalog.tsv') if c['sample'] == s]
    member_seq = panel_members(panel, {x['member'] for x in hits})
    obs = []
    for x in hits:
        for L in FLANKS:
            o = dict(x, molecule=mol.get(x['read'], x['read']), flank=L)
            k = '%s|%d' % (x['hit_id'], L)
            if not flank_written(x, seq_of, L):
                o.update(status='one_sided')
            else:
                pl = placement(alns.get(k + '|L', []), 'L')
                pr = placement(alns.get(k + '|R', []), 'R')
                o.update(left_status=pl[0], left_contig=pl[1], left_pos=pl[2], left_strand=pl[3],
                         right_status=pr[0], right_contig=pr[1], right_pos=pr[2], right_strand=pr[3])
                o['status'] = 'placed' if pl[0] == pr[0] == 'placed' else 'unplaceable_context'
            obs.append(o)
    positions = merge(obs, catalog, seqs, member_seq, od, seq_of, presets[0])
    cols = ['flank', 'structure', 'core_key', 'kind', 'cls', 'raw_cls', 'molecules', 'read_records', 'contig', 'start', 'end',
            'left_contig', 'left_pos', 'left_strand', 'right_contig', 'right_pos', 'right_strand', 'catalog',
            'core_in_assembly', 'members', 'identity_median', 'runs']
    pc.write_tsv(od / 'positions.tsv.gz', positions, cols)
    ocols = ['flank', 'run', 'read', 'molecule', 'structure', 'core_key', 'kind', 'member', 'qs', 'qe', 'strand',
             'coverage', 'identity', 'status', 'left_status', 'left_contig', 'left_pos', 'left_strand', 'right_status',
             'right_contig', 'right_pos', 'right_strand']
    pc.write_tsv(od / 'observations.tsv.gz', obs, ocols)
    subprocess.run(['gzip', '-f', str(fpaf)], check=True)
    fl.unlink()
    (od / 'search.json').write_text(json.dumps(dict(
        version=VERSION, sample=s, runs=[dict(run=r['run'], read_type=r['read_type'], grouping=r['grouping'],
                                              sha256=r['sha256']) for r in runs],
        panel_hits=len(hits), reads_with_hits=len(hit_reads), molecules_with_hits=len(set(mol.values())),
        hit_coverage=HIT_COVERAGE, flanks=FLANKS, minimap2_panel=MM_PANEL, minimap2_flank=MM_FLANK,
        observation_status={str(L): dict(Counter(o['status'] for o in obs if o['flank'] == L)) for L in FLANKS}),
        indent=1) + '\n')


def flank_written(x, seq_of, L):
    n = len(seq_of[x['read']])
    return x['qs'] >= L and n - x['qe'] >= L


def panel_members(panel, names):
    out, name = {}, None
    with open(str(panel)) as h:
        for line in h:
            if line.startswith('>'):
                name = line[1:].strip()
            elif name in names:
                out[name] = line.strip()
    return out


def merge(obs, catalog, seqs, member_seq, od, seq_of, preset):
    """Placed observations -> one row per merged position (per flank length); unplaceable ones -> one row per
    (core, placed side anchor, reason). Molecules are counted once per position."""
    rows = []
    for L in FLANKS:
        groups = defaultdict(list)
        for o in obs:
            if o['flank'] == L and o['status'] == 'placed':
                lc, lp, ls, rc, rp, rs = canonical(o)
                o.update(c_left_pos=lp, c_right_pos=rp)
                groups[(o['structure'], o['core_key'], lc, ls, rc, rs)].append(o)
        placed = []
        for key, os_ in groups.items():
            for cl in cluster(os_):
                placed.append(position_row(L, key, cl, catalog, seqs, member_seq, od, seq_of, preset))
        rows += mark_shifted(placed)
        un = defaultdict(list)
        for o in obs:
            if o['flank'] == L and o['status'] == 'unplaceable_context':
                why = 'left_%s|right_%s' % (o['left_status'], o['right_status'])
                side = 'L' if o['left_status'] == 'placed' else 'R' if o['right_status'] == 'placed' else ''
                anchor = (o['left_contig'], o['left_pos'] // (2 * TOL)) if side == 'L' else (
                    o['right_contig'], o['right_pos'] // (2 * TOL)) if side == 'R' else ('', 0)
                un[(o['structure'], o['core_key'], why, side, anchor)].append(o)
        for (st, ck, why, side, anchor), os_ in un.items():
            rows.append(dict(flank=L, structure=st, core_key=ck, kind=os_[0]['kind'], cls='unplaceable_context:' + why,
                             molecules=len({o['molecule'] for o in os_}), read_records=len(os_), contig=anchor[0],
                             start=anchor[1] * 2 * TOL if side else '', end='', members=len({o['member'] for o in os_}),
                             identity_median=median([o['identity'] for o in os_]),
                             runs=';'.join(sorted({o['run'] for o in os_}))))
        one = Counter((o['structure'], o['core_key'], o['kind']) for o in obs if o['flank'] == L and o['status'] == 'one_sided')
        mols = defaultdict(set)
        for o in obs:
            if o['flank'] == L and o['status'] == 'one_sided':
                mols[(o['structure'], o['core_key'], o['kind'])].add(o['molecule'])
        for (st, ck, kind), n in one.items():
            rows.append(dict(flank=L, structure=st, core_key=ck, kind=kind, cls='one_sided',
                             molecules=len(mols[(st, ck, kind)]), read_records=n))
    return rows


def median(xs):
    xs = sorted(xs)
    return xs[len(xs) // 2] if xs else ''


def cluster(os_):
    """Single linkage on (left position, right position) within TOL."""
    os_ = sorted(os_, key=lambda o: (o['c_left_pos'], o['c_right_pos']))
    out = []
    for o in os_:
        for c in out:
            if any(abs(o['c_left_pos'] - p['c_left_pos']) <= TOL and abs(o['c_right_pos'] - p['c_right_pos']) <= TOL
                   for p in c):
                c.append(o)
                break
        else:
            out.append([o])
    return out


def position_row(L, key, cl, catalog, seqs, member_seq, od, seq_of, preset):
    st, ck, lc, ls, rc, rs = key
    lp = int(median([o['c_left_pos'] for o in cl]))
    rp = int(median([o['c_right_pos'] for o in cl]))
    core_len = int(median([o['qe'] - o['qs'] for o in cl]))
    geo, x, y = geometry(lc, lp, ls, rc, rp, rs, core_len)
    cat, inside = '', ''
    if geo == 'junction':
        cls = 'not_in_assembly:junction'
    else:
        here = [c for c in catalog if c['contig'] == lc and int(c['start0']) >= x - TOL and int(c['end0']) <= y + TOL]
        same = [c for c in here if c['core_key'] == ck]
        cat = ';'.join(sorted(c['occurrence_id'] for c in (same or here)))
        if same:
            cls = 'known_catalog_position'
        elif here:
            cls = 'known_catalog_position_other_core'
        else:
            seg = seqs[lc][max(0, x):max(0, y)]
            if len(seg) >= 0.5 * core_len:
                # does the assembly between the flanks hold the core? the panel sequence (may be from another genus)
                # or the core region of up to 3 of this position's own reads
                mem = Counter(o['member'] for o in cl).most_common(1)[0][0]
                own = sorted(cl, key=lambda o: (-o['identity'], o['read']))[:3]
                inside = ('panel_sequence' if contains(seg, member_seq[mem], od, 'asm20') else
                          'own_reads' if any(contains(seg, seq_of[o['read']][o['qs']:o['qe']], od, preset) for o in own)
                          else '')
            cls = ('catalog_omission' if inside else 'not_in_assembly:other_sequence'
                   if len(seg) >= 0.5 * core_len else 'not_in_assembly:core_absent')
    return dict(flank=L, structure=st, core_key=ck, kind=cl[0]['kind'], cls=cls,
                molecules=len({o['molecule'] for o in cl}), read_records=len(cl), contig=lc if geo == 'proper' else '',
                start=x if geo == 'proper' else '', end=y if geo == 'proper' else '', left_contig=lc, left_pos=lp,
                left_strand=ls, right_contig=rc, right_pos=rp, right_strand=rs, catalog=cat, core_in_assembly=inside,
                members=len({o['member'] for o in cl}), identity_median=median([o['identity'] for o in cl]),
                runs=';'.join(sorted({o['run'] for o in cl})))


def mark_shifted(rows):
    """A non-known position whose assembly interval (or, for a junction, either flank anchor) overlaps a known
    position of the same core or a homologous one in the same sample is the known position seen with shifted
    anchors (a core hit that covers only part of the gene moves a flank): cls shifted_known_position, raw_cls
    keeps the class it had. Found with the negative-control genes of the pilot."""
    known = defaultdict(list)
    for r in rows:
        if r['cls'].startswith('known_catalog_position'):
            known[r['contig']].append((r['start'], r['end']))
    for r in rows:
        if r['cls'].startswith('known_catalog_position'):
            continue
        hit = (any(a < r['end'] + TOL and b > r['start'] - TOL for a, b in known[r['contig']]) if r['contig'] else
               any(a - TOL <= p <= b + TOL for c, p in ((r['left_contig'], r['left_pos']), (r['right_contig'], r['right_pos']))
                   for a, b in known[c]))
        if hit:
            r['raw_cls'], r['cls'] = r['cls'], 'shifted_known_position'
    return rows


def contains(seg, core, od, preset):
    """The assembly segment contains the core: one minimap2 alignment covering >= HIT_COVERAGE of it."""
    if not seg or not core:
        return False
    ref, qry = od / 'seg.fa', od / 'core.fa'
    ref.write_text('>seg\n%s\n' % seg)
    qry.write_text('>core\n%s\n' % core)
    res = subprocess.run(['minimap2', '-x', preset, '-c', str(ref), str(qry)], stdout=subprocess.PIPE,
                         stderr=subprocess.DEVNULL, universal_newlines=True, check=True)
    ref.unlink()
    qry.unlink()
    return any((int(f[3]) - int(f[2])) / int(f[1]) >= HIT_COVERAGE for f in
               (l.split('\t') for l in res.stdout.splitlines()))


def cmd_combine(a):
    samples = pilot_samples(a.out)
    rows, runs = [], {}
    for s in samples:
        od = a.work / 'samples' / s
        for r in pc.read_tsv(od / 'positions.tsv.gz'):
            rows.append(dict(sample=s, **r))
        runs[s] = json.loads((od / 'search.json').read_text())
    catalog = pc.read_tsv(a.work / 'panel' / 'catalog.tsv')
    members = pc.read_tsv(a.work / 'panel' / 'members.tsv')
    panel_cores = {(m['structure'], m['core_key']) for m in members}
    # every (sample, negative-control gene) whose reads were searched: the sample's own gene, or (when its own copy is
    # a catalog occurrence and so left out) the copies of other samples
    neg_genes = {(s, m['core_key']) for m in members if m['member'].startswith('neg|') for s in samples}
    r5 = {}
    for p in pc.read_tsv(ROUND5 / 'interval_map.tsv.gz'):
        r5.setdefault(p['occurrence'], set()).add(p['position_class'])
    summ = dict(version=VERSION, samples=samples, main=dict(flank=MAIN_FLANK, molecules=MAIN_COUNT), grid={})
    for L in FLANKS:
        for N in COUNTS:
            sel = [r for r in rows if int(r['flank']) == L]
            pos = [r for r in sel if not r['cls'].startswith(('unplaceable', 'one_sided', 'shifted'))]
            ok = [r for r in pos if int(r['molecules']) >= N]
            g = {kind: dict(positions_by_class=dict(Counter(r['cls'] for r in ok if r['kind'] == kind)),
                            positions_below_N=dict(Counter(r['cls'] for r in pos if r['kind'] == kind
                                                           and int(r['molecules']) < N)),
                            shifted_known_positions=sum(r['cls'] == 'shifted_known_position' and r['kind'] == kind
                                                        and int(r['molecules']) >= N for r in sel),
                            unplaceable_contexts=dict(Counter(r['cls'] for r in sel if r['kind'] == kind and
                                                              r['cls'].startswith('unplaceable')
                                                              and int(r['molecules']) >= N)))
                 for kind in ('core', 'neg')}
            g.update(recovery(ok, catalog, panel_cores, samples, r5))
            per_gene = Counter((r['sample'], r['core_key']) for r in ok if r['kind'] == 'neg')
            g['negative_control'] = dict(
                genes_searched=len(neg_genes), positions_per_gene=dict(sorted(Counter(
                    per_gene[k] for k in neg_genes).items())),
                extra_positions=[dict(sample=r['sample'], gene=r['core_key'], cls=r['cls'], molecules=r['molecules'],
                                      contig=r['left_contig'], left_pos=r['left_pos'], right_contig=r['right_contig'],
                                      right_pos=r['right_pos'])
                                 for r in ok if r['kind'] == 'neg' and per_gene[(r['sample'], r['core_key'])] > 1])
            summ['grid']['L%d_N%d' % (L, N)] = g
    out = a.out
    pc.write_tsv(out / 'positions.tsv.gz', rows, ['sample'] + [c for c in rows[0] if c != 'sample'])
    main = [r for r in rows if int(r['flank']) == MAIN_FLANK]
    for cls, name in (('catalog_omission', 'catalog_omissions'), ('not_in_assembly', 'not_in_assembly'),
                      ('unplaceable_context', 'unplaceable_contexts')):
        sel = [r for r in main if r['cls'].startswith(cls)]
        if sel:
            pc.write_tsv(out / ('%s.tsv' % name), sorted(sel, key=lambda r: (-int(r['molecules']), r['sample'])),
                         list(sel[0]))
    summ['search'] = {s: {k: v for k, v in j.items() if k in ('panel_hits', 'reads_with_hits', 'molecules_with_hits',
                                                             'observation_status')} for s, j in runs.items()}
    (out / 'summary.json').write_text(json.dumps(summ, indent=1) + '\n')
    print(json.dumps(summ['grid']['L%d_N%d' % (MAIN_FLANK, MAIN_COUNT)], indent=1))


def recovery(ok, catalog, panel_cores, samples, r5):
    """Known catalog positions of panel cores in the pilot samples recovered as a known_catalog_position of the same
    or of a homologous core."""
    got = {o for r in ok if r['cls'].startswith('known_catalog_position') for o in r['catalog'].split(';')}
    tot = Counter()
    for c in catalog:
        if c['sample'] in samples and (c['structure'], c['core_key']) in panel_cores:
            kind = 'negative_control' if c['occurrence_id'].startswith('neg|') else 'core'
            r5c = 'round5_supported' if 'supported' in r5.get(c['occurrence_id'], ()) else 'not_round5_supported'
            tot[(kind, r5c, 'total')] += 1
            tot[(kind, r5c, 'recovered')] += c['occurrence_id'] in got
    return dict(known_recovery={'%s|%s' % k[:2]: '%d/%d' % (tot[k[:2] + ('recovered',)], tot[k[:2] + ('total',)])
                                for k in sorted(tot) if k[2] == 'total'})


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('mode', choices=['freeze', 'panel', 'search', 'combine'])
    p.add_argument('--work', type=Path, default=WORK)
    p.add_argument('--out', type=Path, default=REPO / 'read_site_search_pilot')
    p.add_argument('--sample')
    p.add_argument('--threads', type=int, default=16)
    a = p.parse_args()
    dict(freeze=cmd_freeze, panel=cmd_panel, search=cmd_search, combine=cmd_combine)[a.mode](a)


if __name__ == '__main__':
    main()
