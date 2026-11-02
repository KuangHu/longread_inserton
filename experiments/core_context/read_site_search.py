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

Flank placement (placement): the best assembly alignment of a flank must cover >= PLACED_FRACTION of it AND reach
its inner end (the end next to the core) within INNER_CLIP bp; otherwise the inner coordinate is not observed and the
flank is inner_end_unaligned. A competing alignment within MIN_DELTA_AS elsewhere makes it ambiguous.

Catalog matching (match_catalog), since v0.2: an occurrence matches a merged position only when all three hold
  core relation  same core key, or a homologous core (a shared protein family or homology group, see core_relation);
                 an occurrence of an unrelated core is never a match, whatever its overlap
  orientation    the occurrence strand equals the core's orientation in the assembly as the reads show it (read hit
                 strand and flank strand), converted through the shared family's sign for a homologous core; when no
                 sign can be derived the orientation is undetermined and there is no match
  coverage       the overlap is >= HIT_COVERAGE of the occurrence and >= HIT_COVERAGE of the interval between the flanks
(v0.1 required only that an occurrence of any core lie wholly inside the interval.)

Classes of a merged position (both flanks placed, each without a competing placement within MIN_DELTA_AS):
  known_catalog_position          a catalog occurrence of this core matches
  known_catalog_position_other_core   a catalog occurrence of a homologous core matches
  catalog_mismatch:<reason>       an occurrence of this or a homologous core overlaps but fails orientation or
                                  coverage (reason: orientation, orientation_undetermined, coverage); neither known nor
                                  a candidate
  catalog_omission / homolog_not_in_catalog   no catalog occurrence there, but the assembly sequence between the
                                  flanks contains the core (the panel sequence or the position's own reads, aligned);
                                  split by the annotated CDS there (see annotate_omissions): only a CDS of the core's
                                  own family, or no CDS, is a catalog omission
  not_in_assembly:core_absent     the flanks are adjacent in the assembly (or the sequence between is much shorter
                                  than the core and does not contain it): reads carry the core at a junction the
                                  assembly shows without it
  not_in_assembly:other_sequence  sequence between the flanks of core-like length that does not contain the core
  not_in_assembly:junction        flanks on different contigs, opposite strands, or far apart: the reads join two
                                  assembly regions with the core between
  contig_end_join                 a junction whose flanks sit at the end and at the start of the same contig: the
                                  reads join the contig's end to its start through the core (the assembly breaks
                                  there); not a structure missing from the assembly
  shifted_known_position          overlaps a known position of the sample whose core is the same or homologous, in a
                                  consistent orientation (a partial core hit moved a flank); an unrelated core's
                                  known position nearby does not make a position shifted
not_in_assembly:core_absent positions also report the molecules spanning the empty site (empty_sites): a molecule
counts when one of its records has a whole-assembly alignment (read check genome PAF) over the site +/- EMPTY_SPAN
whose local CIGAR over site +/- ANCHOR passes the read check's indel rule (s3 cigar_ok, limit half the core length);
records are merged into molecules with the same FASTQ-header grouping as the search; a molecule with both a core
record at the position and an empty-site record is counted apart (both_forms_molecules), never on either side.
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

VERSION = 'read_site_search_v0.2'  # v0.2: inner-end rule, catalog match by core relation / orientation / coverage, empty-site CIGAR + molecules, shifted by core relation
WORK = u.PROJ / 'core_context' / 'read_site_search' / 'pilot_v0.1'
REPO = Path(__file__).resolve().parent / 'scale_comamonadaceae_v0.1'
ROUND5 = REPO / 'readcheck_round5'
FLANKS = (500, 1000, 2000)  # bp of read sequence taken next to the core on each side
MAIN_FLANK, MAIN_COUNT = 1000, 3
COUNTS = (1, 2, 3, 5)  # independent molecules a merged position needs
HIT_COVERAGE = 0.8  # fraction of the panel sequence a core hit must cover
PLACED_FRACTION = 0.8  # fraction of a flank its best assembly alignment must cover
INNER_CLIP = 20  # bp; the best flank alignment must reach the flank's inner end within this
TOL = rc3.LOCUS_TOL  # bp; placements of one position agree within this
MIN_DELTA_AS = rc3.MIN_DELTA_AS
END_ZONE = 5000  # bp: a flank this close to a contig end (beyond the core length) counts as at that end
EMPTY_SPAN = 1000  # bp each side: a molecule spans an empty site with one genome alignment over site +/- this
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
                       tp=tags.get('tp', 'P'), cg=tags.get('cg', ''))


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
    coordinate of the flank end next to the core (left flank: its last base; right flank: its first base). It is only
    observed when the best alignment reaches that end of the flank (within INNER_CLIP read bases; the flank is in
    read-forward orientation, so the inner end is the query end for L and the query start for R, either strand);
    otherwise inner_end_unaligned, with the position the alignment would give."""
    if not alns:
        return 'unplaced', '', -1, ''
    best = max(alns, key=lambda g: g['score'])
    if (best['qe'] - best['qs']) < PLACED_FRACTION * best['qlen']:
        return 'unplaced', '', -1, ''

    def inner(g):
        return (g['te'] if g['strand'] == '+' else g['ts']) if side == 'L' else (g['ts'] if g['strand'] == '+' else g['te'])
    if (best['qlen'] - best['qe'] if side == 'L' else best['qs']) > INNER_CLIP:
        return 'inner_end_unaligned', best['t'], inner(best), best['strand']
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
    groups = family_groups(a.work)
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
                         right_status=pr[0], right_contig=pr[1], right_pos=pr[2], right_strand=pr[3],
                         left_inner_clip=inner_clip(alns.get(k + '|L', []), 'L'),
                         right_inner_clip=inner_clip(alns.get(k + '|R', []), 'R'))
                if pl[3]:
                    o['core_strand'] = x['strand'] if pl[3] == '+' else FLIP[x['strand']]
                o['status'] = 'placed' if pl[0] == pr[0] == 'placed' else 'unplaceable_context'
            obs.append(o)
    positions = merge(obs, catalog, seqs, member_seq, od, seq_of, presets[0], groups)
    empty_sites(s, positions, runs, recs)
    cols = ['flank', 'structure', 'core_key', 'kind', 'cls', 'raw_cls', 'molecules', 'read_records', 'contig', 'start', 'end',
            'left_contig', 'left_pos', 'left_strand', 'right_contig', 'right_pos', 'right_strand', 'core_strand',
            'core_len', 'catalog', 'catalog_unmatched', 'shifted_from', 'core_in_assembly'] + EMPTY_COLS + [
            'members', 'identity_median', 'runs']
    pc.write_tsv(od / 'positions.tsv.gz', positions, cols)
    ocols = ['flank', 'run', 'read', 'molecule', 'structure', 'core_key', 'kind', 'member', 'qs', 'qe', 'strand',
             'coverage', 'identity', 'status', 'left_status', 'left_contig', 'left_pos', 'left_strand', 'left_inner_clip',
             'right_status', 'right_contig', 'right_pos', 'right_strand', 'right_inner_clip', 'core_strand']
    pc.write_tsv(od / 'observations.tsv.gz', obs, ocols)
    subprocess.run(['gzip', '-f', str(fpaf)], check=True)
    fl.unlink()
    (od / 'search.json').write_text(json.dumps(dict(
        version=VERSION, sample=s, runs=[dict(run=r['run'], read_type=r['read_type'], grouping=r['grouping'],
                                              sha256=r['sha256']) for r in runs],
        panel_hits=len(hits), reads_with_hits=len(hit_reads), molecules_with_hits=len(set(mol.values())),
        hit_coverage=HIT_COVERAGE, flanks=FLANKS, minimap2_panel=MM_PANEL, minimap2_flank=MM_FLANK,
        observation_status={str(L): dict(Counter(o['status'] for o in obs if o['flank'] == L)) for L in FLANKS},
        inner_clip_of_aligned_flanks={str(L): dict(Counter(
            '%s:%s' % (side, 'none' if o[side + '_inner_clip'] == '' else '0' if o[side + '_inner_clip'] == 0 else
                       '1-%d' % INNER_CLIP if o[side + '_inner_clip'] <= INNER_CLIP else '%d-%d' % (INNER_CLIP + 1, TOL)
                       if o[side + '_inner_clip'] <= TOL else '>%d' % TOL)
            for o in obs if o['flank'] == L and o['status'] != 'one_sided' for side in ('left', 'right'))) for L in FLANKS}),
        indent=1) + '\n')


EMPTY_COLS = ['empty_site_molecules', 'empty_site_read_records', 'empty_site_indel_molecules', 'both_forms_molecules']


def empty_site_reads(alns, site_lo, site_hi, core_len):
    """Genome alignments of one read -> 'empty' (one alignment over [site_lo - EMPTY_SPAN, site_hi + EMPTY_SPAN] whose
    CIGAR over [site_lo - ANCHOR, site_hi + ANCHOR] passes s3 cigar_ok with limit max(MIN_INDEL, core_len / 2)),
    'indel' (spans, but every spanning alignment carries an indel event there: possibly the core or another
    insertion), or '' (does not span)."""
    span = [g for g in alns if g['ts'] <= site_lo - EMPTY_SPAN and g['te'] >= site_hi + EMPTY_SPAN]
    if not span:
        return ''
    limit = max(rc3.MIN_INDEL, core_len / 2.0)
    return 'empty' if any(rc3.cigar_ok(g, site_lo - rc3.ANCHOR, site_hi + rc3.ANCHOR, limit) for g in span) else 'indel'


def empty_site_counts(core_reads, site_reads, mol):
    """Molecule counts at one site: core_reads = records with the core there (search), site_reads = {read: 'empty' |
    'indel'}; mol = {read: molecule}. A molecule with a core record and an empty-site record is both_forms (counted on
    neither side); an indel-only molecule is reported apart."""
    m = lambda r: mol.get(r, r)
    core = {m(r) for r in core_reads}
    empty = {m(r) for r, c in site_reads.items() if c == 'empty'}
    indel = {m(r) for r, c in site_reads.items() if c == 'indel'} - empty - core
    return dict(empty_site_molecules=len(empty - core), empty_site_read_records=sum(
        c == 'empty' and m(r) not in core for r, c in site_reads.items()), empty_site_indel_molecules=len(indel),
        both_forms_molecules=len(empty & core))


def empty_sites(sample, positions, runs, hit_recs):
    """Empty-site molecules of every not_in_assembly:core_absent position (see module doc), from the read check's
    whole-assembly PAFs of the sample. Not required for the class."""
    sites = [r for r in positions if r['cls'] == 'not_in_assembly:core_absent']
    if not sites:
        return
    per = defaultdict(lambda: defaultdict(list))  # site index -> read -> spanning-candidate alignments
    for paf in sorted((u.PROJ / 'core_context' / 'readmap' / sample).glob('*.paf.gz')):
        for g in paf_rows(paf):
            for i, r in enumerate(sites):
                lo, hi = min(r['start'], r['end']), max(r['start'], r['end'])
                if g['t'] == r['contig'] and g['ts'] <= lo - EMPTY_SPAN and g['te'] >= hi + EMPTY_SPAN:
                    per[i][g['q']].append(g)
    calls = {i: {q: empty_site_reads(gs, min(sites[i]['start'], sites[i]['end']), max(sites[i]['start'], sites[i]['end']),
                                     sites[i]['core_len']) for q, gs in per[i].items()} for i in range(len(sites))}
    names = {q for c in calls.values() for q in c}
    recs = list(hit_recs)
    for r in runs:  # molecule ids of the empty-site records (FASTQ headers), merged with the core records
        with gzip.open(r['local_file'], 'rt') as h:
            for i, line in enumerate(h):
                if i % 4 == 0:
                    name = line[1:].split(None, 1)[0]
                    if name in names:
                        recs.append((name, rc3.record_ids(line[1:].strip(), r['grouping'])))
    mol = rc3.molecules(recs)
    for i, r in enumerate(sites):
        r.update(empty_site_counts(r['_reads'], calls[i], mol))


def inner_clip(alns, side):
    """Read bases between the best flank alignment and the flank's inner end ('' when unaligned)."""
    if not alns:
        return ''
    best = max(alns, key=lambda g: g['score'])
    return best['qlen'] - best['qe'] if side == 'L' else best['qs']


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


def merge(obs, catalog, seqs, member_seq, od, seq_of, preset, groups):
    """Placed observations -> one row per merged position (per flank length); unplaceable ones -> one row per
    (core, placed side anchor, reason). Molecules are counted once per position."""
    rows = []
    for L in FLANKS:
        by = defaultdict(list)
        for o in obs:
            if o['flank'] == L and o['status'] == 'placed':
                lc, lp, ls, rc, rp, rs = canonical(o)
                o.update(c_left_pos=lp, c_right_pos=rp)
                by[(o['structure'], o['core_key'], lc, ls, rc, rs, o['core_strand'])].append(o)
        placed = []
        for key, os_ in by.items():
            for cl in cluster(os_):
                placed.append(position_row(L, key, cl, catalog, seqs, member_seq, od, seq_of, preset, groups))
        rows += mark_shifted(placed, groups)
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


FLIP = {'+': '-', '-': '+'}


def core_families(ck, groups=None):
    """Protein families of a core key with their sign (orientation relative to the core); no sign means '+'. A
    negative-control key (NEG_<product>) has the families of the control genes' proteins (groups['neg:' + key])."""
    if ck.startswith('NEG_'):
        return {f: '+' for f in (groups or {}).get('neg:' + ck, ())}
    return {f.rstrip('+-'): '-' if f.endswith('-') else '+' for f in ck.split(',') if f.startswith('PF_')}


def core_relation(a, b, groups):
    """('same', '+'), ('homolog', sign) or ('unrelated', ''). Homolog: the cores share a protein family, or a family
    of each lies in a common homology group (groups: family -> homology groups). sign = '+' when the shared families
    have the same orientation in both cores, '-' when opposite, '' when undetermined (the pairs disagree)."""
    if a == b:
        return 'same', '+'
    fa, fb = core_families(a, groups), core_families(b, groups)
    signs = {'+' if fa[f] == fb[f] else '-' for f in set(fa) & set(fb)}
    if not signs:
        signs = {'+' if fa[x] == fb[y] else '-' for x in fa for y in fb if groups.get(x, set()) & groups.get(y, set())}
    if not signs:
        return 'unrelated', ''
    return 'homolog', signs.pop() if len(signs) == 1 else ''


def family_groups(work):
    """family -> homology groups for the families of the panel and catalog core keys, and 'neg:' + control key ->
    protein families of the control genes (their proteins in the family catalog); cached in panel/."""
    path = work / 'panel' / 'family_groups.tsv'
    if not path.exists():
        prot = {}  # control gene protein -> control key
        for m in pc.read_tsv(work / 'panel' / 'members.tsv'):
            if m['member'].startswith('neg|'):
                s, contig, lo = m['first_occurrence'].split('|')
                g = next(g for g in pc.read_tsv(u.PROJ / 'core_context' / 'scale_comamonadaceae_v0.1' / 'genes' /
                                                ('%s.tsv.gz' % s))
                         if g['type'] == 'CDS' and g['contig'] == contig and g['start0'] == lo)
                prot[g['protein_sha256']] = 'neg:' + m['core_key']
        want = {f for t in ('members.tsv', 'catalog.tsv') for m in pc.read_tsv(work / 'panel' / t)
                for f in core_families(m['core_key'])}
        got, neg = defaultdict(set), set()
        with gzip.open(str(u.PROJ / 'core_context' / 'scale_comamonadaceae_v0.1' / 'families' / 'catalog.tsv.gz'),
                       'rt') as h:
            next(h)
            for line in h:
                f = line.split('\t', 4)
                got[f[1]].add(f[3])
                if f[0] in prot:
                    neg.add((prot[f[0]], f[1]))
        want |= {f for _, f in neg}
        rows = [dict(family=f, homology_group=g) for f in sorted(want) for g in sorted(got[f])]
        rows += [dict(family=k, homology_group=f) for k, f in sorted(neg)]
        pc.write_tsv(path, rows, ['family', 'homology_group'])
    out = defaultdict(set)
    for r in pc.read_tsv(path):
        out[r['family']].add(r['homology_group'])
    return out


def match_catalog(contig, x, y, strand, ck, catalog, groups):
    """Catalog occurrences at the interval [x, y) of contig for a core ck lying on strand there (see module doc):
    (class or None, matched occurrence ids, unmatched notes). Only occurrences overlapping the interval are looked at;
    an unrelated core's occurrence is only noted."""
    lo, hi = min(x, y), max(x, y)
    match, miss, notes = defaultdict(list), [], []
    for c in catalog:
        s0, e0 = int(c['start0']), int(c['end0'])
        ov = min(hi, e0) - max(lo, s0)
        if c['contig'] != contig or ov <= 0:
            continue
        rel, sign = core_relation(ck, c['core_key'], groups)
        if rel == 'unrelated':
            notes.append('%s:unrelated' % c['occurrence_id'])
            continue
        why = (['orientation_undetermined'] if not sign else
               ['orientation'] if c['strand'] != (strand if sign == '+' else FLIP[strand]) else [])
        if not (ov >= HIT_COVERAGE * (e0 - s0) and ov >= HIT_COVERAGE * (hi - lo)):
            why.append('coverage')
        if why:
            miss.append((rel != 'same', -ov, c['occurrence_id'], why))
            notes.append('%s:%s:%s' % (c['occurrence_id'], rel, '+'.join(why)))
        else:
            match[rel].append(c['occurrence_id'])
    if match['same']:
        return 'known_catalog_position', ';'.join(sorted(match['same'])), ';'.join(notes)
    if match['homolog']:
        return 'known_catalog_position_other_core', ';'.join(sorted(match['homolog'])), ';'.join(notes)
    if miss:
        first = min(miss)
        return 'catalog_mismatch:' + first[3][0], first[2], ';'.join(notes)
    return None, '', ';'.join(notes)


def position_row(L, key, cl, catalog, seqs, member_seq, od, seq_of, preset, groups):
    st, ck, lc, ls, rc, rs, cs = key
    lp = int(median([o['c_left_pos'] for o in cl]))
    rp = int(median([o['c_right_pos'] for o in cl]))
    core_len = int(median([o['qe'] - o['qs'] for o in cl]))
    geo, x, y = geometry(lc, lp, ls, rc, rp, rs, core_len)
    cat, inside, notes = '', '', ''
    if geo == 'junction':
        n = len(seqs[lc])
        cls = ('contig_end_join' if lc == rc and ls == rs and lp > n - core_len - END_ZONE and rp < core_len + END_ZONE
               else 'not_in_assembly:junction')
    else:
        cls, cat, notes = match_catalog(lc, x, y, cs, ck, catalog, groups)
        if cls is None:
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
                left_strand=ls, right_contig=rc, right_pos=rp, right_strand=rs, core_strand=cs, core_len=core_len,
                catalog=cat, catalog_unmatched=notes, core_in_assembly=inside,
                members=len({o['member'] for o in cl}), identity_median=median([o['identity'] for o in cl]),
                runs=';'.join(sorted({o['run'] for o in cl})), _reads={o['read'] for o in cl})


def mark_shifted(rows, groups):
    """A non-known position whose assembly interval (or, for a junction, either flank anchor) overlaps a known
    position of the same sample is the known position seen with shifted anchors (a core hit that covers only part
    of the gene moves a flank) only if the known position's core is the same or a homologous one (core_relation)
    and, when the relation gives a sign, in the consistent orientation: cls shifted_known_position, raw_cls keeps
    the class it had, shifted_from names the known position. Found with the negative-control genes of the pilot;
    v0.1 did not check the core. Two known positions of one core matched to the same catalog occurrence are one
    position: the one with fewer molecules (ties: later start) is shifted_known_position as well."""
    known = defaultdict(list)
    first = {}
    for r in sorted((r for r in rows if r['cls'].startswith('known_catalog_position')),
                    key=lambda r: (-r['molecules'], r['start'])):
        k = (r['core_key'], r['catalog'])
        if k in first:
            r['raw_cls'], r['cls'] = r['cls'], 'shifted_known_position'
            r['shifted_from'] = '%s:%s-%s:%s' % (r['contig'], first[k]['start'], first[k]['end'], r['core_key'])
        else:
            first[k] = r
            known[r['contig']].append(r)
    for r in rows:
        if r['cls'].startswith(('known_catalog_position', 'shifted')):
            continue
        near = ([k for k in known[r['contig']] if k['start'] < r['end'] + TOL and k['end'] > r['start'] - TOL]
                if r['contig'] else
                [k for c, p in ((r['left_contig'], r['left_pos']), (r['right_contig'], r['right_pos']))
                 for k in known[c] if k['start'] - TOL <= p <= k['end'] + TOL])
        for k in near:
            rel, sign = core_relation(r['core_key'], k['core_key'], groups)
            if rel != 'unrelated' and (not sign or r['core_strand'] == (k['core_strand'] if sign == '+'
                                                                         else FLIP[k['core_strand']])):
                r['raw_cls'], r['cls'] = r['cls'], 'shifted_known_position'
                r['shifted_from'] = '%s:%s-%s:%s' % (k['contig'], k['start'], k['end'], k['core_key'])
                break
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


def annotate_omissions(rows):
    """Positions where the assembly between the flanks holds the read's core sequence but the catalog lists no
    reviewed occurrence, split by the annotated CDS there:
       catalog_omission:same_family            a CDS of the core's own protein family: a gap of the catalog
       catalog_omission:no_annotated_cds       no CDS lies within the interval or covers >= 80% of it
       homolog_not_in_catalog:same_homology_group   a CDS of another family of the core's homology group, a family
                                               the reviewed catalog does not cover
       homolog_not_in_catalog:other_family     a CDS of an unrelated family (a DNA-level homolog of the core)"""
    om = [r for r in rows if r['cls'] == 'catalog_omission']
    if not om:
        return
    gdir = u.PROJ / 'core_context' / 'scale_comamonadaceae_v0.1' / 'genes'
    genes = {}
    for s in {r['sample'] for r in om}:
        d = defaultdict(list)
        for g in pc.read_tsv(gdir / ('%s.tsv.gz' % s)):
            if g['type'] == 'CDS':
                d[g['contig']].append(g)
        genes[s] = d
    within = {}
    for r in om:
        x, y = int(r['start']), int(r['end'])
        within[id(r)] = [g for g in genes[r['sample']][r['contig']]  # within the interval, or covering >= 80% of it
                         if int(g['start0']) >= x - TOL * 4 and int(g['end0']) <= y + TOL * 4
                         or min(y, int(g['end0'])) - max(x, int(g['start0'])) >= HIT_COVERAGE * (y - x)]
    want = {g['protein_sha256'] for gs in within.values() for g in gs}
    core_fams = {f.rstrip('+-') for r in om for f in r['core_key'].split(',')}
    fam, group_of = {}, defaultdict(set)  # protein -> (family, homology group); core family -> homology groups
    with gzip.open(str(u.PROJ / 'core_context' / 'scale_comamonadaceae_v0.1' / 'families' / 'catalog.tsv.gz'), 'rt') as h:
        next(h)
        for line in h:
            f = line.split('\t', 4)
            if f[0] in want:
                fam[f[0]] = (f[1], f[3])
            if f[1] in core_fams:
                group_of[f[1]].add(f[3])
    core_group = {r['core_key']: set().union(*(group_of[f.rstrip('+-')] for f in r['core_key'].split(',')))
                  for r in om}
    for r in om:
        gs = within[id(r)]
        fams = [fam.get(g['protein_sha256'], ('', '')) for g in gs]
        r['annotation'] = ';'.join('%s:%s:%s' % (g['product'], f[0], f[1]) for g, f in zip(gs, fams))
        if not gs:
            r['cls'] = 'catalog_omission:no_annotated_cds'
        elif any(f[0] and f[0] in r['core_key'] for f in fams):
            r['cls'] = 'catalog_omission:same_family'
        elif any(f[1] and f[1] in core_group[r['core_key']] for f in fams):
            r['cls'] = 'homolog_not_in_catalog:same_homology_group'
        else:
            r['cls'] = 'homolog_not_in_catalog:other_family'


def cmd_combine(a):
    samples = pilot_samples(a.out)
    rows, runs = [], {}
    for s in samples:
        od = a.work / 'samples' / s
        for r in pc.read_tsv(od / 'positions.tsv.gz'):
            rows.append(dict(sample=s, **r))
        runs[s] = json.loads((od / 'search.json').read_text())
    annotate_omissions(rows)
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
    cols = ['sample'] + [c for c in rows[0] if c != 'sample'] + ['annotation']
    pc.write_tsv(out / 'positions.tsv.gz', rows, cols)
    main = [r for r in rows if int(r['flank']) == MAIN_FLANK]
    for cls, name in (('catalog_omission', 'catalog_omissions'), ('not_in_assembly', 'not_in_assembly'),
                      ('catalog_mismatch', 'catalog_mismatches'), ('shifted_known_position', 'shifted_known_positions'),
                      ('unplaceable_context', 'unplaceable_contexts'),
                      ('homolog_not_in_catalog', 'homologs_not_in_catalog'), ('contig_end_join', 'contig_end_joins')):
        sel = [r for r in main if r['cls'].startswith(cls)]
        if sel:
            pc.write_tsv(out / ('%s.tsv' % name), sorted(sel, key=lambda r: (-int(r['molecules']), r['sample'])),
                         cols)
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
