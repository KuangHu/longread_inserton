"""Batch layer of the read-site search scale-up (read_site_search.py v0.3.2 rules unchanged):
physical sites across cores, batch acceptance, and cleanup lists.

sites    Merged positions of the main setting (MAIN_FLANK) of one batch are grouped into physical sites per sample:
         positions of different (homologous) core seeds whose left junctions and right junctions both agree within TOL
         (either orientation) are one site, and a shifted_known_position joins the known position it was shifted
         from. Several seeds hitting one place are never counted as several sites. Each site gets one category,
         from its positions' classes in this order:
           known                       a known_catalog_position (same or homologous core) is there
           in_assembly_not_in_catalog  catalog_omission:* / homolog_not_in_catalog:*: the assembly holds the core there,
                                       the catalog does not list it
           not_in_assembly_candidate   not_in_assembly:core_absent / other_sequence: reads carry the core at a junction
                                       the assembly shows without it (candidate structure in part of the DNA)
           structural_anomaly          not_in_assembly:junction (far apart / opposite strands / other contig) and
                                       contig_end_join
           unresolved                  catalog_mismatch:* only (an occurrence overlaps but fails orientation/coverage)
         and flags: below_molecule_threshold (no position of the site has >= MAIN_COUNT molecules; the site is then
         reported apart, not counted), boundary_unknown (a position with boundary unknown), size_ge_2500, and
         cores_disagree (positions of the site fall in more than one category). Unplaceable contexts (a flank
         ambiguous = repeated placement, unplaced, no anchor, inner end unaligned) and one_sided observations are not
         sites; they are counted per reason.
accept   Completeness by id of one batch: every sample of samples.tsv has search.json with the frozen version, the
         runs and read checksums of samples.tsv, and positions.tsv.gz; runs not in the ready manifest make the sample
         input_not_ready. Then the combined tables, the candidate review outputs, resource use (sacct) and disk.
         Writes <out>/acceptance.json; exit 1 when not complete.
cleanup  After acceptance and archiving only. Lists rebuildable files of the batch scratch directory (intermediate
         alignments and temporary FASTA/BAM of the review); files any listed consumer still needs are kept. Writes
         <out>/cleanup_manifest.tsv; deletes only with --apply, and only paths under --work that are listed.
Usage: python scale_sites.py sites|accept|cleanup --out <batch repo dir> --work <batch scratch dir> [...]
"""
import argparse
import hashlib
import json
import os
import subprocess
import sys
import tarfile
from collections import Counter, defaultdict
from pathlib import Path

import protein_context as pc
import read_site_search as rs

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'scripts'))
import reads_ready  # noqa: E402

CATEGORY = (('known', ('known_catalog_position',)),
            ('in_assembly_not_in_catalog', ('catalog_omission', 'homolog_not_in_catalog')),
            ('not_in_assembly_candidate', ('not_in_assembly:core_absent', 'not_in_assembly:other_sequence')),
            ('structural_anomaly', ('not_in_assembly:junction', 'contig_end_join')),
            ('unresolved', ('catalog_mismatch',)))


def category(cls):
    for name, prefixes in CATEGORY:
        if cls.startswith(prefixes):
            return name
    return ''


def ends(r):
    return ((r['left_contig'], int(r['left_pos'])) if r['left_pos'] not in ('', None) else None,
            (r['right_contig'], int(r['right_pos'])) if r['right_pos'] not in ('', None) else None)


def near(a, b):
    return a and b and a[0] == b[0] and abs(a[1] - b[1]) <= rs.TOL


def same_site(x, y):
    (xl, xr), (yl, yr) = ends(x), ends(y)
    return near(xl, yl) and near(xr, yr) or near(xl, yr) and near(xr, yl)


def group_sites(rows):
    """Union-find of a sample's site positions (rows: classed positions and shifted ones of the main setting)."""
    parent = list(range(len(rows)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def union(i, j):
        parent[find(i)] = find(j)

    # compare only rows with an end in the same or a neighbouring TOL bin of a contig
    keyed = defaultdict(list)
    for i, r in enumerate(rows):
        for e in ends(r):
            if e:
                keyed[(e[0], e[1] // rs.TOL)].append(i)
    for i, r in enumerate(rows):
        cand = set()
        for e in ends(r):
            if e:
                for b in (e[1] // rs.TOL - 1, e[1] // rs.TOL, e[1] // rs.TOL + 1):
                    cand.update(keyed.get((e[0], b), ()))
        for j in cand:
            if j != i and same_site(r, rows[j]):
                union(i, j)
    known_at = {}
    for i, r in enumerate(rows):
        if r['cls'].startswith('known_catalog_position'):
            known_at['%s:%s-%s:%s' % (r['contig'], r['start'], r['end'], r['core_key'])] = i
    for i, r in enumerate(rows):
        if r['cls'] == 'shifted_known_position' and r['shifted_from'] in known_at:
            union(i, known_at[r['shifted_from']])
    groups = defaultdict(list)
    for i in range(len(rows)):
        groups[find(i)].append(rows[i])
    return list(groups.values())


def site_row(sid, s, ps):
    cats = Counter(category(r['cls']) for r in ps if category(r['cls']))
    cat = next((c for c, _ in CATEGORY if cats[c]), 'shifted_only')
    mols = [int(r['molecules']) for r in ps]
    ins = sorted(int(r['insert_len']) for r in ps if r['insert_len'] not in ('', None))
    flags = [f for f, on in (
        ('below_molecule_threshold', max(mols) < rs.MAIN_COUNT),
        ('boundary_unknown', any(r['boundary'] == 'unknown' for r in ps if category(r['cls']) == cat)),
        ('size_ge_2500', any(r['insert_size'].startswith('ge_') for r in ps if category(r['cls']) == cat)),
        ('cores_disagree', len(cats) > 1)) if on]
    lead = max((r for r in ps if category(r['cls']) == cat), key=lambda r: int(r['molecules']), default=ps[0])
    return dict(site=sid, sample=s, category=cat, flags=';'.join(flags), max_molecules=max(mols),
                positions=len(ps), cores=len({r['core_key'] for r in ps}),
                core_keys=';'.join(sorted({r['core_key'] for r in ps})),
                classes=';'.join('%s:%d' % kv for kv in sorted(Counter(r['cls'] for r in ps).items())),
                left_contig=lead['left_contig'], left_pos=lead['left_pos'], right_contig=lead['right_contig'],
                right_pos=lead['right_pos'], lead_core=lead['core_key'], lead_class=lead['cls'],
                lead_boundary=lead['boundary'], insert_len_range='%d-%d' % (ins[0], ins[-1]) if ins else '',
                catalog=';'.join(sorted({o for r in ps for o in r['catalog'].split(';') if o})))


def neg_sites(rows):
    """Negative-control genes (single-copy by annotation): physical sites (>= MAIN_COUNT molecules) per gene, after
    the same merge as for cores, so one locus split into two position clusters counts once; more than one site
    is listed: one_locus_junction_differs_gt_TOL (the sites' intervals overlap: one locus whose clusters disagree on
    a junction by more than TOL) or separate_loci (another copy elsewhere)."""
    per, extra = Counter(), []
    neg = [r for r in rows if r['kind'] == 'neg' and int(r['molecules']) >= rs.MAIN_COUNT
           and (category(r['cls']) or r['cls'] == 'shifted_known_position')]
    for key in sorted({(r['sample'], r['core_key']) for r in neg}):
        g = group_sites([r for r in neg if (r['sample'], r['core_key']) == key])
        g = [x for x in g if any(category(r['cls']) for r in x)]
        per[len(g)] += 1
        if len(g) > 1:
            iv = [(x[0]['left_contig'], min(int(x[0]['left_pos']), int(x[0]['right_pos'])),
                   max(int(x[0]['left_pos']), int(x[0]['right_pos']))) for x in g]
            same = all(c == iv[0][0] and lo < iv[0][2] and hi > iv[0][1] for c, lo, hi in iv)
            extra.append(dict(sample=key[0], gene=key[1], kind='one_locus_junction_differs_gt_TOL' if same
                              else 'separate_loci', sites=['%s:%d-%d' % v for v in iv]))
    return dict(genes_by_site_count={str(k): v for k, v in sorted(per.items())}, genes_with_several_sites=extra)


def cmd_sites(a):
    rows = [r for r in pc.read_tsv(a.out / 'positions.tsv.gz') if int(r['flank']) == rs.MAIN_FLANK]
    core = [r for r in rows if r['kind'] == 'core']
    sites, links = [], []
    for s in sorted({r['sample'] for r in core}):
        ps = [r for r in core if r['sample'] == s and (category(r['cls']) or r['cls'] == 'shifted_known_position')]
        for k, g in enumerate(sorted(group_sites(ps), key=lambda g: min(ends(r)[0] or ('', 0) for r in g))):
            sid = '%s:S%04d' % (s, k)
            sites.append(site_row(sid, s, g))
            links += [dict(site=sid, sample=s, core_key=r['core_key'], cls=r['cls'], molecules=r['molecules'],
                           left_contig=r['left_contig'], left_pos=r['left_pos'], right_contig=r['right_contig'],
                           right_pos=r['right_pos'], boundary=r['boundary'], insert_len=r['insert_len'],
                           insert_size=r['insert_size']) for r in g]
    cols = list(sites[0]) if sites else ['site']
    pc.write_tsv(a.out / 'sites.tsv.gz', sites, cols)
    pc.write_tsv(a.out / 'site_positions.tsv.gz', links, list(links[0]) if links else ['site'])
    counted = [x for x in sites if 'below_molecule_threshold' not in x['flags']]
    new = [x for x in counted if x['category'] == 'not_in_assembly_candidate']
    new_pos = [r for r in links if any(r['site'] == x['site'] for x in new)
               and category(r['cls']) == 'not_in_assembly_candidate']
    unpl = Counter()
    for r in core:
        if r['cls'].startswith('unplaceable') and int(r['molecules']) >= rs.MAIN_COUNT:
            why = r['cls'].split(':', 1)[1]
            unpl['ambiguous_placement (repeat)' if 'ambiguous' in why else
                 'inner_end_unaligned (alignment truncated)' if 'inner_end_unaligned' in why else
                 'no_anchor_or_unplaced'] += 1
    summ = dict(
        version=rs.VERSION, main=dict(flank=rs.MAIN_FLANK, molecules=rs.MAIN_COUNT),
        sites_by_category={c: sum(x['category'] == c for x in counted) for c, _ in CATEGORY},
        sites_below_molecule_threshold=dict(Counter(x['category'] for x in sites
                                                    if 'below_molecule_threshold' in x['flags'])),
        flags_of_counted_sites=dict(Counter(f for x in counted for f in x['flags'].split(';') if f)),
        sites_with_several_cores=sum(x['cores'] > 1 for x in counted),
        positions_merged_into_sites=sum(x['positions'] for x in counted),
        not_in_assembly_candidates=dict(sites=len(new), positions=len(new_pos),
                                        cores=len({r['core_key'] for r in new_pos}),
                                        samples=len({x['sample'] for x in new}),
                                        by_class=dict(Counter(r['cls'] for r in new_pos))),
        negative_control_sites_per_gene=neg_sites(rows),
        unplaceable_contexts_ge_N=dict(unpl),
        one_sided_observation_groups=sum(r['cls'] == 'one_sided' for r in core))
    (a.out / 'sites.json').write_text(json.dumps(summ, indent=1) + '\n')
    print(json.dumps(summ, indent=1))


def sacct(jobs):
    if not jobs:
        return []
    out = subprocess.run(['sacct', '-n', '-P', '-j', ','.join(jobs), '--format=JobID,JobName,State,Elapsed,MaxRSS,'
                          'NodeList,ExitCode'], stdout=subprocess.PIPE, universal_newlines=True).stdout
    rows = []
    for line in out.splitlines():
        f = line.split('|')
        if f[0].endswith(('.batch',)) or '.' not in f[0]:
            rows.append(dict(zip(('job', 'name', 'state', 'elapsed', 'max_rss', 'node', 'exit'), f)))
    return rows


def du(path):
    return int(subprocess.run(['du', '-sb', str(path)], stdout=subprocess.PIPE, universal_newlines=True)
               .stdout.split()[0]) if Path(path).exists() else 0


def cmd_accept(a):
    samples = list(pc.read_tsv(a.out / 'samples.tsv'))
    status, problems = {}, []
    for r in samples:
        s, od = r['sample'], a.work / 'samples' / r['sample']
        runs = r['runs'].split(';')
        if reads_ready.not_ready(runs):
            status[s] = 'input_not_ready'
            continue
        j = od / 'search.json'
        if not j.exists() or not (od / 'positions.tsv.gz').exists():
            status[s] = 'failed_or_missing'
            continue
        sj = json.loads(j.read_text())
        bad = [w for w, b in (('version', sj['version'] != rs.VERSION),
                              ('runs', sorted(x['run'] for x in sj['runs']) != sorted(runs)),
                              ('sha256', sorted(x['sha256'] for x in sj['runs']) != sorted(r['sha256'].split(';'))))
               if b]
        status[s] = 'complete' if not bad else 'mismatch:' + ','.join(bad)
    n = Counter(v.split(':')[0] for v in status.values())
    summ = json.loads((a.out / 'summary.json').read_text()) if (a.out / 'summary.json').exists() else {}
    if summ and sorted(summ['samples']) != sorted(status):
        problems.append('combined samples differ from samples.tsv')
    if summ and summ['version'] != rs.VERSION:
        problems.append('combined version %s' % summ['version'])
    sites = json.loads((a.out / 'sites.json').read_text()) if (a.out / 'sites.json').exists() else {}
    if not sites:
        problems.append('sites.json missing')
    cand = a.out / 'candidates'
    ev = list(pc.read_tsv(cand / 'evidence_summary.tsv')) if (cand / 'evidence_summary.tsv').exists() else None
    st = list(pc.read_tsv(cand / 'structure_summary.tsv')) if (cand / 'structure_summary.tsv').exists() else None
    if ev is None:
        problems.append('candidate evidence not run')
    elif any(c['status'] == 'tested' for c in ev) and st is None:
        problems.append('structure check not run')
    jobs = [x for x in a.jobs.split(',') if x]
    res = sacct(jobs)
    failed_jobs = [x for x in res if x['state'] not in ('COMPLETED',) and '.' not in x['job']]
    main = summ.get('grid', {}).get('L%d_N%d' % (rs.MAIN_FLANK, rs.MAIN_COUNT), {})
    out = dict(
        version=rs.VERSION, batch=str(a.out.name), expected_samples=len(samples), sample_status=dict(n),
        samples_not_complete={s: v for s, v in status.items() if v != 'complete'},
        complete=n['complete'] == len(samples) and not problems and not failed_jobs, problems=problems,
        known_recovery=main.get('known_recovery'), negative_control=dict(
            (k, v) for k, v in main.get('negative_control', {}).items() if k != 'extra_positions'),
        negative_control_extra_positions=len(main.get('negative_control', {}).get('extra_positions', [])),
        sites=sites, candidates_reviewed=dict(
            evidence=dict(Counter(c['status'] for c in ev)) if ev else {},
            structure=dict(Counter((c.get('status'), c.get('target_1_2_genes')) and
                                   '%s|%s' % (c.get('status'), c.get('target_1_2_genes', '')) for c in st)) if st else {}),
        jobs=res, failed_jobs=failed_jobs,
        max_rss_gb=round(max([float(x['max_rss'].rstrip('K')) / 1e6 for x in res if x['max_rss'].endswith('K')]
                             or [0]), 2),
        disk_bytes=dict(work=du(a.work), out=du(a.out)))
    (a.out / 'acceptance.json').write_text(json.dumps(out, indent=1) + '\n')
    print(json.dumps({k: v for k, v in out.items() if k != 'jobs'}, indent=1))
    sys.exit(0 if out['complete'] else 1)


def sha256(path):
    h = hashlib.sha256()
    with open(str(path), 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


# rebuildable files in a batch scratch directory: intermediate alignments, temporary FASTA/BAM of the review
REBUILDABLE = (('samples/*/segments.paf.gz', 'boundary-stage alignment (rebuilt by search)'),
               ('samples/*/flanks.paf.gz', 'placement-stage alignment (rebuilt by search)'),
               ('samples/*/panel_*.paf.gz', 'reads-to-panel alignment (rebuilt by search)'),
               ('structure/*/round*.bam', 'consensus round alignment'),
               ('structure/*/round*.bam.bai', 'consensus round alignment index'),
               ('structure/*/draft*.fa', 'consensus draft'),
               ('structure/*/reads.fa', 'reads extracted from the verified FASTQ'),
               ('structure/*/*_q.fa', 'blast query'), ('structure/*/*_s.fa', 'blast subject'),
               ('structure/proteins.faa', 'blastp subject (rebuilt from the protein pool)'))
KEEP_NOTE = ('kept: samples/*/search.json, positions.tsv.gz, observations.tsv.gz (per-molecule evidence); panel/; '
             'candidates/ (allele windows and read-to-allele / genome alignments of the candidates); '
             'structure/*/final.fa, final.bam(.bai) (consensus and its read support)')


def cmd_cleanup(a):
    acc = json.loads((a.out / 'acceptance.json').read_text())
    if not acc.get('complete'):
        sys.exit('batch not accepted; no cleanup')
    arch = a.work / 'archive.tar'
    if not (a.out / 'archive.json').exists():
        sys.exit('no archive recorded; run with --archive first')
    rec = json.loads((a.out / 'archive.json').read_text())
    if not arch.exists() or sha256(arch) != rec['sha256']:
        sys.exit('archive missing or checksum differs: %s' % arch)
    work = a.work.resolve()
    keep_samples = set(filter(None, a.consumers.split(',')))  # samples a later stage (e.g. positive test) reuses
    rows = []
    for pat, why in REBUILDABLE:
        for f in sorted(work.glob(pat)):
            f = f.resolve()
            assert str(f).startswith(str(work) + os.sep), f
            s = f.relative_to(work).parts[1] if f.relative_to(work).parts[0] == 'samples' else ''
            if s in keep_samples and pat.startswith('samples/*/panel_'):
                rows.append(dict(path=str(f), bytes=f.stat().st_size, reason=why, action='keep_consumer'))
                continue
            rows.append(dict(path=str(f), bytes=f.stat().st_size, reason=why, action='delete'))
    freed = 0
    if a.apply:
        for r in rows:
            if r['action'] == 'delete':
                Path(r['path']).unlink()
                freed += int(r['bytes'])
                r['action'] = 'deleted'
    pc.write_tsv(a.out / 'cleanup_manifest.tsv', rows, ['path', 'bytes', 'reason', 'action'])
    out = dict(work=str(work), listed=len(rows), bytes_listed=sum(int(r['bytes']) for r in rows if
                                                                   r['action'] != 'keep_consumer'),
               bytes_freed=freed, applied=a.apply, kept_for_consumers=sorted(keep_samples), kept=KEEP_NOTE,
               work_bytes_after=du(work))
    (a.out / 'cleanup.json').write_text(json.dumps(out, indent=1) + '\n')
    print(json.dumps(out, indent=1))


def cmd_archive(a):
    """tar of the batch scratch results that are kept (no reads), with its sha256 recorded in the repo."""
    work = a.work.resolve()
    arch = work / 'archive.tar'
    rebuild = {f.resolve() for pat, _ in REBUILDABLE for f in work.glob(pat)}
    files = sorted(f for f in work.rglob('*') if f.is_file() and f.resolve() not in rebuild and f != arch)
    with tarfile.open(str(arch), 'w') as t:
        for f in files:
            t.add(str(f), arcname=str(f.relative_to(work)))
    rec = dict(path=str(arch), sha256=sha256(arch), bytes=arch.stat().st_size, files=len(files))
    (a.out / 'archive.json').write_text(json.dumps(rec, indent=1) + '\n')
    print(json.dumps(rec, indent=1))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('mode', choices=['sites', 'accept', 'archive', 'cleanup'])
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--work', type=Path, required=True)
    p.add_argument('--jobs', default='', help='accept: comma-separated SLURM job ids of the batch')
    p.add_argument('--consumers', default='', help='cleanup: samples whose rebuildable files a later stage reuses')
    p.add_argument('--apply', action='store_true')
    a = p.parse_args()
    rs.set_work(a.work)
    dict(sites=cmd_sites, accept=cmd_accept, archive=cmd_archive, cleanup=cmd_cleanup)[a.mode](a)


if __name__ == '__main__':
    main()
