"""Add existing long-read evidence to the strict U-X-D pairs as separate columns. Pass / fail is not changed.

Two independent sources, never merged into one label:

1. Interval-level read check, round 5 v0.2 (readcheck_round5_v02, ucd_readcheck v0.2). Per side, a result is taken
   only when the tested interval is the same as this row's X:
   - same sample (same manifest assembly), same contig, same lo / hi (0-based half-open, forward strand);
   - and the round-5 record of this occurrence has the same strand.
   If the same occurrence was tested on another interval, nothing is inherited: the side is tested_other_interval,
   with those intervals listed. The raw interval class, verdict, count mode and counts are copied as they are.
     rc_status  supported | locus_ambiguous | insufficient_coverage | conflict   (round-5 interval class, unchanged)
                tested_other_interval   occurrence tested, but on a different interval
                not_tested              the assembly has usable long-read runs, this interval was not tested
                no_reads                no usable long-read run listed for the assembly's BioSample
2. Read-site search (read_site_search v0.3.2), batches validation1-2 and scale1-3 (133 real samples) and the pilot
   (development samples, labelled as such). Rule as in its own recovery count: main flank 1,000 bp, >= 3 molecules,
   a known_catalog_position* position that lists this occurrence.
     rss_status recovered | seen_below_threshold | not_recovered | not_in_search_panel | sample_not_searched
   not_recovered is not 'unsupported': the search may miss positions (recovery of known positions was ~95%).
"""
import argparse
import gzip
import json
import time
from collections import Counter, defaultdict
from pathlib import Path

import ucd_batch as ub
import ucd_compare as uc

REPO = uc.REPO
R5 = REPO / 'readcheck_round5_v02'
RSS_REPO = REPO
RSS_WORK = Path('/global/scratch/users/kh36969/longread_callout/core_context/read_site_search')
RSS_BATCHES = (('pilot', 'pilot_v0.1', 'read_site_search_pilot', 'development'),
               ('validation1', 'validation1', 'read_site_search_validation1', 'real'),
               ('validation2', 'validation2', 'read_site_search_validation2', 'real'),
               ('scale1', 'scale1', 'read_site_search_scale1', 'real'),
               ('scale2', 'scale2', 'read_site_search_scale2', 'real'),
               ('scale3', 'scale3', 'read_site_search_scale3', 'real'))
POOL = REPO.parent.parent / 'fna_pool' / 'scaleup_stage1_20261026' / 'comamonadaceae' / 'pool.tsv'
MAIN_FLANK, MAIN_COUNT = 1000, 3
EVID = (['rc_%s_%s' % (s, f) for s in 'ab' for f in ('status', 'test_id', 'verdict', 'count_mode', 'P_matched',
                                                       'A_matched', 'locus_ambiguous', 'other_intervals')] +
        ['rc_pair'] +
        ['rss_%s_%s' % (s, f) for s in 'ab' for f in ('status', 'batch', 'molecules')] + ['rss_pair'])
KEEP = ['pair_id', 'interval_kind', 'structure', 'core_key', 'pair_scope', 'a_occurrence', 'a_sample', 'a_assembly',
        'a_contig', 'a_strand', 'a_x_interval', 'b_occurrence', 'b_sample', 'b_assembly', 'b_contig', 'b_strand',
        'b_x_interval', 'x_len_a', 'x_len_b', 'chk_insertion_identity', 'final_status', 'code_version']


def load_r5():
    pos = {}
    for r in uc.read_tsv(R5 / 'positions.tsv.gz'):
        pos[(r['sample'], r['contig'], r['lo'], r['hi'])] = r
    by_occ = defaultdict(list)
    for r in uc.read_tsv(R5 / 'interval_map.tsv.gz'):
        by_occ[r['occurrence']].append(r['test_id'])
    tid = {r['test_id']: r for r in pos.values()}
    strand = {}
    for r in uc.read_tsv(R5 / 'pairs.tsv.gz'):
        for s in 'ab':
            strand[r[s + '_occurrence']] = r[s + '_strand']
    samples = set(Path(R5 / 'samples.txt').read_text().split())
    return pos, by_occ, tid, strand, samples


def load_rss():
    searched, panel, found = {}, defaultdict(set), defaultdict(list)
    for name, wdir, rdir, group in RSS_BATCHES:
        rd = RSS_REPO / rdir
        if not (rd / 'samples.tsv').exists():
            continue
        for r in uc.read_tsv(rd / 'samples.tsv'):
            searched[r['sample']] = '%s (%s)' % (name, group)
        for c in uc.read_tsv(RSS_WORK / wdir / 'panel' / 'catalog.tsv'):
            panel[c['sample']].add(c['occurrence_id'])
        for r in uc.read_tsv(rd / 'positions.tsv.gz'):
            if int(r['flank']) != MAIN_FLANK or not r['cls'].startswith('known_catalog_position'):
                continue
            for o in r['catalog'].split(';'):
                if o:
                    found[(r['sample'], o)].append(int(r['molecules']))
    return searched, panel, found


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--rows', required=True, help='strict (or pass) table of ucd_combine')
    p.add_argument('--out-table', required=True, help='full rows + evidence (scratch)')
    p.add_argument('--out-compact', required=True, help='key columns + evidence (git)')
    p.add_argument('--summary', required=True)
    a = p.parse_args()
    pos, by_occ, tid, r5strand, r5samples = load_r5()
    searched, panel, found = load_rss()
    usable = {r['biosample']: int(r['n_usable_runs'] or 0) for r in uc.read_tsv(POOL)}
    cnt = defaultdict(Counter)
    strand_mismatch = 0
    with gzip.open(a.rows, 'rt') as f, gzip.open(a.out_table, 'wt') as fo, gzip.open(a.out_compact, 'wt') as fc:
        hdr = f.readline().rstrip('\n').split('\t')
        fo.write('\t'.join(hdr + EVID) + '\n')
        fc.write('\t'.join(KEEP + EVID) + '\n')
        for line in f:
            r = dict(zip(hdr, line.rstrip('\n').split('\t')))
            e = {}
            for s in 'ab':
                smp, occ, ctg = r[s + '_sample'], r[s + '_occurrence'], r[s + '_contig']
                lo, hi = r[s + '_x_interval'].split('-')
                hit = pos.get((smp, ctg, lo, hi))
                if hit is not None and occ not in hit['occurrences'].split(';'):
                    hit = None  # same coordinates tested for another occurrence: not this occurrence's result
                if hit is not None and r5strand.get(occ, r[s + '_strand']) != r[s + '_strand']:
                    strand_mismatch += 1
                    hit = None
                other = [t for t in by_occ.get(occ, []) if hit is None or t != hit['test_id']]
                if hit is not None:
                    st = hit['position_class']
                    e.update({'rc_%s_test_id' % s: hit['test_id'], 'rc_%s_verdict' % s: hit['verdict'],
                              'rc_%s_count_mode' % s: hit['count_mode'], 'rc_%s_P_matched' % s: hit['P_matched'],
                              'rc_%s_A_matched' % s: hit['A_matched'],
                              'rc_%s_locus_ambiguous' % s: hit['locus_ambiguous']})
                elif other:
                    st = 'tested_other_interval'
                elif smp in r5samples or usable.get(smp, 0) > 0:
                    st = 'not_tested'
                else:
                    st = 'no_reads'
                e['rc_%s_status' % s] = st
                e['rc_%s_other_intervals' % s] = ';'.join(
                    '%s:%s-%s:%s' % (t, tid[t]['lo'], tid[t]['hi'], tid[t]['position_class']) for t in other
                    if t in tid)
                if smp not in searched:
                    rs = 'sample_not_searched'
                elif occ not in panel.get(smp, ()):
                    rs = 'not_in_search_panel'
                else:
                    m = found.get((smp, occ), [])
                    rs = 'recovered' if m and max(m) >= MAIN_COUNT else ('seen_below_threshold' if m else
                                                                         'not_recovered')
                    e['rss_%s_molecules' % s] = max(m) if m else 0
                e['rss_%s_status' % s] = rs
                e['rss_%s_batch' % s] = searched.get(smp, '')
            e['rc_pair'] = '|'.join(sorted((e['rc_a_status'], e['rc_b_status'])))
            e['rss_pair'] = '|'.join(sorted((e['rss_a_status'], e['rss_b_status'])))
            fo.write(line.rstrip('\n') + '\t' + '\t'.join(str(e.get(c, '')) for c in EVID) + '\n')
            fc.write('\t'.join(str(r.get(c, '')) for c in KEEP) + '\t' + '\t'.join(str(e.get(c, '')) for c in EVID)
                     + '\n')
            scope = r['pair_scope']
            cnt['pairs'][scope] += 1
            for s in 'ab':
                cnt['rc_side_status'][e['rc_%s_status' % s]] += 1
                cnt['rss_side_status'][e['rss_%s_status' % s]] += 1
            cnt['rc_pair'][e['rc_pair']] += 1
            cnt['rc_pair|' + scope][e['rc_pair']] += 1
            cnt['rss_pair'][e['rss_pair']] += 1
            cnt['final_status_unchanged'][r['final_status']] += 1
    summ = dict(version='ucd_read_evidence_v0.1', rows=a.rows, rows_sha256=ub.sha_file(a.rows),
                sources=dict(round5_v02={str(x): ub.sha_file(x) for x in (R5 / 'positions.tsv.gz',
                                                                         R5 / 'interval_map.tsv.gz',
                                                                         R5 / 'pairs.tsv.gz')},
                             read_site_search=[b[0] for b in RSS_BATCHES], pool=str(POOL),
                             manifest_sha256=ub.sha_file(REPO / 'manifest.tsv')),
                rule=dict(rc='exact sample + contig + lo/hi + strand; no inheritance across intervals',
                          rss='flank %d, >= %d molecules, known_catalog_position*' % (MAIN_FLANK, MAIN_COUNT)),
                strand_mismatch_rejected=strand_mismatch,
                counts={k: dict(v.most_common()) for k, v in cnt.items()},
                out_table=dict(path=a.out_table, sha256=ub.sha_file(a.out_table)),
                out_compact=dict(path=a.out_compact, sha256=ub.sha_file(a.out_compact)),
                at=time.strftime('%Y-%m-%d %H:%M:%S'))
    Path(a.summary).write_text(json.dumps(summ, indent=1) + '\n')
    print(json.dumps({k: summ['counts'][k] for k in ('pairs', 'rc_side_status', 'rss_side_status',
                                                     'final_status_unchanged')}, indent=1))


if __name__ == '__main__':
    main()
