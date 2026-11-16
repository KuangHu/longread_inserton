"""v0.1 -> v0.2 class of every main-setting position (flank 1000, panel cores), and the round-5-supported recovery
misses of v0.2 by reason. A v0.1 position is matched to a v0.2 position of the same sample and core whose canonical
flank ends are within TOL; otherwise the v0.2 observations of that core with a flank end within TOL of either v0.1
end give its fate (their status). Writes v01_to_v02.tsv and recovery_misses_v02.tsv here; prints the tables.
Usage: python compare_v01_v02.py   (from experiments/core_context, after combine; v0.1 per-sample tables in the
work dir's v0.1_backup/)
"""
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import protein_context as pc  # noqa: E402
import read_site_search as rs  # noqa: E402

HERE = Path(__file__).resolve().parent
L, N = rs.MAIN_FLANK, rs.MAIN_COUNT


def rows_of(path, sample=None):
    for r in pc.read_tsv(path):
        if sample:
            r['sample'] = sample
        if int(r['flank']) == L and r['kind'] == 'core':
            yield r


def coarse(cls):
    return cls.split(':')[0] if cls.startswith(('unplaceable', 'homolog', 'catalog_omission')) else cls


def main():
    samples = rs.pilot_samples(HERE)
    old = [r for s in samples for r in rows_of(rs.WORK / 'v0.1_backup' / s / 'positions.tsv.gz', s)]
    new = list(rows_of(HERE / 'positions.tsv.gz'))
    placed = defaultdict(list)
    for r in new:
        if r['left_pos'] not in ('', '-1'):
            placed[(r['sample'], r['core_key'], r['left_contig'])].append(r)
    obs = {}
    out, trans = [], Counter()
    for r in old:
        if r['cls'].startswith(('unplaceable', 'one_sided')) or int(r['molecules']) < N:
            continue
        lp, rp = int(r['left_pos']), int(r['right_pos'])
        m = [x for x in placed[(r['sample'], r['core_key'], r['left_contig'])] if x['right_contig'] == r['right_contig']
             and abs(int(x['left_pos']) - lp) <= rs.TOL and abs(int(x['right_pos']) - rp) <= rs.TOL]
        if m:
            x = max(m, key=lambda x: int(x['molecules']))
            fate = x['cls'] + ('' if int(x['molecules']) >= N else ' (<N)')
        else:
            if r['sample'] not in obs:
                obs[r['sample']] = [o for o in pc.read_tsv(rs.WORK / 'samples' / r['sample'] / 'observations.tsv.gz')
                                    if int(o['flank']) == L]
            st = Counter('%s|%s' % (o['left_status'], o['right_status']) for o in obs[r['sample']]
                         if o['core_key'] == r['core_key'] and o['status'] != 'one_sided' and any(
                             o[k + '_contig'] == c and o[k + '_pos'] not in ('', '-1') and abs(int(o[k + '_pos']) - p) <= rs.TOL
                             for k in ('left', 'right') for c, p in ((r['left_contig'], lp), (r['right_contig'], rp))))
            fate = 'observations:' + (st.most_common(1)[0][0] if st else 'none')
        t = (r['cls'] if r['cls'] != 'shifted_known_position' else 'shifted_known_position', fate)
        trans[(coarse(t[0]), coarse(fate.split(' ')[0]) if not fate.startswith('observations') else fate)] += 1
        out.append(dict(sample=r['sample'], core_key=r['core_key'], left_contig=r['left_contig'], left_pos=lp,
                        right_contig=r['right_contig'], right_pos=rp, molecules=r['molecules'], v01_cls=r['cls'],
                        v02=fate))
    pc.write_tsv(HERE / 'v01_to_v02.tsv', out, list(out[0]))
    for (a, b), n in sorted(trans.items(), key=lambda x: (x[0][0], -x[1])):
        print('%-45s -> %-60s %5d' % (a, b, n))
    # recovery misses (round-5 supported catalog positions not recovered at L/N)
    import json
    r5 = defaultdict(set)
    for p in pc.read_tsv(rs.ROUND5 / 'interval_map.tsv.gz'):
        r5[p['occurrence']].add(p['position_class'])
    cat = [c for c in pc.read_tsv(rs.WORK / 'panel' / 'catalog.tsv') if c['sample'] in samples]
    cores = {(m['structure'], m['core_key']) for m in pc.read_tsv(rs.WORK / 'panel' / 'members.tsv')}
    got = {o for r in new if r['cls'].startswith('known_catalog_position') and int(r['molecules']) >= N
           for o in r['catalog'].split(';')}
    miss = []
    for c in cat:
        if (c['structure'], c['core_key']) not in cores or 'supported' not in r5.get(c['occurrence_id'], ()) \
                or c['occurrence_id'] in got:
            continue
        s0, e0 = int(c['start0']), int(c['end0'])
        near = [r for r in new if r['sample'] == c['sample'] and (
            r['contig'] == c['contig'] and r['end'] != '' and int(r['start']) < e0 + rs.TOL and int(r['end']) > s0 - rs.TOL
            or not r['contig'] and any(r[k + '_contig'] == c['contig'] and r[k + '_pos'] not in ('', '-1') and
                                       s0 - rs.TOL <= int(r[k + '_pos']) <= e0 + rs.TOL for k in ('left', 'right')))]
        same = [r for r in near if r['core_key'] == c['core_key']]
        mine = [r for r in same if c['occurrence_id'] in r['catalog'].split(';')]
        if mine:
            why = 'found_below_N' if all(int(r['molecules']) < N for r in mine) else 'found_as_' + mine[0]['cls']
        elif same:
            why = 'same_core_nearby:' + ','.join(sorted({r['cls'] for r in same}))
        elif near:
            why = 'other_cores_only:' + ','.join(sorted({coarse(r['cls']) for r in near}))
        else:
            why = 'nothing_nearby'
        miss.append(dict(c, reason=why, nearby=';'.join('%s:%s:%s:%s' % (r['core_key'], r['cls'], r['molecules'],
                                                                         r['catalog_unmatched'][:80]) for r in near[:6])))
    pc.write_tsv(HERE / 'recovery_misses_v02.tsv', miss, list(miss[0]) if miss else ['occurrence_id'])
    print(len(miss), json.dumps(Counter(m['reason'] for m in miss).most_common(), indent=0))


if __name__ == '__main__':
    main()
