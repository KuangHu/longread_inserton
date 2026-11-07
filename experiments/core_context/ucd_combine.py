"""Combine ucd_batch runs into one result per interval kind, with completeness, counts and the strict length table.

  rerun-report  rows of a rerun vs the earlier run:
                - affected pairs: verdicts before / after, new passes;
                - control pairs (unaffected, rerun only to compare): must be unchanged in every compared column.
  combine       final rows per interval kind (streamed) from:
                - --reuse: earlier all-rows files whose pairs are not in --replace-ids, tagged --reuse-version;
                - --replace: a rerun's rows for exactly the --replace-ids;
                - --runs: further runs, taken whole.
                Checked per pair_id x interval_kind against --ids: no missing, no duplicated and no unexpected rows,
                and every row has a status. Writes, per kind:
                - all rows;
                - the pass table (source_check pass and ucd_check ucd_supported);
                - the strict table: pass, and the X actually extracted on both sides has 0 < length < 2500 bp,
                  computed from the Check-1 interval and sequence (an old lt2500 label is not used; exactly 2500 bp
                  fails).
                Counts are pairs, distinct occurrences, samples and cores, kept apart, never mixed units.
"""
import argparse
import gzip
import json
import time
from collections import Counter, defaultdict
from pathlib import Path

import ucd_batch as ub
import ucd_compare as uc

STRICT_MAX = 2500
COMPARE = [c for c in ub.ROW_COLS if c not in ('code_version',)]


def rows_of(path):
    with gzip.open(path, 'rt') as f:
        hdr = f.readline().rstrip('\n').split('\t')
        for line in f:
            yield dict(zip(hdr, line.rstrip('\n').split('\t')))


def ids_file(p):
    return ub.ids_of(p)


def cmd_rerun_report(a):
    old_dir, new_dir, out = Path(a.old_work), Path(a.new_work), Path(a.out)
    aff, ctl = ids_file(a.affected), ids_file(a.control)
    rep = {}
    for kind in ub.KINDS:
        old = {r['pair_id']: r for r in rows_of(old_dir / ('ucd_all_rows_%s.tsv.gz' % kind)) if r['pair_id'] in aff | ctl}
        new = {r['pair_id']: r for r in rows_of(new_dir / ('rows_%s.tsv.gz' % kind))}
        assert set(new) == aff | ctl, 'rerun rows do not match the requested IDs'
        trans = Counter('%s / %s -> %s / %s' % (old[i]['source_check'], old[i]['ucd_check'], new[i]['source_check'],
                                                new[i]['ucd_check']) for i in aff)
        diffs = Counter()
        changed = []
        for i in ctl:
            d = [c for c in COMPARE if old[i].get(c, '') != new[i].get(c, '')]
            for c in d:
                diffs[c] += 1
            if d:
                changed.append(dict(pair_id=i, columns=d))
        rep[kind] = dict(
            affected_pairs=len(aff),
            affected_source_check_before=dict(Counter(old[i]['source_check'] for i in aff)),
            affected_source_check_after=dict(Counter(new[i]['source_check'] for i in aff)),
            affected_transitions=dict(sorted(trans.items())),
            affected_new_pass=sum(new[i]['final_status'] == 'pass' and old[i]['final_status'] != 'pass' for i in aff),
            affected_remaining_source_reasons=dict(Counter(new[i]['source_reasons'] for i in aff
                                                           if new[i]['source_check'] != 'pass')),
            control_pairs=len(ctl), control_rows_identical=len(ctl) - len(changed),
            control_columns_differing=dict(diffs), control_changed_examples=changed[:10])
    (out / 'rerun_report.json').write_text(json.dumps(rep, indent=1) + '\n')
    # one row per affected pair and interval kind: before / after
    tab = []
    for kind in ub.KINDS:
        old = {r['pair_id']: r for r in rows_of(old_dir / ('ucd_all_rows_%s.tsv.gz' % kind)) if r['pair_id'] in aff}
        new = {r['pair_id']: r for r in rows_of(new_dir / ('rows_%s.tsv.gz' % kind)) if r['pair_id'] in aff}
        for i in sorted(aff):
            o, n = old[i], new[i]
            tab.append(dict(pair_id=i, interval_kind=kind, pair_scope=n['pair_scope'], structure=n['structure'],
                            core_key=n['core_key'], a_occurrence=n['a_occurrence'], b_occurrence=n['b_occurrence'],
                            old_source_check=o['source_check'], old_source_reasons=o['source_reasons'],
                            old_ucd_check=o['ucd_check'], old_final_status=o['final_status'],
                            new_source_check=n['source_check'], new_source_reasons=n['source_reasons'],
                            new_ucd_check=n['ucd_check'], new_final_status=n['final_status'],
                            new_x_len_a=n['x_len_a'], new_x_len_b=n['x_len_b']))
    cols = list(tab[0]) if tab else ['pair_id']
    uc.write_tsv(out / 'rerun_affected.tsv', tab, cols)
    print(json.dumps({k: {x: v[x] for x in ('affected_transitions', 'affected_new_pass', 'control_rows_identical',
                                            'control_columns_differing')} for k, v in rep.items()}, indent=1))


class Tally:
    """Pairs, distinct occurrences, samples and cores of one category."""

    def __init__(self):
        self.p, self.o, self.s, self.c = 0, set(), set(), set()

    def add(self, r):
        self.p += 1
        self.o.update((r['a_occurrence'], r['b_occurrence']))
        self.s.update((r['a_sample'], r['b_sample']))
        self.c.add((r['structure'], r['core_key']))

    def out(self):
        return dict(pairs=self.p, occurrences=len(self.o), samples=len(self.s), cores=len(self.c))


def strict_len(r):
    """(ok, [failing sides]) from the Check-1 X interval and the X length actually extracted."""
    bad, err = [], []
    for s in 'ab':
        iv, n = r.get(s + '_x_interval', ''), r.get('x_len_' + s, '')
        if not iv or n == '':
            err.append(s)
            continue
        lo, hi = map(int, iv.split('-'))
        if hi - lo != int(n):
            err.append(s)
        elif not (0 < int(n) < STRICT_MAX):
            bad.append(s)
    return bad, err


def norm_reason(x):
    side, _, rest = x.partition(':')
    return rest.split('@')[0] if rest else x


def cmd_combine(a):
    out, work = Path(a.out), Path(a.work)
    out.mkdir(parents=True, exist_ok=True)
    work.mkdir(parents=True, exist_ok=True)
    expect = ids_file(a.ids)
    replace = ids_file(a.replace_ids) if a.replace_ids else set()
    sources = []
    if a.reuse:
        sources.append(('reuse', Path(a.reuse), 'ucd_all_rows_%s.tsv.gz', lambda r: r['pair_id'] not in replace,
                        a.reuse_version))
    if a.replace:
        sources.append(('replace', Path(a.replace), 'rows_%s.tsv.gz', lambda r: r['pair_id'] in replace, None))
    for w in filter(None, a.runs.split(',')):
        sources.append(('run', Path(w), 'rows_%s.tsv.gz', lambda r: True, None))
    res, comp_all, files = {}, {}, {}
    for kind in ub.KINDS:
        seen, dup, extra, nostatus = set(), 0, 0, 0
        T = defaultdict(Tally)
        cnt = defaultdict(Counter)
        strict_ex_pairs, strict_ex_sides, strict_ex_iv, len_err = 0, 0, set(), 0
        paths = {n: work / ('%s_%s.tsv.gz' % (n, kind)) for n in ('all_rows', 'pass', 'strict_lt2500')}
        fo = {n: gzip.open(p, 'wt') for n, p in paths.items()}
        for f in fo.values():
            f.write('\t'.join(ub.ROW_COLS) + '\n')
        for name, d, pat, keep, ver in sources:
            n_in = 0
            for r in rows_of(d / (pat % kind)):
                if not keep(r):
                    continue
                n_in += 1
                if ver and not r.get('code_version'):
                    r['code_version'] = ver
                pid = r['pair_id']
                if pid not in expect:
                    extra += 1
                    continue
                if pid in seen:
                    dup += 1
                    continue
                seen.add(pid)
                if not r.get('final_status') or not r.get('source_check') or not r.get('ucd_check'):
                    nostatus += 1
                line = '\t'.join(str(r.get(c, '')) for c in ub.ROW_COLS) + '\n'
                fo['all_rows'].write(line)
                scope = r['pair_scope']
                cnt['code_version'][r['code_version']] += 1
                cnt['source_check'][r['source_check']] += 1
                if r['source_check'] == 'fail':
                    for x in r['source_reasons'].split(';'):
                        cnt['source_failure_reason (sides)'][norm_reason(x)] += 1
                if r['source_check'] == 'not_run':
                    for x in r['source_reasons'].split(';'):
                        cnt['not_run_reason (sides)'][norm_reason(x)] += 1
                    T['not_run:' + ('candidate_interval_not_available' if 'candidate_interval_not_available' in
                                    r['source_reasons'] else 'input_error')].add(r)
                T['ucd_check:' + r['ucd_check']].add(r)
                T['ucd_check:%s|%s' % (r['ucd_check'], scope)].add(r)
                if r['final_status'] == 'pass':
                    fo['pass'].write(line)
                    for g in ('pass', 'pass|' + scope, 'pass|%s|%s' % (scope, r['structure'])):
                        T[g].add(r)
                    bad, err = strict_len(r)
                    if err:
                        len_err += 1
                    elif not bad:
                        fo['strict_lt2500'].write(line)
                        for g in ('strict', 'strict|' + scope, 'strict|%s|%s' % (scope, r['structure'])):
                            T[g].add(r)
                    else:
                        strict_ex_pairs += 1
                        strict_ex_sides += len(bad)
                        strict_ex_iv |= {(r[s + '_occurrence'], r[s + '_x_interval']) for s in bad}
                        T['strict_excluded'].add(r)
            cnt['rows_by_source']['%s:%s' % (name, d)] += n_in
        for f in fo.values():
            f.close()
        missing = len(expect - seen)
        comp_all[kind] = dict(expected=len(expect), observed=len(seen), missing=missing, duplicated=dup,
                              unexpected=extra, without_status=nostatus, x_length_inconsistent=len_err,
                              complete=not (missing or dup or extra or nostatus or len_err))
        res[kind] = dict(counts={k: dict(v.most_common()) for k, v in cnt.items()},
                         tallies={k: v.out() for k, v in sorted(T.items())},
                         strict=dict(rule='pass and 0 < X length < %d bp on both sides (from the extracted interval)' %
                                          STRICT_MAX,
                                     before=T['pass'].out(), after=T['strict'].out(),
                                     excluded_pairs=strict_ex_pairs, excluded_side_intervals=strict_ex_sides,
                                     excluded_distinct_side_intervals=len(strict_ex_iv)))
        for n, p in paths.items():
            files[p.name] = dict(path=str(p), bytes=p.stat().st_size, sha256=ub.sha_file(p))
    summ = dict(version='ucd_combine_v0.1', ids=a.ids, ids_sha256=ub.sha_file(a.ids), expected_pairs=len(expect),
                replace_ids=a.replace_ids, replaced_pairs=len(replace), reuse=a.reuse, reuse_version=a.reuse_version,
                replace=a.replace, runs=a.runs, completeness=comp_all, by_kind=res, files=files,
                at=time.strftime('%Y-%m-%d %H:%M:%S'))
    (out / 'combined_summary.json').write_text(json.dumps(summ, indent=1) + '\n')
    print(json.dumps(dict(completeness=comp_all, strict={k: v['strict'] for k, v in res.items()}), indent=1))
    if not all(c['complete'] for c in comp_all.values()):
        raise SystemExit('combined result incomplete or inconsistent')


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = p.add_subparsers(dest='cmd', required=True)
    x = sp.add_parser('rerun-report')
    x.add_argument('--out', required=True)
    x.add_argument('--old-work', required=True)
    x.add_argument('--new-work', required=True)
    x.add_argument('--affected', required=True)
    x.add_argument('--control', required=True)
    x = sp.add_parser('combine')
    x.add_argument('--out', required=True)
    x.add_argument('--work', required=True)
    x.add_argument('--ids', required=True, help='expected pair IDs (catalogue)')
    x.add_argument('--reuse', default='', help='earlier work dir with ucd_all_rows_<kind>.tsv.gz')
    x.add_argument('--reuse-version', default='')
    x.add_argument('--replace', default='', help='rerun work dir with rows_<kind>.tsv.gz')
    x.add_argument('--replace-ids', default='')
    x.add_argument('--runs', default='', help='comma-separated run work dirs taken whole')
    a = p.parse_args()
    dict(**{'rerun-report': cmd_rerun_report, 'combine': cmd_combine})[a.cmd](a)


if __name__ == '__main__':
    main()
