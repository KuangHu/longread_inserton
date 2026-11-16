"""Concordance of two DNA-entry pair runs on the same occurrences (e.g. v0.1.2 pairs vs v0.2 enumerate/align):
every source pair (library_id, a_occurrence, b_occurrence) of the earlier run is matched in the later run and its
status and X intervals are compared. Reports identical / changed (by kind) / missing, and pairs only in the later run
(previously deferred). Used to verify that the de-duplicated alignment does not drop or alter candidates."""
import argparse
import gzip
import json
from collections import Counter
from pathlib import Path


def stream(paths):
    for p in paths:
        with gzip.open(p, 'rt') as f:
            h = f.readline().rstrip('\n').split('\t')
            for line in f:
                yield dict(zip(h, line.rstrip('\n').split('\t')))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--old', required=True, help='dna_pairs_all.tsv.gz of the earlier run')
    ap.add_argument('--new-work', required=True)
    ap.add_argument('--gen', default='0')
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    old = {}
    for r in stream([a.old]):
        old[(r['library_id'], r['a_occurrence'], r['b_occurrence'])] = (r['status'], r['a_x'], r['b_x'])
    c, ex = Counter(), []
    seen = set()
    new_files = sorted(Path(a.new_work).glob('align/g%s/k*/pairs.tsv.gz' % a.gen))
    for r in stream(new_files):
        k = (r['library_id'], r['a_occurrence'], r['b_occurrence'])
        o = old.get(k)
        if o is None:
            c['only_new:%s' % r['status']] += 1
            continue
        seen.add(k)
        n = (r['status'], r['a_x'], r['b_x'])
        if n == o:
            c['identical:%s' % o[0]] += 1
        elif o[0] != n[0]:
            c['status_changed:%s->%s' % (o[0], n[0])] += 1
            if len(ex) < 50:
                ex.append(dict(key=k, old=o, new=n, boundary=(r['boundary_up'], r['boundary_down'])))
        else:
            c['x_changed:%s' % o[0]] += 1
            if len(ex) < 50:
                ex.append(dict(key=k, old=o, new=n))
    c['missing_in_new'] = len(old) - len(seen)
    res = dict(old_sources=len(old), new_files=len(new_files), counts=dict(sorted(c.items())), examples=ex[:20])
    Path(a.out).write_text(json.dumps(res, indent=1, default=list) + '\n')
    print(json.dumps({k: v for k, v in res.items() if k != 'examples'}, indent=1))


if __name__ == '__main__':
    main()
