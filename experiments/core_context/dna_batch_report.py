"""Small per-batch report of the DNA entry (one tag): input scale, discovery, pairing (computed / cached / deferred),
the two checks (new / reused / retired; status counts), strict < 2500 bp passes (pairs, occurrences, assemblies;
same / cross assembly), resources (from /usr/bin/time files of this batch's jobs) and failed / retried units.
Counts of pairs, occurrences and assemblies are kept apart; candidates checked are not passes."""
import argparse
import glob
import gzip
import json
import re
from collections import Counter
from pathlib import Path

import ucd_compare as uc


def jl(p):
    p = Path(p)
    return json.loads(p.read_text()) if p.exists() else None


def secs(t):
    p = [float(x) for x in t.split(':')]
    return p[0] * 3600 + p[1] * 60 + p[2] if len(p) == 3 else p[0] * 60 + p[1]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    for x in ('--work', '--manifest', '--tag', '--gens', '--collect', '--ucd', '--final', '--jobs', '--out'):
        ap.add_argument(x, required=True)
    a = ap.parse_args()
    w = Path(a.work)
    gens = [int(g) for g in a.gens.split(',')]
    man = uc.read_tsv(a.manifest)
    chunks = [r for r in uc.read_tsv(w / 'chunks.tsv') if r['tag'] == a.tag]
    batch_samples = {s for r in chunks for s in r['samples'].split(';')}
    all_samples = {s for r in uc.read_tsv(w / 'chunks.tsv') for s in r['samples'].split(';')}
    seeds = [jl(w / 'seeds' / r['chunk'] / 'done.json') for r in chunks]
    lib = [jl(p) for p in sorted(glob.glob(str(w / 'library' / 'lib_*.json')))]
    lib_tag = [x for x in lib if x and x.get('tag') == a.tag]
    occ = [jl(p) for p in glob.glob(str(w / 'occ' / 'lib*' / '*' / 'done.json'))]
    R = dict(tag=a.tag, input=dict(batch_assemblies=len(batch_samples), cumulative_assemblies=len(all_samples),
                                   manifest_assemblies=len(man),
                                   sha256_mismatch=sorted({s for x in seeds if x for s in x['sha256_mismatch']})),
             seeds=dict(chunks=len(chunks), done=sum(x is not None for x in seeds),
                        seed_loci=sum(x['seed_loci'] for x in seeds if x),
                        representatives=sum(x['representatives'] for x in seeds if x)),
             library=dict(shards_total=len(lib), this_batch=lib_tag,
                          sequences_total=sum(x['new_library_sequences'] for x in lib if x)),
             occurrences=dict(tasks_done=len(occ), hits=sum(x['hits'] for x in occ)))
    en = {}
    for g in gens:
        e = [jl(p) for p in glob.glob(str(w / 'enum' / ('g%d' % g) / 'b*' / 'done.json'))]
        al = [jl(p) for p in glob.glob(str(w / 'align' / ('g%d' % g) / 'k*' / 'done.json'))]
        en[g] = dict(enumerate_tasks=len(e), pairs_enumerated=sum(x['pairs_enumerated'] for x in e),
                     pairs_deferred_this_gen=sum(x['pairs_deferred'] for x in e),
                     align_tasks=len(al), sources=sum(x['sources'] for x in al),
                     keys_computed=sum(x['keys_computed'] for x in al), keys_cached=sum(x['keys_cached'] for x in al),
                     blast_calls=sum(x['blast_calls'] for x in al))
    R['generations'] = en
    c = jl(Path(a.collect) / 'collect.json')
    if c:
        c.pop('params', None)
        R['collect'] = {k: v for k, v in c.items() if 'sha256' not in k}
    f = jl(Path(a.ucd) / 'freeze.json') or {}
    R['checks_input'] = dict(catalogue_pairs=f.get('catalogue_pairs'), checked_now=f.get('input_pairs'),
                             pair_status=f.get('pair_status'), chunks=f.get('chunks'))
    R['checks_completeness'] = jl(Path(a.ucd) / 'completeness.json')
    # final rows
    st, src, ucd, fin = Counter(), Counter(), Counter(), Counter()
    socc, sasm, sscope = set(), set(), Counter()
    with gzip.open(Path(a.final) / 'all_rows_candidate_insertion.tsv.gz', 'rt') as fh:
        h = fh.readline().rstrip('\n').split('\t')
        for line in fh:
            r = dict(zip(h, line.rstrip('\n').split('\t')))
            src[r['source_check']] += 1
            ucd[r['ucd_check']] += 1
            fin[r['final_status']] += 1
    with gzip.open(Path(a.final) / 'strict_lt2500_candidate_insertion.tsv.gz', 'rt') as fh:
        h = fh.readline().rstrip('\n').split('\t')
        n = 0
        for line in fh:
            r = dict(zip(h, line.rstrip('\n').split('\t')))
            n += 1
            sscope[r['pair_scope']] += 1
            for s in 'ab':
                socc.add((r[s + '_sample'], r[s + '_contig'], r[s + '_x_interval']))
                sasm.add(r[s + '_sample'])
    R['final'] = dict(rows=sum(fin.values()), source_check=dict(src), ucd_check=dict(ucd), final_status=dict(fin),
                      strict_lt2500_pairs=n, strict_by_scope=dict(sscope), strict_occurrences=len(socc),
                      strict_assemblies=len(sasm), discovery_source='dna',
                      note='cumulative over all batches processed in this work dir')
    tj = jl(Path(a.final) / 'target_700_2500_summary.json')
    R['target_700_2500'] = tj['main'] if tj else None  # main table; final['strict_*'] is the wide < 2500 table
    # resources
    res = Counter()
    peak = {}
    for p in glob.glob(str(w / 'time' / '*.time')):
        stage = Path(p).name.split('_')[0]
        t = Path(p).read_text()
        m = re.search(r'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): (\S+)', t)
        r = re.search(r'Maximum resident set size \(kbytes\): (\d+)', t)
        if m:
            res[stage] += secs(m.group(1))
        if r:
            peak[stage] = max(peak.get(stage, 0), int(r.group(1)))
    R['resources_cumulative'] = {k: dict(task_hours=round(v / 3600, 2), peak_rss_gb=round(peak.get(k, 0) / 1e6, 1))
                                 for k, v in res.items()}
    jobs = [l.rstrip('\n').split('\t') for l in open(a.jobs)] if Path(a.jobs).exists() else []
    R['submissions'] = Counter(j[1] for j in jobs)
    R['retries'] = {k: v - 1 for k, v in R['submissions'].items() if v > 1}
    pend = sum(int(r['pairs']) for g in gens[-1:] for p in glob.glob(str(w / 'enum' / ('g%d' % g) / 'b*' /
                                                                        'progress.tsv.gz'))
               for r in uc.read_tsv(p) if r['status'] == 'deferred')
    R['unfinished'] = dict(deferred_pairs_pending=pend, complete=pend == 0)
    Path(a.out).write_text(json.dumps(R, indent=1, default=dict) + '\n')
    print(json.dumps(R, indent=1, default=dict))


if __name__ == '__main__':
    main()
