"""Download plan for the core-context read check (run from the repo root, no network).

Usable single-file runs in the pool that are not in the ready manifest are ordered greedily by the number of new cores
that would get an A or B pair (ties: new A/B pairs, then smaller download), then cut into batches of about 10 GB.
Runs that add neither a core nor an A/B pair are left out. Writes plan.tsv and batch<n>_runs.tsv (pool runs.tsv rows)."""
import sys, csv
from collections import defaultdict
sys.path.insert(0,'experiments/core_context')
import ucd_readcheck as u
rows=list(u.pc.read_tsv(u.POOL/'runs.tsv'))
haveS=set(u.run_info())
ready=u.reads_ready.load(u.READS/'ready_manifest.tsv')
miss=defaultdict(list)
for r in rows:
    if r['usable']=='True' and r['run_accession'] not in ready and ';' not in r['fastq_ftp']:
        miss[r['sample_accession']].append(r)
out=u.Path('experiments/core_context/scale_comamonadaceae_v0.1')
pairs=list(u.pc.read_tsv(out/'match_both_lt2500_pairs_b1-2-3.tsv.gz'))
by=defaultdict(list)
for i,p in enumerate(pairs):
    by[p['sample_a']].append(i); by[p['sample_b']].append(i)
S=set(haveS)
def ab(p,S): return u.tier_of(p['pair_scope'],p['sample_a'] in S,p['sample_b'] in S) in ('A','B')
covered={(p['structure'],p['core_key']) for p in pairs if ab(p,S)}
gb={s:sum(int(r['fastq_bytes']) for r in rs)/1e9 for s,rs in miss.items()}
order=[]; left=set(miss)
while left:
    best=None
    for s in left:
        T=S|{s}
        nc={(pairs[i]['structure'],pairs[i]['core_key']) for i in by[s] if ab(pairs[i],T)}-covered
        npairs=sum(ab(pairs[i],T) and not ab(pairs[i],S) for i in by[s])
        k=(len(nc),npairs,-gb[s])
        if best is None or k>best[0]: best=(k,s,nc)
    k,s,nc=best
    if k[0]==0 and k[1]==0: break
    order.append((s,k[0],k[1],gb[s])); S.add(s); covered|=nc; left.remove(s)
print('samples with gain',len(order),'GB',round(sum(o[3] for o in order),1),'new cores',sum(o[1] for o in order),'skipped (no gain)',len(left),'GB',round(sum(gb[s] for s in left),1))
batches=[];cur=[];size=0
for o in order:
    if cur and size+o[3]>10: batches.append(cur);cur=[];size=0
    cur.append(o);size+=o[3]
batches.append(cur)
d=out/'download_batches'
cols=['batch','rank','sample_accession','run_accession','instrument_platform','fastq_bytes','new_cores','new_AB_pairs']
rank=0
with open(str(d/'plan.tsv'),'w') as h:
    w=csv.writer(h,delimiter='\t',lineterminator='\n'); w.writerow(cols)
    for b,bl in enumerate(batches,1):
        runs=[]
        for s,nc,npr,g in bl:
            rank+=1
            for r in miss[s]:
                w.writerow([b,rank,s,r['run_accession'],r['instrument_platform'],r['fastq_bytes'],nc,npr]); runs.append(r)
        with open(str(d/('batch%d_runs.tsv'%b)),'w') as f:
            ww=csv.DictWriter(f,list(rows[0].keys()),delimiter='\t',lineterminator='\n'); ww.writeheader(); ww.writerows(runs)
        print('batch',b,'samples',len(bl),'runs',len(runs),'GB',round(sum(x[3] for x in bl),2),'new cores',sum(x[1] for x in bl),'new A/B pairs',sum(x[2] for x in bl))
