"""Local junction diagnostic for one split-only DEL; does not change calls."""
import collections,json,pathlib,sys
import pysam
from site_recovery.core import aligned_slice,identity,usable_read,Config,MATCH,REF_CONSUME,QUERY_CONSUME
root=pathlib.Path(sys.argv[1]);d=next(d for d in json.loads((root/'audit/candidate_audit.json').read_text()) if d['candidate']['id']=='Sniffles2.DEL.4SB6');c=d['candidate'];wanted={g['read'] for g in d['split_geometry']};byname=collections.defaultdict(list)
with pysam.AlignmentFile(str(root/'calls/alignments.bam'),'rb') as b:
 for r in b.fetch(until_eof=True):
  if r.query_name in wanted:byname[r.query_name].append(r)
def qp(r,pos):
 rp=r.reference_start;qq=r.cigartuples[0][1] if r.cigartuples[0][0]==5 else 0
 for op,n in r.cigartuples:
  if op in MATCH and rp<=pos<rp+n:return qq+pos-rp
  if op in REF_CONSUME:rp+=n
  if op in QUERY_CONSUME:qq+=n
 return None
examples=[]
with pysam.FastaFile(str(root/'reference/marine_mag_pool.fasta')) as ref,pysam.AlignmentFile(str(root/'calls/alignments.bam'),'rb') as bam:
 for g in d['split_geometry']:
  if abs(g['query_gap'])>10:continue
  rr=byname[g['read']];a,b=g['left'],g['right'];left=next(r for r in rr if r.reference_start==a[0] and r.reference_end==a[1]);right=next(r for r in rr if r.reference_start==b[0] and r.reference_end==b[1]);primary=next(r for r in rr if not r.is_secondary and not r.is_supplementary)
  if any(op==5 for op,n in primary.cigartuples):continue
  q1=qp(left,a[1]-100);q2=qp(right,b[0]+99)
  if q1 is None or q2 is None:continue
  query=primary.query_sequence[q1:q2+1].upper();empty=ref.fetch(c['contig'],a[1]-100,a[1]).upper()+ref.fetch(c['contig'],b[0],b[0]+100).upper()
  l=aligned_slice(left,a[1]-100,a[1]);r=aligned_slice(right,b[0],b[0]+100)
  examples.append(dict(read=g['read'],start0=a[1],end0=b[0],query_gap=g['query_gap'],junction_query_length=len(query),junction_identity=identity(query,empty),left_anchor_identity=identity(l or '',empty[:100]),right_anchor_identity=identity(r or '',empty[100:])))
 groups=collections.Counter((e['start0'],e['end0']) for e in examples);start,end=groups.most_common(1)[0][0];occ=[]
 for r in bam.fetch(c['contig'],start-100,end+100):
  if not usable_read(r,Config()):continue
  s=aligned_slice(r,start-100,end+100)
  if s and identity(s,ref.fetch(c['contig'],start-100,end+100).upper())>=.95:occ.append(r.query_name)
res=dict(candidate=c,split_junction_checks=examples,dominant_reference_interval=[start,end],occupied_reference_like_reads=sorted(set(occ)),caution='Audit evidence only. Junction tested across 100-bp anchors; no full element classification or formal rescue performed.')
(root/'audit/split_example_local_check.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res,indent=2))
