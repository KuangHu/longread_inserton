"""Inspect selected failure examples without rescuing/promoting any event."""
import collections,itertools,json,pathlib,statistics,sys
import pysam
from site_recovery.core import cigar_gaps,identity,aligned_slice,usable_read,Config
root=pathlib.Path(sys.argv[1]);audit=json.loads((root/'audit/candidate_audit.json').read_text());cases=[]
selected={'Sniffles2.INS.3SB9','Sniffles2.INS.FSB9','Sniffles2.INS.10SBD'}
with pysam.AlignmentFile(str(root/'calls/alignments.bam'),'rb') as bam,pysam.FastaFile(str(root/'reference/marine_mag_pool.fasta')) as ref:
 for d in audit:
  c=d['candidate']
  if c['id'] not in selected:continue
  targets={(g['read'],g['start'],g['end']) for g in d['passing_gaps']};seqs={};states=[]
  for r in bam.fetch(c['contig'],c['start0']-1000,c['start0']+1000):
   if not usable_read(r,Config()):continue
   for g in cigar_gaps(r):
    if (r.query_name,g.start0,g.end0) in targets:seqs[r.query_name]=g.sequence
   pos=c['start0']+1;s=aligned_slice(r,pos-100,pos+100)
   if s and identity(s,ref.fetch(c['contig'],pos-100,pos+100).upper())>=.95:states.append(r.query_name)
  pairs=[dict(read_a=a,read_b=b,length_a=len(seqs[a]),length_b=len(seqs[b]),global_identity=identity(seqs[a],seqs[b])) for a,b in itertools.combinations(sorted(seqs),2)]
  cases.append(dict(candidate=c,passing_gaps=d['passing_gaps'],pairwise_insert_identity=pairs,locally_reference_like_read_names=sorted(set(states)),caution='Diagnostic only; no canonical breakpoint normalization, complete allele reconstruction or IS classification performed.'))
split_summary=[]
for d in audit:
 c=d['candidate']
 if c['svlen']<1000 or c['status']!='unresolved' or c.get('split_geometry_compatible_names',0)<3:continue
 opposite={};net={}
 for g in d['split_geometry']:
  opposite.setdefault(g['read'],[]).append(g['query_gap'] if c['svtype']=='DEL' else g['reference_gap']);net.setdefault(g['read'],[]).append(g['estimated_svlen'])
 vals=[min(v,key=abs) for v in opposite.values()]
 split_summary.append(dict(id=c['id'],svtype=c['svtype'],svlen=c['svlen'],split_read_names=len(vals),reads_opposite_gap_abs_le10=sum(abs(v)<=10 for v in vals),median_opposite_unaligned_gap=statistics.median(vals),note='Unaligned opposite gaps require local realignment; they do not by themselves prove replacement.'))
(root/'audit/example_cases.json').write_text(json.dumps(cases,indent=2)+'\n');(root/'audit/split_gap_summary.json').write_text(json.dumps(split_summary,indent=2)+'\n')
print('Unresolved long split geometry >=3:',len(split_summary),'of which >=3 names with opposite gap <=10:',sum(d['reads_opposite_gap_abs_le10']>=3 for d in split_summary))
for c in cases:
 print(c['candidate']['id'],'ref-like',len(c['locally_reference_like_read_names']),'gaps',[(g['start'],g['length']) for g in c['passing_gaps']]);print([(p['length_a'],p['length_b'],round(p['global_identity'],4)) for p in c['pairwise_insert_identity']])
