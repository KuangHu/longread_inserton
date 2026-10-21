"""Read-only audit of all original INS/DEL calls; never changes thresholds/results."""
import collections,csv,json,pathlib,sys
import pysam
from site_recovery.core import Config,candidates_from_vcf,cigar_gaps,usable_read,aligned_slice,identity
root=pathlib.Path(sys.argv[1]);out=root/'audit';out.mkdir(exist_ok=True);cfg=Config()
old={x['source_call']['source_id']:x for x in map(json.loads,(root/'recovered/records.jsonl').open())}
candidates=[c for c in candidates_from_vcf(root/'targeted/calls.vcf') if c.svtype in ('INS','DEL')]
infos={};needed=set()
with pysam.VariantFile(str(root/'targeted/calls.vcf')) as v:
 for x in v:
  if x.id not in old or x.info['SVTYPE'] not in ('INS','DEL'):continue
  names=set(x.info.get('RNAMES',()));needed.update(names)
  infos[x.id]=dict(names=names,support=x.info.get('SUPPORT'),stdev_pos=x.info.get('STDEV_POS'),stdev_len=x.info.get('STDEV_LEN'),alt_length=len(x.alts[0])-len(x.ref) if x.alts and not x.alts[0].startswith('<') else None)
by_name=collections.defaultdict(list)
with pysam.AlignmentFile(str(root/'calls/alignments.bam'),'rb') as b:
 for r in b.fetch(until_eof=True):
  if r.query_name in needed:by_name[r.query_name].append(r)
rows=[];details=[]
with pysam.AlignmentFile(str(root/'calls/alignments.bam'),'rb') as bam,pysam.FastaFile(str(root/'reference/marine_mag_pool.fasta')) as ref:
 for c in candidates:
  n=ref.get_reference_length(c.chrom);op=1 if c.svtype=='INS' else 2;slop=max(5,round(c.svlen*.25));info=infos[c.source_id]
  sets=collections.defaultdict(set);groups=collections.defaultdict(set);anchor_examples=[];match_gaps=[]
  for r in bam.fetch(c.chrom,max(0,c.start0-600),min(n,c.start0+c.svlen+601 if op==2 else c.start0+601)):
   sets['local_names'].add(r.query_name)
   for g in cigar_gaps(r):
    length=g.qend-g.qstart if op==1 else g.end0-g.start0
    if g.op!=op or abs(length-c.svlen)>slop or abs(g.start0-c.start0)>500:continue
    sets['compatible_gap_any_alignment_500'].add(r.query_name)
    if not usable_read(r,cfg):continue
    sets['compatible_gap_primary_mq20_500'].add(r.query_name)
    if abs(g.start0-c.start0)>100:continue
    sets['compatible_gap_primary_mq20_100'].add(r.query_name)
    if g.start0<100 or g.end0+100>n:sets['anchor_contig_edge'].add(r.query_name);continue
    l=aligned_slice(r,g.start0-100,g.start0);rr=aligned_slice(r,g.end0,g.end0+100)
    if l is None or rr is None:
     sets['anchor_structural_failure'].add(r.query_name)
     if len(anchor_examples)<3:anchor_examples.append(dict(read=r.query_name,cigar=r.cigarstring,problem='anchor_structure',start=g.start0,end=g.end0))
     continue
    li=identity(l,ref.fetch(c.chrom,g.start0-100,g.start0).upper());ri=identity(rr,ref.fetch(c.chrom,g.end0,g.end0+100).upper())
    sets['anchor_structurally_spanned'].add(r.query_name)
    if min(li,ri)<.95:
     sets['anchor_identity_below_95'].add(r.query_name)
     if len(anchor_examples)<3:anchor_examples.append(dict(read=r.query_name,problem='anchor_identity',left_identity=li,right_identity=ri,start=g.start0,end=g.end0))
    if min(li,ri)>=.90:sets['anchor_pass_90_diagnostic_only'].add(r.query_name)
    if min(li,ri)>=.95:
     sets['anchor_pass_95'].add(r.query_name);groups[(g.start0,g.end0)].add(r.query_name)
     match_gaps.append(dict(read=r.query_name,start=g.start0,end=g.end0,length=length))
  split=[]
  for name in info['names']:
   rr=by_name[name]
   if any(r.has_tag('SA') for r in rr):sets['caller_support_names_with_SA'].add(name)
   if any(r.is_supplementary for r in rr):sets['caller_support_names_with_supplementary'].add(name)
   if any(r.reference_name==c.chrom and usable_read(r,cfg) for r in rr):sets['caller_support_names_usable_primary_here'].add(name)
   if not rr:sets['caller_support_names_missing_in_bam'].add(name)
   segs=[]
   for r in rr:
    if r.is_unmapped or r.is_secondary or r.reference_name!=c.chrom or r.mapping_quality<20:continue
    cig=r.cigartuples or [];hardleft=cig[0][1] if cig and cig[0][0]==5 else 0
    segs.append((r.reference_start,r.reference_end,r.query_alignment_start+hardleft,r.query_alignment_end+hardleft,r.is_reverse))
   for a in segs:
    for b in segs:
     if a[0]>=b[0] or a[4]!=b[4] or abs(a[1]-c.start0)>500:continue
     rg=b[0]-a[1];qg=b[2]-a[3];est=qg-rg if op==1 else rg-qg
     if abs(est-c.svlen)<=slop and min(rg,qg)>=-100 and a[1]-a[0]>=100 and b[1]-b[0]>=100:
      sets['split_geometry_compatible_names'].add(name)
      split.append(dict(read=name,left=list(a),right=list(b),query_gap=qg,reference_gap=rg,estimated_svlen=est))
  rec=old[c.source_id];row=dict(id=c.source_id,svtype=c.svtype,svlen=c.svlen,contig=c.chrom,start0=c.start0,end0=c.end0,status=rec['status'],reason=';'.join(rec['reasons']),caller_support=info['support'],caller_stdev_pos=info['stdev_pos'],caller_stdev_len=info['stdev_len'],caller_alt_length=info['alt_length'],best_exact_group=max(map(len,groups.values()),default=0),exact_groups=len(groups))
  row.update({key:len(value) for key,value in sets.items()});rows.append(row)
  details.append(dict(candidate=row,names_by_stage={k:sorted(v) for k,v in sets.items()},anchor_examples=anchor_examples,passing_gaps=match_gaps,split_geometry=split))
fields=list(dict.fromkeys(k for r in rows for k in r))
with (out/'candidate_audit.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(rows)
(out/'candidate_audit.json').write_text(json.dumps(details,indent=2)+'\n')
summary={}
for label,lo,hi in [('50-499',50,500),('500-999',500,1000),('1000-1999',1000,2000),('2000-4999',2000,5000),('5000-9999',5000,10000)]:
 rr=[r for r in rows if lo<=r['svlen']<hi]
 summary[label]=dict(candidates=len(rr),statuses=dict(collections.Counter(r['status'] for r in rr)),any_compatible_primary_cigar=sum(r.get('compatible_gap_primary_mq20_100',0)>0 for r in rr),three_compatible_primary_cigar=sum(r.get('compatible_gap_primary_mq20_100',0)>=3 for r in rr),three_anchored95=sum(r.get('anchor_pass_95',0)>=3 for r in rr),three_split_geometry=sum(r.get('split_geometry_compatible_names',0)>=3 for r in rr))
(out/'audit_summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
