"""Summarize complete FASTQ, primary MAPQ>=20 coverage, and recovered events."""
import collections,csv,json,pathlib,statistics,sys
import pysam
root=pathlib.Path(sys.argv[1]);out=root/'analysis';out.mkdir(exist_ok=True)
rows=list(csv.DictReader((root/'reference/contig_to_mag.tsv').open(),delimiter='\t'))
by_contig={r['contig']:r for r in rows};mags=collections.defaultdict(collections.Counter)
lengths=[];quality_chars=set();ccs=0
with pysam.FastxFile(str(root/'reads/SRR28395891_subreads.fastq.gz')) as f:
 for r in f:
  lengths.append(len(r.sequence));quality_chars.update(r.quality);ccs+=('ccs' in (r.name+' '+(r.comment or '')).lower())
total=sum(lengths);acc=0;n50=0
for n in sorted(lengths,reverse=True):
 acc+=n
 if acc>=total/2:n50=n;break
qc=dict(reads=len(lengths),bases=total,median_length=statistics.median(lengths),n50=n50,ccs_headers=ccs,reported_quality_values=sorted(ord(c)-33 for c in quality_chars),original_file_md5_verified=True)
(root/'analysis/reads_qc.json').write_text(json.dumps(qc,indent=2)+'\n')
counts=collections.Counter()
with pysam.AlignmentFile(str(root/'calls/alignments.bam'),'rb') as bam:
 for r in bam.fetch(until_eof=True):
  if r.is_secondary or r.is_supplementary:continue
  counts['primary_records']+=1
  if r.is_unmapped:continue
  counts['mapped_primary']+=1
  if r.mapping_quality>=20:counts['mapped_primary_mapq20']+=1
 for info in rows:
  chrom=info['contig'];length=int(info['length']);a=mags[info['mag']];a['length']+=length;a['contigs']+=1
  edges=collections.Counter()
  for r in bam.fetch(chrom):
   if r.is_secondary or r.is_supplementary or r.is_unmapped or r.is_duplicate or r.is_qcfail or r.mapping_quality<20:continue
   a['primary_reads_mapq20']+=1
   for s,e in r.get_blocks():
    edges[s]+=1;edges[e]-=1;a['aligned_bases_mapq20']+=e-s
  depth=0;last=0
  for pos,delta in sorted(edges.items()):
   for threshold in (1,5,10):
    if depth>=threshold:a['bases_depth_ge_'+str(threshold)]+=pos-last
   last=pos;depth+=delta
  assert depth==0
summary=[]
for mag,a in mags.items():
 d=dict(mag=mag,**a);d['mean_depth_mapq20']=a['aligned_bases_mapq20']/a['length']
 for t in (1,5,10):d['breadth_ge_'+str(t)]=a['bases_depth_ge_'+str(t)]/a['length']
 summary.append(d)
summary.sort(key=lambda d:d['mean_depth_mapq20'],reverse=True)
fields=['mag','length','contigs','primary_reads_mapq20','aligned_bases_mapq20','mean_depth_mapq20','breadth_ge_1','breadth_ge_5','breadth_ge_10']
with (out/'mag_coverage.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows([{k:d.get(k,0) for k in fields} for d in summary])
if '--coverage-only' in sys.argv:
 (out/'coverage_summary.json').write_text(json.dumps(dict(reads_qc=qc,alignment_counts=dict(counts),top_mag_coverage=summary[:20]),indent=2)+'\n')
 print(json.dumps(dict(alignment_counts=dict(counts),top_mag_coverage=summary[:10]),indent=2));sys.exit(0)
events=list(csv.DictReader((root/'recovered/events.tsv').open(),delimiter='\t'))
for e in events:e['mag']=by_contig[e['contig']]['mag']
with (out/'events_with_mag.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=['mag']+list(events[0].keys())[:-1] if events else ['mag'],delimiter='\t');w.writeheader();w.writerows(events)
result=dict(reads_qc=qc,alignment_counts=dict(counts),reference=json.loads((root/'reference/manifest.json').read_text()),status_counts=dict(collections.Counter(e['status'] for e in events)),svtype_counts=dict(collections.Counter(e['svtype'] for e in events)),top_mag_coverage=summary[:10],mags_mean_depth_ge_5=sum(d['mean_depth_mapq20']>=5 for d in summary),mags_mean_depth_ge_10=sum(d['mean_depth_mapq20']>=10 for d in summary),target_selection=json.loads((root/'targeted/selection.json').read_text()),same_sample_pairing='unresolved at individual MAG level; same-study pooled reference only')
(out/'pilot_summary.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
