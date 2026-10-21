"""Exploratory IR pairing; NOT the full Palidis/Palidis4 pipeline or IS validation.

Use native pal-MEM + native clip_reads + native Bowtie2 mapping, then require
unclipped, ungapped exact matches. Pair reverse-complement repeats on a contig.
No protein annotation is fabricated. These are structural candidates only.
"""
import argparse, collections, csv, json
from pathlib import Path
import pysam

p=argparse.ArgumentParser();p.add_argument('root',type=Path);a=p.parse_args()
root=a.root;out=root/'runs/palidis_ir_mapping'
aliases={v:k for k,v in json.loads((root/'inputs/alias_map.json').read_text()).items()}
rc=lambda s:s.translate(str.maketrans('ACGTN','TGCAN'))[::-1]
positions=collections.defaultdict(lambda:[set(),set()]);stats=collections.Counter()
with pysam.AlignmentFile(str(out/'irs.sam')) as sam:
    for r in sam:
        stats['sam_records']+=1
        if r.is_unmapped:continue
        stats['mapped_records']+=1
        if r.is_secondary or r.is_supplementary:continue
        if r.get_tag('NM')!=0 or any(op not in (0,7) for op,n in r.cigartuples):continue
        stats['exact_full_length_records']+=1
        s=r.query_sequence.upper();rev=rc(s)
        if s==rev:continue
        key=min(s,rev);side=int(s!=key)
        positions[(r.reference_name,key)][side].add((r.reference_start,r.reference_end))
ref=pysam.FastaFile(str(root/'inputs/aliases.fasta'));pairs={}
for (chrom,seq),(left,right) in positions.items():
    for a in left:
        for b in right:
            lo,hi=sorted([a,b]);size=hi[1]-lo[0]
            if 500<=size<=3000 and lo[1]<hi[0]:
                pairs.setdefault((chrom,lo[0],hi[1]),[]).append(len(seq))
rows=[]
with (out/'structural_candidates.fasta').open('w') as f:
    for i,((chrom,start,end),lengths) in enumerate(sorted(pairs.items()),1):
        ident=f'IRPAIR_{i}'
        rows.append(dict(id=ident,contig=aliases[chrom],start0=start,end0=end,length=end-start,max_exact_ir_length=max(lengths),status='IR_pair_only_not_validated_IS'))
        for tag,s in [('candidate',ref.fetch(chrom,start,end)),('upstream',ref.fetch(chrom,max(0,start-500),start)),('downstream',ref.fetch(chrom,end,min(ref.get_reference_length(chrom),end+500)))]:
            f.write(f'>{ident}|{tag}\n{s}\n')
with (out/'structural_candidates.tsv').open('w') as f:
    fields=['id','contig','start0','end0','length','max_exact_ir_length','status']
    w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(rows)
stats['distinct_structural_intervals']=len(rows)
stats['overlapping_intervals_not_deduplicated']=True
stats['full_palidis_pipeline_completed']=False
stats['classification']='Exploratory structural IR pairs, not validated insertion events or IS elements'
(out/'summary.json').write_text(json.dumps(stats,indent=2)+'\n');print(json.dumps(stats,indent=2))
