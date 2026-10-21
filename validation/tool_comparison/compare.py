"""Compare native calls and audit full-CIGAR evidence; no ground-truth claims.

The audit is diagnostic, not a new insertion/IS caller. Split-only evidence is
not evaluated. Empty-state evidence is reported but never required for INS.
"""
import argparse, collections, csv, itertools, json, re
from pathlib import Path
import pysam
from site_recovery.core import Config, aligned_slice, anchored_gap, cigar_gaps, identity, usable_read

p=argparse.ArgumentParser();p.add_argument('root',type=Path);p.add_argument('pilot',type=Path);p.add_argument('out',type=Path);a=p.parse_args()
a.out.mkdir(parents=True,exist_ok=True)
aliases={v:k for k,v in json.loads((a.root/'inputs/alias_map.json').read_text()).items()}
ref=pysam.FastaFile(str(a.root/'inputs/targets.fasta'))
bam=pysam.AlignmentFile(str(a.root/'inputs/targets.bam'))
cfg=Config();calls=[]

def add_vcf(label,path,sniffles=False):
    for r in pysam.VariantFile(str(path)):
        info=dict(r.info)
        typ=info.get('SVTYPE');length=abs(info.get('SVLEN',0))
        if typ not in ('INS','DEL') or not 500<=length<=50000:continue
        chrom=aliases.get(r.contig,r.contig)
        # Sniffles uses a VCF anchor base; DeBreak's writer emits unpadded alleles.
        start=r.start+int(sniffles);end=start+(length if typ=='DEL' else 0)
        seq=''
        if typ=='INS' and r.alts and set(r.alts[0].upper())<=set('ACGTN'):
            seq=r.alts[0][len(r.ref):] if sniffles else r.alts[0]
        if typ=='DEL':seq=ref.fetch(chrom,start,end)
        calls.append(dict(tool=label,id=r.id,contig=chrom,start0=start,end0=end,svtype=typ,length=length,
            native_filter=';'.join(r.filter),native_precise=bool(info.get('PRECISE')),
            native_support=info.get('SUPPORT',info.get('SUPPREAD',0)),sequence=seq.upper(),
            upstream=ref.fetch(chrom,max(0,start-500),start),downstream=ref.fetch(chrom,end,min(ref.get_reference_length(chrom),end+500)),
            flank_source='reference_derived_context',native_empty_reads=None))

add_vcf('Sniffles2',a.pilot/'targeted/calls.vcf',True)
for label,name in [('DeBreak_default','debreak'),('DeBreak_no_depth_cap','debreak_no_depth_cap')]:
    add_vcf(label,a.root/f'runs/{name}/output/debreak.vcf')
tldr=a.root/'runs/tldr_bounded/calls.table.txt'
if not tldr.exists():tldr=a.root/'runs/tldr_no_detail/calls.table.txt'
for label,table in [('TLDR',tldr),('TLDR_focal',a.root/'runs/tldr_focal/calls.table.txt')]:
    if not table.exists():continue
    for r in csv.DictReader(table.open(),delimiter='\t'):
        if not 500<=int(r['LengthIns'])<=50000:continue
        cons=r['Consensus'];segs=sorted(re.finditer('[a-z]+',cons),key=lambda m:len(m.group()),reverse=True)
        m=segs[0] if segs else None
        seq=m.group().upper() if m else '';up=cons[:m.start()].upper() if m else '';dn=cons[m.end():].upper() if m else ''
        calls.append(dict(tool=label,id=r['UUID'],contig=aliases[r['Chrom']],start0=int(r['Start']),end0=int(r['End']),svtype='INS',length=int(r['LengthIns']),native_filter=r['Filter'],native_precise=None,native_support=int(r['UsedReads']),sequence=seq,upstream=up,downstream=dn,flank_source='occupied_consensus',native_empty_reads=r['EmptyReads'],native_span_reads=r['SpanReads'],consensus_partition_length_ok=len(seq)==int(r['LengthIns'])))

details=[]
for c in calls:
    chrom,start,end=c['contig'],c['start0'],c['end0'];n=ref.get_reference_length(chrom)
    gap_hits={};occupied={};empty={};geometric=set()
    for r in bam.fetch(chrom,max(0,start-600),min(n,max(start,end)+601)):
        if not usable_read(r,cfg):continue
        for g in cigar_gaps(r):
            if g.op!=(1 if c['svtype']=='INS' else 2):continue
            size=(g.qend-g.qstart) if g.op==1 else g.end0-g.start0
            if abs(g.start0-start)>100 or abs(size-c['length'])>max(5,.25*c['length']):continue
            if c['svtype']=='DEL' and abs(g.end0-end)>100:continue
            geometric.add(r.query_name)
            if not anchored_gap(r,g,ref,cfg):continue
            si=identity(c['sequence'],g.sequence) if g.op==1 and c['sequence'] else None
            info=dict(read=r.query_name,start0=g.start0,end0=g.end0,length=size,insert_identity=si)
            gap_hits[r.query_name]=info
            if g.op==1 and si is not None and si>=.95:occupied[r.query_name]=info
            elif g.op==2:empty[r.query_name]=info
        if c['svtype']=='INS':
            lo,hi=start-100,start+100
            if 0<=lo<hi<=n:
                s=aligned_slice(r,lo,hi)
                if s and identity(s,ref.fetch(chrom,lo,hi).upper())>=.95:empty[r.query_name]=dict(read=r.query_name)
        else:
            lo,hi=start-100,end+100
            if 0<=lo<hi<=n:
                s=aligned_slice(r,lo,hi)
                if s and identity(s,ref.fetch(chrom,lo,hi).upper())>=.95:occupied[r.query_name]=dict(read=r.query_name)
    # Compatible read-to-reference boundaries are evidence, not a single-base truth.
    starts=[x['start0'] for x in occupied.values() if 'start0' in x]
    supported=len(occupied)>=3 and len(c['sequence'])>=500
    c.update(sequence_length=len(c['sequence']),geometric_gap_reads=len(geometric),anchored_gap_reads=len(gap_hits),sequence_supported_occupied_reads=len(occupied),reference_like_or_empty_gap_reads=len(empty),occupied_start_min0=min(starts) if starts else None,occupied_start_max0=max(starts) if starts else None,audit_three_occupied_reads=supported,audit_three_occupied_and_gap_reads=supported and len(gap_hits)>=3,
        audit_label='occupied_sequence_supported_not_IS_classified' if supported else 'candidate_not_promoted')
    details.append(dict(tool=c['tool'],id=c['id'],anchored_gaps=list(gap_hits.values()),occupied_support=list(occupied.values()),empty_support=list(empty.values())))

summary={}
for label in dict.fromkeys(c['tool'] for c in calls):
    rr=[c for c in calls if c['tool']==label]
    summary[label]={}
    for typ in ('INS','DEL'):
        ss=[c for c in rr if c['svtype']==typ]
        summary[label][typ]=dict(candidates=len(ss),native_pass=sum(c['native_filter']=='PASS' for c in ss),sequence_present=sum(bool(c['sequence']) for c in ss),three_occupied_read_support=sum(c['audit_three_occupied_reads'] for c in ss),three_occupied_and_gap_support=sum(c['audit_three_occupied_and_gap_reads'] for c in ss),native_pass_and_three_occupied=sum(c['native_filter']=='PASS' and c['audit_three_occupied_reads'] for c in ss),length_bins={f'{lo}-{hi-1}':sum(lo<=c['length']<hi for c in ss) for lo,hi in [(500,1000),(1000,2000),(2000,5000),(5000,50001)]})
overlaps=[]
labels=[k for k in summary if not k.endswith('_focal')]
for l1,l2 in itertools.combinations(labels,2):
    for typ in ('INS','DEL'):
        aa=[c for c in calls if c['tool']==l1 and c['svtype']==typ];bb=[c for c in calls if c['tool']==l2 and c['svtype']==typ]
        pairs=[(x,y) for x in aa for y in bb if x['contig']==y['contig'] and abs(x['start0']-y['start0'])<=100 and min(x['length'],y['length'])/max(x['length'],y['length'])>=.75]
        overlaps.append(dict(tool_a=l1,tool_b=l2,svtype=typ,total_a=len(aa),total_b=len(bb),matched_a=len({x['id'] for x,y in pairs}),matched_b=len({y['id'] for x,y in pairs}),matching_rule='same type/contig; start distance <=100 bp; length ratio >=0.75; NOT truth precision/recall'))
fields=list(dict.fromkeys(k for c in calls for k in c if k not in ('sequence','upstream','downstream')))
with (a.out/'candidates.tsv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(calls)
with (a.out/'candidate_sequences.fasta').open('w') as f:
    for c in calls:
        for tag,key in [('insertion','sequence'),('upstream','upstream'),('downstream','downstream')]:
            if c[key]:f.write(f'>{c["tool"]}|{c["id"]}|{tag}|{c["audit_label"]}\n{c[key]}\n')
with (a.out/'read_evidence.jsonl').open('w') as f:
    for d in details:f.write(json.dumps(d)+'\n')
(a.out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
(a.out/'overlap.json').write_text(json.dumps(overlaps,indent=2)+'\n')
print(json.dumps(summary,indent=2))
