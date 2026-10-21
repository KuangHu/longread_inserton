"""Reproduce focal-anchor and IR-interval checks without promoting them to IS."""
import argparse,csv,json
from pathlib import Path
import pysam
from site_recovery.core import usable_read,Config,aligned_slice,identity
p=argparse.ArgumentParser();p.add_argument('root',type=Path);p.add_argument('pilot',type=Path);a=p.parse_args()
rc=lambda s:s.translate(str.maketrans('ACGT','TGCA'))[::-1]
f=pysam.FastaFile(str(a.root/'inputs/targets.fasta'));b=pysam.AlignmentFile(str(a.root/'inputs/targets.bam'))
out=[]
for c in csv.DictReader((a.root/'runs/palidis_ir_mapping/structural_candidates.tsv').open(),delimiter='\t'):
    chrom=c['contig'];start,end=int(c['start0']),int(c['end0']);n=int(c['max_exact_ir_length'])
    seq=f.fetch(chrom,start,end);names=[]
    if start>=100 and end+100<=f.get_reference_length(chrom):
        target=f.fetch(chrom,start-100,end+100)
        for r in b.fetch(chrom,start-100,end+100):
            if not usable_read(r,Config()):continue
            q=aligned_slice(r,start-100,end+100)
            if q and identity(q,target)>=.95:names.append(r.query_name)
    out.append(dict(candidate=c['id'],contig=chrom,start0=start,end0=end,terminal_repeat_length=n,exact_inverted_terminals=seq[:n]==rc(seq[-n:]),occupied_spanning_read_names=sorted(set(names)),is_classification='unclassified; no insertion/empty polymorphism required or demonstrated'))
(a.root/'runs/palidis_ir_mapping/read_check.json').write_text(json.dumps(out,indent=2)+'\n')
anchors=list(pysam.FastxFile(str(a.root/'runs/dive_focal_compat/output/focal_assemb_anchors.fasta')))
diag=[]
for v in pysam.VariantFile(str(a.pilot/'targeted/calls.vcf')):
    if v.id!='Sniffles2.INS.3SB9':continue
    ins=v.alts[0][len(v.ref):]
    for anchor in anchors:
        s=anchor.sequence
        diag.append(dict(anchor_id=anchor.name,anchor_length=len(s),focal_insertion_length=len(ins),matches=[dict(strand=d,offset0=ins.find(q)) for d,q in [('+',s),('-',rc(s))] if q in ins],meaning='Anchor within a previously nominated insertion does not recover its ends or full sequence. Focal reads selected around this known locus.'))
(a.root/'runs/dive_focal_compat/output/anchor_diagnostic.json').write_text(json.dumps(diag,indent=2)+'\n')
print('IR intervals checked:',len(out),'DIVE anchors checked:',len(diag))
