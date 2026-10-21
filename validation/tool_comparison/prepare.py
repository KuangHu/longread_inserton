"""Prepare frozen target alignments and a labelled, focused diagnostic read set."""
import argparse, gzip, hashlib, json
from pathlib import Path
import pysam

p = argparse.ArgumentParser()
p.add_argument('pilot', type=Path)
p.add_argument('out', type=Path)
a = p.parse_args()
a.out.mkdir(parents=True, exist_ok=True)
selection = json.loads((a.pilot/'targeted/selection.json').read_text())
targets = set(selection['contigs'])
src = pysam.AlignmentFile(str(a.pilot/'calls/alignments.bam'))
ref = pysam.FastaFile(str(a.pilot/'reference/marine_mag_pool.fasta'))
header = src.header.to_dict()
header['SQ'] = [s for s in header['SQ'] if s['SN'] in targets]
newheader = pysam.AlignmentHeader.from_dict(header)
focal = next(n for n in targets if n.startswith('Prochlorococcus-A-LRMED-150__') and n.endswith('-C11'))
focal_names = set()
target_names = set()
records = 0
with pysam.AlignmentFile(str(a.out/'targets.bam'), 'wb', header=newheader) as dst:
    for sq in header['SQ']:
        n = sq['SN']
        for read in src.fetch(n):
            dst.write(pysam.AlignedSegment.fromstring(read.to_string(), newheader))
            records += 1
            if not read.is_secondary and not read.is_supplementary and read.mapping_quality >= 20:
                target_names.add(read.query_name)
                if n == focal and read.reference_start < 110000 and read.reference_end > 80000:
                    focal_names.add(read.query_name)
pysam.index(str(a.out/'targets.bam'))
with (a.out/'targets.fasta').open('w') as f:
    for s in header['SQ']:
        f.write('>'+s['SN']+'\n'+ref.fetch(s['SN'])+'\n')
pysam.faidx(str(a.out/'targets.fasta'))
(a.out/'contigs.txt').write_text(''.join(s['SN']+'\n' for s in header['SQ']))
counts = {k: {'reads':0,'bases':0} for k in ('targets','focal')}
with gzip.open(a.out/'targets.fastq.gz','wt') as fo, gzip.open(a.out/'focal.fastq.gz','wt') as ff, (a.out/'targets_reads.fasta').open('w') as fa, (a.out/'focal_reads.fasta').open('w') as ffa:
    with pysam.FastxFile(str(a.pilot/'reads/SRR28395891_subreads.fastq.gz')) as fq:
        for r in fq:
            if r.name in target_names:
                fo.write(f'@{r.name}\n{r.sequence}\n+\n{r.quality}\n')
                fa.write(f'>{r.name}\n{r.sequence}\n')
                counts['targets']['reads'] += 1
                counts['targets']['bases'] += len(r.sequence)
                if r.name in focal_names:
                    ff.write(f'@{r.name}\n{r.sequence}\n+\n{r.quality}\n')
                    ffa.write(f'>{r.name}\n{r.sequence}\n')
                    counts['focal']['reads'] += 1
                    counts['focal']['bases'] += len(r.sequence)
manifest = {'source_run':'SRR28395891','target_contigs':len(targets),'target_bases':sum(s['LN'] for s in header['SQ']), 'bam_records':records, 'read_counts':counts,
 'focal':{'contig':focal,'start0':80000,'end0':110000,'purpose':'diagnostic region around previously audited 1.95-kb insertion; NOT an unbiased discovery benchmark; whole original reads retained'},
 'alignment_policy':'Original competitive alignments/MAPQ against all 418 MAGs preserved; target header subset only. FASTQ recruitment uses primary MAPQ>=20 reads.',
 'provenance':'Study-wide MAG pool; not proven same-sample per MAG. No lab strain added.'}
manifest['sha256'] = {n:hashlib.sha256((a.out/n).read_bytes()).hexdigest() for n in ['targets.fasta','targets.bam','targets.fastq.gz','focal.fastq.gz']}
(a.out/'input_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(json.dumps(manifest,indent=2))
