"""Call the ten highest-coverage eligible MAGs using full-reference alignments.
Usage: python call_targets.py DATA_DIR SNIFFLES_BINARY RECOVERY_BINARY
Selection uses coverage and pre-existing assembly QC, never observed SV counts.
"""
import csv,hashlib,json,pathlib,subprocess,sys,time
import pysam
root=pathlib.Path(sys.argv[1]).resolve();out=root/'targeted';out.mkdir(exist_ok=False)
qc={r['mag']:r for r in json.loads((root/'metadata/mag_rank.json').read_text())}
bad={r['mag'] for r in json.loads((root/'reference/table_validation.json').read_text())['mismatches']}
coverage=list(csv.DictReader((root/'analysis/mag_coverage.tsv').open(),delimiter='\t'))
eligible=[r for r in coverage if r['mag'] not in bad and float(r['mean_depth_mapq20'])>=5 and float(r['breadth_ge_1'])>=0.7 and qc[r['mag']]['completeness']>=75 and qc[r['mag']]['contamination']<=5]
selected=sorted(eligible,key=lambda r:(-float(r['mean_depth_mapq20']),r['mag']))[:10];names={r['mag'] for r in selected}
assert names,'No eligible MAGs'
contigs=[r['contig'] for r in csv.DictReader((root/'reference/contig_to_mag.tsv').open(),delimiter='\t') if r['mag'] in names]
selection=dict(criteria='top 10 by MAPQ>=20 mean depth; depth>=5; breadth>=0.7; completeness>=75%; contamination<=5%; no archive/table discrepancy',selected=selected,eligible_mags=len(eligible),contigs=contigs)
(out/'selection.json').write_text(json.dumps(selection,indent=2)+'\n');print('Selected',len(names),'MAGs,',len(contigs),'contigs',flush=True)
bam=root/'calls/alignments.bam';ref=root/'reference/marine_mag_pool.fasta';vcf=out/'calls.vcf'
pysam.quickcheck(str(bam))
with pysam.AlignmentFile(str(bam),'rb') as b:assert b.has_index()
cmd=[sys.argv[2],'--input',str(bam),'--vcf',str(vcf),'--reference',str(ref),'--threads','4','--all-contigs','--mosaic','--mosaic-include-germline','--mosaic-af-min','0.05','--max-svlen-mosaic','100000','--minsupport','3','--minsvlen','50','--mapq','20','--output-rnames']
for c in contigs:cmd.extend(['--contig',c])
manifest=dict(status='running',started_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),command=cmd,reference_scope='study-wide MAG pool; per-MAG source sample unresolved',selection=selection['criteria'])
def save(): (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
save()
try:
 with (out/'sniffles.log').open('w') as f:subprocess.run(cmd,check=True,stdout=f,stderr=f)
 recover=[sys.argv[3],'recover','--reference',str(ref),'--bam',str(bam),'--vcf',str(vcf),'--platform','hifi','--sample','SRR28395891','--out',str(root/'recovered')]
 manifest['recover_command']=recover;save();subprocess.run(recover,check=True)
 manifest.update(status='complete',completed_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()));save()
except BaseException as e:
 manifest.update(status='failed',error=repr(e));save();raise
