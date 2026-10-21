"""Resource-bounded alignment stage for reproduction; requires fresh calls/ output."""
import hashlib,json,pathlib,subprocess,sys,time
import pysam
root=pathlib.Path(sys.argv[1]).resolve();tool=sys.argv[2];out=root/'calls'
if out.exists() and any(out.iterdir()):raise FileExistsError('calls/ must be new or empty')
out.mkdir(exist_ok=True)
ref=root/'reference/marine_mag_pool.fasta';reads=root/'reads/SRR28395891_subreads.fastq.gz'
source=json.loads(reads.with_suffix(reads.suffix+'.json').read_text())
def digest(p,algorithm='sha256'):
 h=hashlib.new(algorithm)
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
assert source['original_file_md5_verified'] and digest(reads,'md5')==source['source']['archive_md5']
cmd=[tool,'-K','50M','-a','-x','map-hifi','-t','4',str(ref),str(reads)]
m=dict(status='running',stage='alignment_only',command=cmd,inputs=[dict(path=str(p),sha256=digest(p)) for p in (ref,reads)],started_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()))
def save():(out/'manifest.json').write_text(json.dumps(m,indent=2)+'\n')
save()
try:
 sam=out/'alignments.unsorted.sam';bam=out/'alignments.bam'
 with sam.open('w') as f,(out/'minimap2.log').open('w') as log:subprocess.run(cmd,check=True,stdout=f,stderr=log)
 pysam.sort('-@','4','-o',str(bam),str(sam));pysam.index(str(bam));pysam.quickcheck(str(bam));sam.unlink()
 m.update(status='complete',completed_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()));save()
except BaseException as e:m.update(status='failed',error=repr(e));save();raise
