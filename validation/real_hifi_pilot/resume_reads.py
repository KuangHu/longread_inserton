"""Resume this pilot's interrupted FASTQ; publish only after full archive MD5.
This is pilot-specific, not a replacement for the generic download tool.
"""
import concurrent.futures,hashlib,json,pathlib,sys,time,urllib.request
root=pathlib.Path(sys.argv[1]);d=root/'reads';d.mkdir(exist_ok=True)
u='https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRR28395891&result=read_run&fields=run_accession,sample_accession,fastq_ftp,fastq_bytes,fastq_md5&format=json'
meta=json.load(urllib.request.urlopen(u,timeout=45))[0]
size=int(meta['fastq_bytes']);expected=meta['fastq_md5'];url='https://'+meta['fastq_ftp']
assert ';' not in url
p=d/'SRR28395891_subreads.fastq.gz';part=p.with_suffix(p.suffix+'.part')
if p.exists():raise FileExistsError(p)
start=part.stat().st_size if part.exists() else 0
assert 0<=start<=size
chunk=8*1024*1024;intervals=[(s,min(s+chunk,size)-1) for s in range(start,size,chunk)]
def get(interval):
 s,e=interval
 for trial in range(3):
  try:
   req=urllib.request.Request(url,headers={'Range':f'bytes={s}-{e}'})
   with urllib.request.urlopen(req,timeout=60) as r:
    assert r.status==206 and r.headers.get('Content-Range')==f'bytes {s}-{e}/{size}'
    b=r.read(e-s+2)
    assert len(b)==e-s+1
   return s,b
  except Exception:
   if trial==2:raise
with part.open('ab') as f,concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
 for s,b in pool.map(get,intervals):
  assert f.tell()==s
  f.write(b);f.flush();print('bytes',f.tell(),'/',size,flush=True)
h=hashlib.md5()
with part.open('rb') as f:
 for b in iter(lambda:f.read(chunk),b''):h.update(b)
assert part.stat().st_size==size and h.hexdigest()==expected,'Full checksum mismatch: .part retained'
result=dict(source=dict(run='SRR28395891',url=url,source_size=size,archive_md5=expected,ena_metadata=meta,metadata_url=u),original_file_md5_verified=True,bytes=size,resumed_at_byte=start,range_chunks=len(intervals),retrieved_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()))
part.replace(p);p.with_suffix(p.suffix+'.json').write_text(json.dumps(result,indent=2)+'\n');print('Full size and MD5 verified',flush=True)
