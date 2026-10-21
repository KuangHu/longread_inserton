"""Fetch published study reference data; validate Zenodo file MD5.
Needs openpyxl only to read the authors' supplement. No raw reads are downloaded here.
Usage: python fetch_reference.py DATA_DIRECTORY
"""
import hashlib,io,json,pathlib,sys,urllib.request,zipfile
import openpyxl
root=pathlib.Path(sys.argv[1]);d=root/'metadata';d.mkdir(parents=True,exist_ok=True)
def fetch(url,path,md5=None,size=None):
 if path.exists():
  if md5:
   h=hashlib.md5()
   with path.open('rb') as f:
    for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
   if h.hexdigest()!=md5:raise ValueError(f'Existing file checksum mismatch: {path}')
  return
 p=path.with_suffix(path.suffix+'.part');h=hashlib.md5();n=0
 with urllib.request.urlopen(url,timeout=60) as r,p.open('wb') as f:
  while b:=r.read(1024*1024):
   f.write(b);h.update(b);n+=len(b)
 if (md5 and h.hexdigest()!=md5) or (size and n!=size):raise ValueError('Download size/checksum mismatch')
 p.replace(path)
fetch('https://zenodo.org/api/records/14650439',d/'zenodo.json')
m=json.loads((d/'zenodo.json').read_text());entry=next(f for f in m['files'] if f['key']=='MAG-LRMED.tar')
fetch(entry['links']['self'],d/'MAG-LRMED.tar',entry['checksum'].removeprefix('md5:'),entry['size'])
fetch('https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12379486/supplementaryFiles',d/'supplements.zip')
with zipfile.ZipFile(d/'supplements.zip') as z:
 w=openpyxl.load_workbook(io.BytesIO(z.read(next(n for n in z.namelist() if n.endswith('.xlsx')))),read_only=True,data_only=True)
 rows=[r for r in list(w['Table_S2'].values)[3:] if r[1]]
 rank={r[0]:r[1] for r in list(w['Table_S3'].values)[3:] if isinstance(r[1],(int,float))}
 a=sorted([dict(mag=r[0],accession=r[1],length=r[3],contigs=r[4],completeness=r[7],contamination=r[8],abundance=rank.get(r[0],0)) for r in rows],key=lambda x:x['abundance'],reverse=True)
 assert len(a)==418
 (d/'mag_rank.json').write_text(json.dumps(a,indent=2)+'\n')
print('Reference download verified; per-MAG sample provenance remains unresolved.')
