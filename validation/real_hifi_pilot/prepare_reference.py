"""Build a namespaced competitive reference from the verified Zenodo MAG archive.
Usage: python prepare_reference.py /absolute/path/to/hifi_marine_pilot
The reference is study-wide; this does not establish per-MAG sample provenance.
"""
import csv,hashlib,json,pathlib,sys,tarfile
root=pathlib.Path(sys.argv[1]);out=root/'reference';out.mkdir(exist_ok=True)
p=root/'metadata/MAG-LRMED.tar'
h=hashlib.md5()
with p.open('rb') as f:
 for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
assert h.hexdigest()=='51d2fd500477514ce776374983b974c1'
allowed={r['mag'] for r in json.loads((root/'metadata/mag_rank.json').read_text())}
rows=[];seen=set();excluded=[]
with tarfile.open(p) as t,(out/'marine_mag_pool.fasta').open('w') as f:
 for m in t:
  if not m.isfile() or not m.name.endswith('.fasta'):continue
  mag=pathlib.Path(m.name).stem
  if mag not in allowed:excluded.append(mag);continue
  name=None;seq=[]
  def save():
   if name is None:return
   s=''.join(seq).upper();ident=mag+'__'+name
   assert ident not in seen and s and set(s)<=set('ACGTRYSWKMBDHVN')
   seen.add(ident);f.write('>'+ident+'\n'+s+'\n');rows.append((ident,mag,name,len(s)))
  for line in t.extractfile(m):
   line=line.decode().strip()
   if line.startswith('>'):save();name=line[1:].split()[0];seq=[]
   else:seq.append(line)
  save()
with (out/'contig_to_mag.tsv').open('w') as f:
 w=csv.writer(f,delimiter='\t');w.writerow(['contig','mag','original_contig','length']);w.writerows(rows)
h=hashlib.sha256()
with (out/'marine_mag_pool.fasta').open('rb') as f:
 for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
meta=dict(mags=len(set(r[1] for r in rows)),contigs=len(rows),bases=sum(r[3] for r in rows),excluded_not_in_paper_table=excluded,reference_scope='study-wide MAG pool; exact per-MAG sample provenance unresolved',zenodo_record='14650439',archive_md5='51d2fd500477514ce776374983b974c1',reference_sha256=h.hexdigest())
(out/'manifest.json').write_text(json.dumps(meta,indent=2)+'\n');print(meta)
from collections import Counter
lengths=Counter();counts=Counter()
for _,mag,_,length in rows:lengths[mag]+=length;counts[mag]+=1
mismatch=[dict(mag=r['mag'],archive_length=lengths[r['mag']],table_length=r['length'],archive_contigs=counts[r['mag']],table_contigs=r['contigs']) for r in json.loads((root/'metadata/mag_rank.json').read_text()) if lengths[r['mag']]!=r['length'] or counts[r['mag']]!=r['contigs']]
(out/'table_validation.json').write_text(json.dumps(dict(compared_mags=len(lengths),mismatches=mismatch),indent=2)+'\n')
