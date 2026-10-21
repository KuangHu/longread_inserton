"""Name-only aliases avoid upstream DeBreak's chr-prefix and '__' assumptions."""
import json, sys
from pathlib import Path
import pysam
p=Path(sys.argv[1])
b=pysam.AlignmentFile(str(p/'targets.bam'))
h=b.header.to_dict()
aliases={s['SN']:f'chr{i+1}' for i,s in enumerate(h['SQ'])}
for s in h['SQ']:s['SN']=aliases[s['SN']]
nh=pysam.AlignmentHeader.from_dict(h)
with pysam.AlignmentFile(str(p/'aliases.bam'),'wb',header=nh) as out:
    for r in b:
        fields=r.to_string().split('\t')
        fields[2]=aliases.get(fields[2],fields[2])
        fields[6]=aliases.get(fields[6],fields[6])
        for i,x in enumerate(fields[11:],11):
            if x.startswith('SA:Z:'):
                parts=[]
                for sa in x[5:].split(';'):
                    if not sa:continue
                    z=sa.split(',');z[0]=aliases.get(z[0],z[0]);parts.append(','.join(z))
                fields[i]='SA:Z:'+';'.join(parts)+';'
        out.write(pysam.AlignedSegment.fromstring('\t'.join(fields),nh))
pysam.index(str(p/'aliases.bam'))
ref=pysam.FastaFile(str(p/'targets.fasta'))
with (p/'aliases.fasta').open('w') as f:
    for old,new in aliases.items():f.write(f'>{new}\n{ref.fetch(old)}\n')
pysam.faidx(str(p/'aliases.fasta'))
(p/'alias_map.json').write_text(json.dumps(aliases,indent=2)+'\n')
print({'contig_aliases':len(aliases),'sequences_coordinates_mapq':'unchanged'})
