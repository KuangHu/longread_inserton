"""Copy small reproducibility outputs into the git-backed project, not raw reads."""
import argparse,json,shutil
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('root',type=Path);p.add_argument('out',type=Path);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
def copy(src,dst):
    if not src.is_file():return
    target=a.out/dst;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(src,target)
for name in ['input_manifest.json','alias_map.json']:copy(a.root/'inputs'/name,Path(name))
copy(a.root/'provenance.json',Path('provenance.json'))
for run in (a.root/'runs').iterdir():
    if not run.is_dir():continue
    for name in ['run.json','run.log','started.json','interruption.json']:copy(run/name,Path('runs')/run.name/name)
for run in ['debreak','debreak_no_depth_cap']:
    for name in ['debreak.vcf','log.txt','insertion-merged','deletion-merged']:
        copy(a.root/'runs'/run/'output'/name,Path('native')/run/name)
for run in ['tldr_bounded','tldr_no_detail','tldr_focal']:copy(a.root/'runs'/run/'calls.table.txt',Path('native')/run/'calls.table.txt')
for name in ['summary.json','structural_candidates.tsv','structural_candidates.fasta','read_check.json']:
    copy(a.root/'runs/palidis_ir_mapping'/name,Path('native/palidis_ir_module')/name)
for path in (a.root/'runs/dive_focal_compat/output').glob('*'):
    if path.is_file() and path.stat().st_size<1000000:copy(path,Path('native/dive_focal')/path.name)
for path in (a.root/'comparison').iterdir():
    if path.is_file():copy(path,Path('comparison')/path.name)
print('Collected',len(list(a.out.rglob('*'))),'paths')
