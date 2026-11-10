"""Recover the 2026-10-24 local H16/H850 smoke without gap calls or reads.

Only uses already downloaded pilot data and tools. Windows uses its existing WSL
Ubuntu runtime. This script does not download anything or change the old pilot.
python reproduce_local_smoke.py --pilot /path/to/cross_assembly_pilot --out /new/path
"""
import argparse
import json
import os
import subprocess
import time
from pathlib import Path

import core_context as cc


def linux(path):
    p = Path(path).resolve().as_posix()
    return '/mnt/' + p[0].lower() + p[2:]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--pilot', type=Path, required=True); p.add_argument('--out', type=Path, required=True)
    a = p.parse_args(); a.pilot = a.pilot.resolve(); a.out = a.out.resolve()
    a.out.mkdir(parents=True, exist_ok=False)
    bundle = a.out/'bundle'; bundle.mkdir()
    fasta = a.pilot/'data/assemblies.fna'
    genomes = cc.read_fasta(fasta)
    mapping = []
    for name, sequence in genomes.items():
        sample, contig = name.split('__', 1)
        mapping.append(dict(graph_name=name, sample=sample, contig=contig, length=len(sequence), sha256=cc.sha(sequence)))
    (bundle/'assemblies.fna').write_bytes(fasta.read_bytes())
    cc.write_tsv(bundle/'contigs.tsv', mapping, ['graph_name','sample','contig','length','sha256'])
    binary = a.pilot/'bin/pangraph'
    if os.name == 'nt':
        prefix = ['wsl', '-d', 'Ubuntu', '--', linux(binary)]
        fna, graph_path = linux(bundle/'assemblies.fna'), linux(bundle/'graph.json')
    else:
        prefix = [str(binary)]; fna, graph_path = str(bundle/'assemblies.fna'), str(bundle/'graph.json')
    version = subprocess.run(prefix+['--version'], capture_output=True, text=True, check=True).stdout.strip()
    command = prefix+['build','-j','4','--verify','--no-progress-bar',fna,'-o',graph_path]
    start = time.time()
    with open(bundle/'build.stderr.txt','w') as err:
        subprocess.run(command, stderr=err, check=True)
    cc.write_json(bundle/'graph_provenance.json', dict(tool_version=version, command=command,
        inputs=[dict(path=str(fasta), sha256=cc.digest(fasta))],
        graph_sha256=cc.digest(bundle/'graph.json'), joint_fasta_sha256=cc.digest(bundle/'assemblies.fna'),
        contigs_sha256=cc.digest(bundle/'contigs.tsv'), graph_verify_passed=True,
        binary_sha256=cc.digest(binary), build_seconds=round(time.time()-start,3)))
    summary = cc.call(bundle, a.out/'results')
    # Compare AFTER discovery; this control never enters the caller.
    seeds = json.loads((a.pilot/'results/seeds.json').read_text())
    control = next(x for x in seeds if x['seed_id']=='CONTROL_H16_R33')
    occurrences = json.loads((a.out/'results/occurrences.json').read_text())
    source = [x for x in occurrences if x['sample']==control['filled_sample'] and x['contig']==control['filled_contig']
              and max(x['start0'],control['start0']) < min(x['end0'],control['end0'])]
    matches = []
    for x in source:
        overlap = min(x['end0'],control['end0'])-max(x['start0'],control['start0'])
        if overlap/min(x['core_length'],control['length']) >= .9:
            matches.append(x)
    recovered = {x['family_id'] for x in matches}
    xs = [x for x in occurrences if x['family_id'] in recovered]
    old = json.loads((a.pilot/'results/reviewed_contexts.json').read_text())
    old = next(x for x in old if x['seed_id']=='CONTROL_H16_R33')['loci']
    compared = []
    for r in old:
        hits = [x for x in xs if x['sample']==r['sample'] and x['contig']==r['contig']
                and abs(x['start0']-r['start0']) <= 20 and abs(x['end0']-r['end0']) <= 20]
        compared.append(dict(sample=r['sample'],contig=r['contig'],start0=r['start0'],end0=r['end0'],
            recovered_occurrences=[x['occurrence_id'] for x in hits]))
    report = dict(discovery_inputs=['joint assembly FASTA only'], gap_calls_used=False, reads_used=False,
        seed_control_used_for_discovery=False, control_used_for_posthoc_regression=control['seed_id'],
        control_families=sorted(recovered), control_occurrences=len(xs),
        old_control_loci=len(old), old_control_loci_recovered=sum(bool(x['recovered_occurrences']) for x in compared),
        comparison=compared, summary=summary)
    cc.write_json(a.out/'results/smoke_regression.json', report)
    cc.write_tsv(a.out/'results/control_occurrences.tsv',xs,['family_id','occurrence_id','sample','contig','start0','end0','strand','core_length'])
    # Plain SVG: genomic intervals and core direction, no implied read validation.
    from html import escape
    svg = ['<svg xmlns="http://www.w3.org/2000/svg" width="1120" height="%d">' % (110+50*len(xs)),
        '<rect width="100%" height="100%" fill="white"/>',
        '<text x="20" y="28" font-size="18">R33 control recovered by gap-independent graph discovery</text>',
        '<text x="20" y="52" font-size="12">Assembly only; graph boundaries provisional; boxes not to scale; 1-2 CDS not yet checked.</text>']
    for i,x in enumerate(xs):
        y=85+i*50
        label='%s %s:%d-%d (%s)'%(x['sample'],x['contig'],x['start0']+1,x['end0'],x['strand'])
        svg.append('<text x="20" y="%d" font-size="12">%s</text>'%(y+18,escape(label)))
        for left,width,color,txt in [(450,170,'#b6d1df','upstream 500 bp'),(635,210,'#edc06f','core %d bp'%x['core_length']),
                                      (860,170,'#b4d4b9','downstream 500 bp')]:
            svg.append('<rect x="%d" y="%d" width="%d" height="30" fill="%s"/>'%(left,y,width,color))
            svg.append('<text x="%d" y="%d" font-size="12">%s</text>'%(left+8,y+19,txt))
    svg.append('</svg>'); (a.out/'results/control.svg').write_text('\n'.join(svg))
    print(json.dumps({k:v for k,v in report.items() if k!='comparison'},indent=2))
    if len(old)!=6 or report['old_control_loci_recovered']!=6:
        raise RuntimeError('Old six-location smoke not recovered; inspect smoke_regression.json')


if __name__=='__main__': main()
