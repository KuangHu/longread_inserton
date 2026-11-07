"""Download check for the read-site search scale-up (2026-11-02; run from the repo root, no network).
Every run of the pool runs.tsv not in the ready manifest, with why it can or cannot feed the main search: the search
needs reads of usable ONT/HiFi runs AND the same sample's assembly in the reviewed manifest (catalog, genes), because
flanks are placed on the sample's own assembly (other assemblies are auxiliary only). Writes site_search_download_check.tsv."""
import csv, gzip, sys
from collections import Counter
sys.path.insert(0, 'experiments/core_context')
import ucd_readcheck as u
D = u.Path('experiments/core_context/scale_comamonadaceae_v0.1')
ready = u.reads_ready.load(u.READS / 'ready_manifest.tsv')
man = {r['sample'] for r in u.pc.read_tsv(D / 'manifest.tsv')}
asm, excl = set(), {}
for a in csv.DictReader(gzip.open(str(u.POOL / 'assemblies.tsv.gz'), 'rt'), delimiter='\t'):
    asm.add(a['biosample'])
for e in csv.DictReader(gzip.open(str(u.POOL / 'excluded.tsv.gz'), 'rt'), delimiter='\t'):
    excl[e['biosample']] = e['reason']
out = []
for r in u.pc.read_tsv(u.POOL / 'runs.tsv'):
    if r['run_accession'] in ready:
        continue
    s = r['sample_accession']
    why = ('run_not_usable:' + r['unusable_reasons'] if r['usable'] != 'True' else
           'no_assembly_in_reviewed_manifest:' + ('assembly_excluded:' + excl[s] if s in excl else
                                                  'assembly_in_taxon_table' if s in asm else 'no_assembly_for_biosample')
           if s not in man else 'read_type_not_ont_or_hifi' if r['read_type'] not in ('ont', 'pacbio_hifi_ccs') else
           'downloadable_for_search')
    out.append(dict(sample=s, run=r['run_accession'], read_type=r['read_type'], usable=r['usable'],
                    fastq_bytes=r['fastq_bytes'], in_manifest=s in man, decision=why))
with open(str(D / 'download_batches' / 'site_search_download_check.tsv'), 'w') as h:
    w = csv.DictWriter(h, list(out[0]), delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(out)
c = Counter(o['decision'].split(':')[0] + (':' + o['decision'].split(':')[-1] if o['decision'].startswith('no_assembly') else '') for o in out)
for k, v in c.most_common():
    gb = sum(int(o['fastq_bytes'] or 0) for o in out if o['decision'].startswith(k.split(':')[0]) and (':' not in k or o['decision'].endswith(k.split(':')[-1]))) / 1e9
    print(v, k, round(gb, 1), 'GB', len({o['sample'] for o in out if o['decision'].startswith(k.split(':')[0]) and (':' not in k or o['decision'].endswith(k.split(':')[-1]))}), 'samples')
