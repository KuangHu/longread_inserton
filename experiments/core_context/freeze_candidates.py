"""Freeze a small batch of representative U-core-D candidates from the committed S6 v0.1 review of the three
subsets, and list the long-read data that already exists for their assemblies. Nothing is downloaded or mapped.
Evidence stays assembly_only; production_ready=false; a candidate is not an independent insertion event and its
ends are candidate conserved ends at most.

Candidate = one single-gene core, two copies in ONE assembly (different coordinates) that are both supported
positions (site group or single-observation site), whose S6 host context (N=5 per side) is different on both
sides, and whose host windows share no physical gene (copies close together can have interlocking windows, e.g. one
copy's downstream genes are the other's upstream genes; such a pair does not show two contexts). Admission uses only
these gene-order / coordinate rules. Unknown boundaries are allowed; empty sites, gaps, copy number,
product and the flank DNA re-check are not admission criteria.

Main batch: the top PER_SUBSET cores per subset by the ranking below (one pair per core). Read-available supplement:
the same rule and ranking restricted to assemblies with a usable linked long-read run, top PER_READ_ASSEMBLY cores per
such assembly not already in the main batch (ids UCDR..). Read availability only chooses this supplement; it is not
evidence.
Ranking (one pair per core):
  1. copies with both boundary sides resolved (2, 1, 0)
  2. flank DNA v2 outer windows not_shared on both sides (auxiliary, ranking only)
  3. more copies of the core in that assembly
  4. identifiers
Long-read listing per assembly: linked runs in the frozen pool metadata (runs.tsv, assembly_run_links), usability,
base count and depth estimate (base count / assembly length), whether the FASTQ is already on local disk with a
verified MD5, and whether any existing alignment covers the candidate (the earlier gap-line genome PAFs hold only
reads preselected at gap alleles, so they are not used as coverage).
"""
import argparse
import json
from collections import Counter, defaultdict
from pathlib import Path

import protein_context as pc

VERSION = 'ucd_candidates_v0.1'
PER_SUBSET = 4
PER_READ_ASSEMBLY = 2
SUPPORTED = {'site_group', 'single_observation_site'}
POOL = Path(__file__).resolve().parent.parent / 'fna_pool' / 'scaleup_stage1_20261026' / 'comamonadaceae'
READS = Path('/global/scratch/users/kh36969/longread_callout/data/reads')


def rank(c):
    return (-c['copies_both_ends_resolved'], not c['dna_v2_not_shared_both'], -c['copies_in_assembly'],
            c['subset'], c['core_key'], c['occurrence_a'], c['occurrence_b'])


def eligible(pair, occ, shared_genes):
    return (pair['structure'] == 'single_gene' and pair['context_class'] == 'context_different_both_sides'
            and occ[pair['occurrence_a']]['s6_status'] in SUPPORTED and occ[pair['occurrence_b']]['s6_status'] in SUPPORTED
            and shared_genes == 0)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--subsets', type=Path, required=True, help='subsets_20261027 directory')
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(exist_ok=True)
    runs = defaultdict(list)
    for r in pc.read_tsv(POOL / 'runs.tsv'):
        runs[r['sample_accession']].append(r)

    cands, supp, pool_counts = [], [], {}
    manifests = {}
    for sub in ('comamonas', 'delftia', 'paracidovorax'):
        s6 = a.subsets / sub / 'protein_entry' / 's6'
        occ = {o['occurrence_id']: o for o in pc.read_tsv(s6 / 's6_occurrences.tsv.gz')}
        cores = {c['core_key']: c for c in pc.read_tsv(s6 / 's6_cores.tsv') if c['structure'] == 'single_gene'}
        dna = {(r['occurrence_a'], r['occurrence_b']): r for r in pc.read_tsv(s6 / 'review' / 'flank_dna_v2.tsv.gz')
               if r['source'] == 'same_assembly'}
        man = {m['sample']: m for m in pc.read_tsv(a.subsets / (sub + '.samples.tsv'))}
        manifests[sub] = man
        per_asm = Counter((o['core_key'], o['sample']) for o in occ.values() if o['structure'] == 'single_gene')
        best, best_asm = {}, {}
        n = 0
        for r in pc.read_tsv(s6 / 'review' / 's6_same_assembly_pairs.tsv.gz'):
            d = dna[(r['occurrence_a'], r['occurrence_b'])]
            if not eligible(r, occ, int(d['host_genes_shared_physically'])):
                continue
            n += 1
            x, y = occ[r['occurrence_a']], occ[r['occurrence_b']]
            res = lambda o: o['boundary_up'] == o['boundary_down'] == 'resolved'
            c = cores[r['core_key']]
            row = dict(subset=sub, core_key=r['core_key'], products=c['products'], product_hint=c['product_hint'],
                copies_per_assembly=c['copies_per_assembly'], sample=r['sample'], assembly=man[r['sample']]['assembly'],
                organism_label=man[r['sample']]['organism'], copies_in_assembly=per_asm[(r['core_key'], r['sample'])],
                occurrence_a=r['occurrence_a'], occurrence_b=r['occurrence_b'])
            for k, o in (('a', x), ('b', y)):
                row.update({'loc_' + k: '%s:%d-%s(%s)' % (o['contig'], int(o['start0']) + 1, o['end0'], o['strand']),
                            's6_' + k: o['site_group'] or 'single_observation_site',
                            'coding_span_' + k: o['coding_span'], 'boundary_' + k: '%s/%s' % (o['boundary_up'], o['boundary_down']),
                            'shared_estimate_' + k: o['estimated_span'],
                            'host_up_' + k: o['host_up'], 'host_down_' + k: o['host_down']})
            row.update(copy_gap_bp=d['copy_gap_bp'] or 'different_contigs', copies_both_ends_resolved=res(x) + res(y),
                       dna_v2='%s/%s' % (d['dna_up'], d['dna_down']), dna_v2_not_shared_both=d['dna_up'] == d['dna_down'] == 'not_shared',
                       dna_v2_other_orientation_match=d['other_orientation_match'] == 'True')
            k = (sub, r['core_key'])
            if k not in best or rank(row) < rank(best[k]):
                best[k] = row
            k = (sub, r['core_key'], r['sample'])
            if any(x['usable'] == 'True' for x in runs.get(r['sample'], [])) and (k not in best_asm or rank(row) < rank(best_asm[k])):
                best_asm[k] = row
        pool_counts[sub] = dict(eligible_pairs=n, eligible_cores=len(best))
        main_ = sorted(best.values(), key=rank)[:PER_SUBSET]
        for c in main_:
            c['batch'] = 'main'
        cands += main_
        for smp in sorted({k[2] for k in best_asm}):
            xs = [v for k, v in best_asm.items() if k[2] == smp and all(c['core_key'] != v['core_key'] for c in main_)]
            sup = sorted(xs, key=rank)[:PER_READ_ASSEMBLY]
            for c in sup:
                c['batch'] = 'read_available_supplement'
            supp.extend(sup)
            pool_counts[sub]['eligible_cores_in_read_assembly_' + smp] = len({k[1] for k in best_asm if k[2] == smp})

    n = 0
    for c in cands:
        n += 1
        c['candidate'] = 'UCD%02d' % n
    for i, c in enumerate(supp):
        c['candidate'] = 'UCDR%02d' % (i + 1)
    cands += supp
    for c in cands:
        c['label'] = 'assembly_only; candidate; not an independent insertion event'
    cols = ['candidate', 'batch'] + [k for k in cands[0] if k not in ('candidate', 'batch')]
    pc.write_tsv(a.out / 'candidates.tsv', cands, cols)

    reads = []
    for sub, man in manifests.items():
        for smp, m in sorted(man.items()):
            rs = runs.get(smp, [])
            base = {r['run_accession']: r for r in rs}
            local = {}
            for run in base:
                v = READS / run / 'verify.json'
                if v.exists():
                    j = json.loads(v.read_text())
                    local[run] = 'local_%s_md5_%s' % (j['status'], j['md5'] == base[run]['fastq_md5'])
                else:
                    local[run] = 'not_local'
            reads.append(dict(subset=sub, sample=smp, assembly=m['assembly'], total_bp=m['total_bp'],
                candidates=';'.join(c['candidate'] for c in cands if c['sample'] == smp),
                linked_runs=len(rs), usable_runs=sum(r['usable'] == 'True' for r in rs),
                runs=';'.join('%s:%s:%s:%s' % (r['run_accession'], r['read_type'], r['instrument_model'],
                                               'usable' if r['usable'] == 'True' else 'unusable(' + r['unusable_reasons'] + ')') for r in rs),
                usable_bases=sum(int(r['base_count'] or 0) for r in rs if r['usable'] == 'True'),
                depth_estimate=round(sum(int(r['base_count'] or 0) for r in rs if r['usable'] == 'True') / int(m['total_bp']), 1),
                fastq_local=';'.join('%s:%s' % kv for kv in local.items()),
                alignment_at_candidates='none (gap-line genome PAFs hold only reads preselected at gap alleles)' if rs else ''))
    pc.write_tsv(a.out / 'long_read_availability.tsv', reads, list(reads[0]))
    summary = dict(version=VERSION, evidence_source='assembly_only', production_ready=False, per_subset=PER_SUBSET,
                   rule=__doc__, pool=pool_counts, candidates=len(cands) - len(supp), read_available_supplement=len(supp),
                   candidates_with_usable_long_reads=sum(any(r['usable_runs'] for r in reads if r['sample'] == c['sample']) for c in cands),
                   downloads_started=False)
    (a.out / 'freeze.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
