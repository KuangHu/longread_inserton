"""Per-assembly inventory of the wave-1 FNA pools: which assemblies only have the old FNA / gap pipeline (s1-s4), and
which are in the U-core-D catalogue and have completed the two current checks (ucd_batch).

  python pool_inventory.py --out scale_comamonadaceae_v0.1/pool_inventory   (reads the final ucd_batch rows)

Status per assembly (pool row = BioSample + chosen accession.version):
  ucd_checks_done        in the U-core-D catalogue manifest and a side of >= 1 pair whose checks ran
  ucd_catalogue_no_pair  in the catalogue manifest, but in no catalogue pair (no core shared with another occurrence)
  ucd_pool_not_included  in a family whose catalogue exists, but left out of its manifest (reason given)
  old_pipeline_only      s1-s4 (download, joint FNA / skani, gap caller, locus annotation) only; no U-core-D catalogue
"""
import argparse
import gzip
import json
from collections import Counter, defaultdict
from pathlib import Path

import ucd_compare as uc

POOLS = uc.REPO.parent.parent / 'fna_pool' / 'scaleup_stage1_20261026'
FINAL = uc.SCRATCH / 'ucd_batch_v0.2' / 'final'
CATALOGUES = {'comamonadaceae': uc.REPO / 'manifest.tsv'}


def side_counts():
    c = defaultdict(Counter)
    for kind in uc.KINDS:
        with gzip.open(FINAL / ('all_rows_%s.tsv.gz' % kind), 'rt') as f:
            h = f.readline().rstrip('\n').split('\t')
            ia, ib, isc, ifs = h.index('a_sample'), h.index('b_sample'), h.index('source_check'), h.index('final_status')
            for line in f:
                v = line.split('\t')
                for smp in {v[ia], v[ib]}:
                    c[smp]['%s_rows' % kind] += 1
                    if v[isc] in ('pass', 'fail'):
                        c[smp]['%s_checked' % kind] += 1
                    if v[ifs] == 'pass':
                        c[smp]['%s_pass' % kind] += 1
    with gzip.open(FINAL / 'strict_lt2500_candidate_insertion.tsv.gz', 'rt') as f:
        h = f.readline().rstrip('\n').split('\t')
        ia, ib = h.index('a_sample'), h.index('b_sample')
        for line in f:
            v = line.split('\t')
            for smp in {v[ia], v[ib]}:
                c[smp]['strict_candidate_pass'] += 1
    return c


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--out', required=True)
    a = p.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    sc = side_counts()
    rows, summ = [], {}
    for fam in sorted(d.name for d in POOLS.iterdir() if (d / 'pool.tsv').exists()):
        d = POOLS / fam
        acc = {r['accession']: r for r in uc.read_tsv(d / 's1' / 'accessions.tsv.gz')}
        s2 = {r['sample_id'] for r in uc.read_tsv(d / 's2' / 'catalog.tsv')}
        s3 = Counter()
        for r in uc.read_tsv(d / 's3_full' / 'pairs.tsv.gz'):
            s3[r['sample_a']] += 1
            s3[r['sample_b']] += 1
        s4 = Counter(r['carrier'] for r in uc.read_tsv(d / 's4_s3_full' / 's4_loci.tsv.gz'))
        man = {r['sample']: r for r in uc.read_tsv(CATALOGUES[fam])} if fam in CATALOGUES else None
        for r in uc.read_tsv(d / 'pool.tsv'):
            bs, ac = r['biosample'], r['download_accession']
            s1 = acc.get(ac, {})
            x = dict(family=fam, biosample=bs, accession=ac, genbank_accession=r['genbank_accession'],
                     refseq_accession=r['refseq_accession'], assembly_level=r['assembly_level'],
                     atypical=r['atypical'], n_usable_long_read_runs=r['n_usable_runs'],
                     s1_status=s1.get('status', 'missing'), annotation_status=s1.get('annotation_status', ''),
                     s2_catalog=bs in s2, s3_gap_pairs=s3[bs], s4_loci=s4[bs])
            c = sc.get(bs, Counter())
            if man is None:
                st = 'old_pipeline_only'
            elif bs not in man:
                st = 'ucd_pool_not_included'
                x['ucd_note'] = ('no_gff3_or_faa' if s1.get('annotation_status') != 'present' else
                                 'atypical' if r['atypical'] else 'not_in_manifest')
            else:
                m = man[bs]
                x.update(ucd_manifest_assembly=m['assembly'], ucd_batch=m['batch'],
                         ucd_assembly_matches_pool=m['assembly'] == ac)
                st = 'ucd_checks_done' if c['cds_core_checked'] else 'ucd_catalogue_no_pair'
            x.update(status=st, **{k: c[k] for k in ('cds_core_rows', 'cds_core_checked', 'cds_core_pass',
                                                     'candidate_insertion_checked', 'candidate_insertion_pass',
                                                     'strict_candidate_pass')})
            rows.append(x)
        fr = [x for x in rows if x['family'] == fam]
        summ[fam] = dict(pool=len(fr), status=dict(Counter(x['status'] for x in fr)),
                         s1_verified=sum(x['s1_status'] == 'verified' for x in fr),
                         s2=sum(x['s2_catalog'] for x in fr), s3_assemblies_in_gap_pairs=sum(x['s3_gap_pairs'] > 0
                                                                                            for x in fr),
                         s4_assemblies_with_loci=sum(x['s4_loci'] > 0 for x in fr),
                         notes=dict(Counter(x.get('ucd_note', '') for x in fr if x.get('ucd_note'))))
    cols = ['family', 'biosample', 'accession', 'genbank_accession', 'refseq_accession', 'assembly_level', 'atypical',
            'n_usable_long_read_runs', 's1_status', 'annotation_status', 's2_catalog', 's3_gap_pairs', 's4_loci',
            'status', 'ucd_note', 'ucd_manifest_assembly', 'ucd_batch', 'ucd_assembly_matches_pool', 'cds_core_rows',
            'cds_core_checked', 'cds_core_pass', 'candidate_insertion_checked', 'candidate_insertion_pass',
            'strict_candidate_pass']
    uc.write_tsv(out / 'assemblies.tsv.gz', rows, cols)
    man = uc.read_tsv(CATALOGUES['comamonadaceae'])
    summ['_catalogue_manifest'] = dict(path=str(CATALOGUES['comamonadaceae']), file_lines=len(man) + 1,
                                      assemblies=len(man), distinct_samples=len({m['sample'] for m in man}),
                                      distinct_accessions=len({m['assembly'] for m in man}),
                                      note='1,683 is the line count of manifest.tsv including its header; the '
                                           'catalogue has 1,682 assemblies')
    (out / 'summary.json').write_text(json.dumps(summ, indent=1) + '\n')
    print(json.dumps(summ, indent=1))


if __name__ == '__main__':
    main()
