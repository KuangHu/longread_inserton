"""Direct core comparison of the two copies of every frozen U-core-D candidate (ucd_candidates_v0.1), used as the
regression set for scale-up. Evidence stays assembly_only; production_ready=false.

The core review measured boundaries only on pairs that protein entry called different_context_candidate (K=3), so a
frozen pair may have no direct row there (UCDR03). Here every candidate pair is compared directly with the core
review's own functions (core_review.blast / chain / extensions, unchanged):
  core protein   identical sequence (sha256) or blastp identity and min coverage
  core DNA       identical, or blastn identity and min coverage of the oriented coding spans
  pair boundary  the two copies' coding span +/- W bp aligned; per side the shared extension and its cover status,
                 measured only if that end is covered in both copies (unknown otherwise, never 0); pair-level only
  context        K=3 near relation (protein entry rule), S6 host side relations and supports, flank DNA v2 and the
                 physical-overlap flags, the per-copy boundary status of the core review (all partners)
  in_review      whether core_review has a direct row for this pair
"""
import argparse
import json
from argparse import Namespace
from collections import defaultdict
from pathlib import Path

import core_review as cr
import protein_context as pc


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--subsets', type=Path, required=True)
    p.add_argument('--work-root', type=Path, required=True, help='<root>/<subset>/protein work dirs')
    p.add_argument('--mmseqs', default='/global/scratch/users/kh36969/longread_callout/env/s4/bin/mmseqs')
    a = p.parse_args()
    out_dir = a.subsets / 'ucd_candidates_v0.1'
    cands = pc.read_tsv(out_dir / 'candidates.tsv')
    rows = []
    for sub in sorted({c['subset'] for c in cands}):
        work = a.work_root / sub / 'protein'
        pe = a.subsets / sub / 'protein_entry'
        genes = {g['gene_uid']: g for g in pc.read_tsv(work / 'genes.tsv.gz')}
        assign, _, _ = pc.assignments(Namespace(work=work, mmseqs=a.mmseqs, threads=1))
        prot, name = {}, None
        for line in open(work / 'proteins.faa'):
            if line.startswith('>'):
                name = line[1:].strip()
            else:
                prot[name] = line.strip()
        occs = {o['occurrence_id']: o for o in pc.read_tsv(pe / 'protein_occurrences.tsv.gz')}
        near = {frozenset((r['occurrence_a'], r['occurrence_b'])): r['relation'] for r in pc.read_tsv(pe / 'protein_pairs.tsv.gz')
                if r['scope'] == 'target_lt2500'}
        s6 = {frozenset((r['occurrence_a'], r['occurrence_b'])): r for r in pc.read_tsv(pe / 's6' / 's6_pairs.tsv.gz')}
        dna = {frozenset((r['occurrence_a'], r['occurrence_b'])): r for r in pc.read_tsv(pe / 's6' / 'review' / 'flank_dna_v2.tsv.gz')
               if r['source'] == 'same_assembly'}
        rv_pairs = {frozenset((r['occurrence_a'], r['occurrence_b'])) for r in pc.read_tsv(pe / 'core_review' / 'review_pairs.tsv.gz')}
        rv_occ = {r['occurrence_id']: r for r in pc.read_tsv(pe / 'core_review' / 'review_occurrences.tsv.gz')}
        fasta = {m['sample']: m['fasta'] for m in pc.read_tsv(a.subsets / (sub + '.samples.tsv'))}
        seqs = {}
        for c in [c for c in cands if c['subset'] == sub]:
            A, B = occs[c['occurrence_a']], occs[c['occurrence_b']]
            for o in (A, B):
                if o['sample'] not in seqs:
                    seqs[o['sample']] = pc.read_fna(fasta[o['sample']])
                g = genes[o['genes'].split(';')[0]]
                o['ext'], o['c0'], o['c1'] = cr.window(seqs[o['sample']][o['contig']], int(o['start0']), int(o['end0']),
                                                       o['strand'] == '+', cr.W, g['circular_contig'] == 'True')
                o['core_dna'] = o['ext'][o['c0']:o['c1']]
                o['prot'] = prot[genes[o['genes']]['protein_sha256']]
            pb = cr.best(cr.blast('blastp', {'a': A['prot'], 'b': B['prot']}, None)[('a', 'b')])
            db = cr.best(cr.blast('blastn', {'a': A['core_dna'], 'b': B['core_dna']}, 'blastn')[('a', 'b')])
            hs = cr.blast('blastn', {'a': A['ext'], 'b': B['ext']}, 'blastn')[('a', 'b')]
            blk = cr.chain([(h['qs'], h['qe'], h['ss'], h['se']) for h in hs], A['c0'], A['c1'], B['c0'], B['c1'])
            ext = cr.extensions(blk, A['c0'], A['c1'], len(A['ext']), B['c0'], B['c1'], len(B['ext'])) if blk else \
                dict(up_cover='no_core_alignment', down_cover='no_core_alignment', a_up=None, b_up=None, a_down=None, b_down=None)
            k = frozenset((A['occurrence_id'], B['occurrence_id']))
            s, d = s6[k], dna.get(k, {})
            same_fam = assign['family'][A_sha] if (A_sha := genes[A['genes']]['protein_sha256']) else ''
            rows.append(dict(candidate=c['candidate'], batch=c['batch'], subset=sub, core_key=c['core_key'], sample=A['sample'],
                occurrence_a=A['occurrence_id'], occurrence_b=B['occurrence_id'],
                same_family=same_fam == assign['family'].get(genes[B['genes']]['protein_sha256'], '-'),
                identical_protein=A['prot'] == B['prot'],
                protein_identity=pb[0] if pb else '', protein_coverage_min=round(min(pb[1], pb[2]), 4) if pb else '',
                identical_core_dna=A['core_dna'].upper() == B['core_dna'].upper(),
                core_dna_identity=db[0] if db else '', core_dna_coverage_min=round(min(db[1], db[2]), 4) if db else '',
                coding_span_a=len(A['core_dna']), coding_span_b=len(B['core_dna']),
                pair_up_cover=ext['up_cover'], pair_ext_up_a=ext['a_up'], pair_ext_up_b=ext['b_up'],
                pair_down_cover=ext['down_cover'], pair_ext_down_a=ext['a_down'], pair_ext_down_b=ext['b_down'],
                near_relation_k3=near.get(k, 'not_in_pair_table'),
                s6_verdict=s['verdict'], s6_up=s['up_side'], s6_up_support=s['up_support'], s6_down=s['down_side'],
                s6_down_support=s['down_support'],
                copy_gap_bp=d.get('copy_gap_bp', ''), host_genes_shared_physically=d.get('host_genes_shared_physically', ''),
                windows_can_overlap=d.get('windows_can_overlap', ''), dna_v2=('%s/%s' % (d['dna_up'], d['dna_down'])) if d else '',
                dna_v2_other_orientation_max=max(float(d['same_side_reverse_max']), float(d['cross_side_max'])) if d else '',
                boundary_copy_a='%s/%s' % (rv_occ[A['occurrence_id']]['up_status'], rv_occ[A['occurrence_id']]['down_status'])
                if A['occurrence_id'] in rv_occ else 'not_reviewed_no_different_partner',
                boundary_copy_b='%s/%s' % (rv_occ[B['occurrence_id']]['up_status'], rv_occ[B['occurrence_id']]['down_status'])
                if B['occurrence_id'] in rv_occ else 'not_reviewed_no_different_partner',
                pair_in_core_review=k in rv_pairs))
    pc.write_tsv(out_dir / 'regression_core_check.tsv', rows, list(rows[0]))
    main_rows = [r for r in rows if r['batch'] == 'main']
    summary = dict(evidence_source='assembly_only', production_ready=False, pairs=len(rows),
                   main=dict(n=len(main_rows), identical_core_dna=sum(r['identical_core_dna'] for r in main_rows),
                             core_dna_identity_min=min(r['core_dna_identity'] for r in main_rows),
                             core_dna_coverage_min=min(r['core_dna_coverage_min'] for r in main_rows),
                             protein_identity_min=min(r['protein_identity'] for r in main_rows),
                             near_relation=dict(defaultdict(int, {}))),
                   not_in_core_review=[r['candidate'] for r in rows if not r['pair_in_core_review']])
    for r in main_rows:
        summary['main']['near_relation'][r['near_relation_k3']] = summary['main']['near_relation'].get(r['near_relation_k3'], 0) + 1
    (out_dir / 'regression_core_check.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
