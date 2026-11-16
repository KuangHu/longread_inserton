"""Evidence tiers and a core list from the scale review (after scale_review boundary), plus the first long-read
validation batch from reads that are already local. Evidence stays assembly_only; production_ready=false; counts are
pairs, occurrences or cores, never sites or insertion events. Nothing is downloaded or mapped.

Per reviewed pair (same-assembly U-core-D structural pairs and cross-assembly pairs, both S6 grouping tiers), three
layers kept apart:
  context    different on both sides (every reviewed pair; scale_review s6)
  core match core_match_class: identical_core_dna / homolog_full_coverage = sufficient; partial_core_match and
             no_core_dna_hit are not (partial matches stay pending review)
  length     both occurrences' estimated span, from the boundary of that pair's scope: lt2500 for both, else ge2500 for
             one, else unknown (pending review). A short coding core never stands in for a resolved span.
evidence_state: partial_core_match_pending | no_core_dna_hit | match_boundary_unknown_pending | match_span_ge2500 |
                match_both_lt2500.

Core list: every core with at least one match_both_lt2500 pair, with its pairs per scope / copy class / S6 grouping,
distinct occurrences and assemblies of those pairs, S6 status of those occurrences, the pairs of the core still
pending, and the assemblies of those occurrences with verified local long reads (data/reads/<run>/verify.json,
status verified and MD5 equal to the archive MD5 in the pool runs.tsv).

Validation batch 1: every occurrence of a match_both_lt2500 pair whose assembly has verified local reads, with its
coordinates, estimated span, partners in that tier and the local runs. Long-read checks themselves are a later step.
"""
import argparse
import json
from collections import Counter, defaultdict
from pathlib import Path

import protein_context as pc
import scale_entry as se

VERSION = 'scale_candidates_v0.1'
POOL = Path(__file__).resolve().parent.parent / 'fna_pool' / 'scaleup_stage1_20261026' / 'comamonadaceae'
READS = Path('/global/scratch/users/kh36969/longread_callout/data/reads')
SUFFICIENT = {'identical_core_dna', 'homolog_full_coverage'}
STATES = ('match_both_lt2500', 'match_span_ge2500', 'match_boundary_unknown_pending', 'partial_core_match_pending',
          'no_core_dna_hit')


def evidence_state(match, span_a, span_b):
    if match == 'partial_core_match':
        return 'partial_core_match_pending'
    if match not in SUFFICIENT:
        return 'no_core_dna_hit'
    if 'unknown' in (span_a, span_b):
        return 'match_boundary_unknown_pending'
    if 'ge2500' in (span_a, span_b):
        return 'match_span_ge2500'
    return 'match_both_lt2500'


def local_reads(runs_tsv, reads_dir):
    """{sample: [run descriptors]} for runs with a verified local FASTQ whose MD5 equals the archive MD5."""
    out = defaultdict(list)
    for r in pc.read_tsv(runs_tsv):
        v = reads_dir / r['run_accession'] / 'verify.json'
        if not v.exists():
            continue
        j = json.loads(v.read_text())
        if j.get('status') == 'verified' and j.get('md5') == r['fastq_md5']:
            out[r['sample_accession']].append('%s:%s:%s:%s' % (r['run_accession'], r['read_type'], r['instrument_model'],
                                                               'usable' if r['usable'] == 'True' else 'unusable'))
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--work', type=Path, required=True)
    p.add_argument('--batches', required=True)
    p.add_argument('--repo-out', type=Path, required=True)
    p.add_argument('--runs', type=Path, default=POOL / 'runs.tsv')
    p.add_argument('--reads', type=Path, default=READS)
    a = p.parse_args()
    tag = 'b' + '-'.join(a.batches.split(','))
    reads = local_reads(a.runs, a.reads)
    pair_states = Counter()
    pair_states_by = Counter()
    core_rows, val_rows, good_pairs = [], [], []
    occ_states = defaultdict(Counter)
    missing = []
    for k in range(se.SHARDS):
        bd = a.work / 'review' / 'boundary' / tag / ('shard_%03d' % k)
        sd = a.work / 'review' / 's6' / tag / ('shard_%03d' % k)
        if not (bd / 'run.json').exists():
            missing.append(k); continue
        span = {(r['occurrence_id'], r['pair_scope']): r for r in pc.read_tsv(bd / 'occurrences.tsv.gz')}
        s6 = {r['occurrence_id']: r for r in pc.read_tsv(sd / 'occurrences.tsv.gz')}
        cores = {(r['structure'], r['core_key']): r for r in pc.read_tsv(sd / 'cores.tsv.gz')}
        per_core = defaultdict(list)
        for r in pc.read_tsv(bd / 'pairs.tsv.gz'):
            sa = span[(r['occurrence_a'], r['pair_scope'])]['span_class']
            sb = span[(r['occurrence_b'], r['pair_scope'])]['span_class']
            st = evidence_state(r['core_match_class'], sa, sb)
            pair_states[(r['structure'], r['pair_scope'], st)] += 1
            pair_states_by[(r['structure'], r['pair_scope'], r['copy_class'] or 'same_assembly', r['s6_grouping'], st)] += 1
            for x in (r['occurrence_a'], r['occurrence_b']):
                occ_states[(r['structure'], x)][st] += 1
            per_core[(r['structure'], r['core_key'])].append((r, st))
        for key in sorted(per_core):
            rows = per_core[key]
            good = [r for r, st in rows if st == 'match_both_lt2500']
            if not good:
                continue
            c = cores[key]
            occs = sorted({x for r in good for x in (r['occurrence_a'], r['occurrence_b'])})
            smps = sorted({s6[x]['sample'] for x in occs})
            with_reads = [s for s in smps if reads.get(s)]
            best = min(good, key=lambda r: (r['core_match_class'] != 'identical_core_dna', -float(r['core_dna_identity'] or 0),
                                            r['pair_scope'] != 'same_assembly', r['occurrence_a'], r['occurrence_b']))
            n = Counter(st for _, st in rows)
            core_rows.append(dict(
                structure=key[0], core_key=key[1], product=c['products'].split(' | ')[0], product_hint=c['product_hint'],
                k3_different_pairs=c['k3_different_pairs'],
                pairs_match_both_lt2500=len(good),
                same_assembly=sum(r['pair_scope'] == 'same_assembly' for r in good),
                cross_single_x_single=sum(r['copy_class'] == 'single_x_single' for r in good),
                cross_single_x_multi=sum(r['copy_class'] == 'single_x_multi' for r in good),
                cross_multi_x_multi=sum(r['copy_class'] == 'multi_x_multi' for r in good),
                s6_supported_pairs=sum(r['s6_grouping'] == 'supported_positions' for r in good),
                identical_core_dna_pairs=sum(r['core_match_class'] == 'identical_core_dna' for r in good),
                core_dna_identity_min=min(float(r['core_dna_identity']) for r in good if r['core_dna_identity'] != ''),
                occurrences=len(occs), assemblies=len(smps),
                occurrence_s6_status=','.join('%s=%d' % kv for kv in sorted(Counter(s6[x]['s6_status'] for x in occs).items())),
                pairs_reviewed=len(rows), pairs_span_ge2500=n['match_span_ge2500'],
                pairs_boundary_unknown_pending=n['match_boundary_unknown_pending'],
                pairs_partial_core_match_pending=n['partial_core_match_pending'],
                best_pair='%s~%s' % (best['occurrence_a'], best['occurrence_b']), best_pair_scope=best['pair_scope'],
                best_pair_core_match=best['core_match_class'], best_pair_core_dna_identity=best['core_dna_identity'],
                assemblies_with_local_reads=len(with_reads), local_read_assemblies=';'.join(with_reads)))
            partners = defaultdict(list)
            for r in good:
                partners[r['occurrence_a']].append(r); partners[r['occurrence_b']].append(r)
                good_pairs.append(dict({f: r[f] for f in ('structure', 'core_key', 'pair_scope', 'copy_class', 's6_grouping',
                                                         'sample_a', 'sample_b', 'occurrence_a', 'occurrence_b',
                                                         'core_match_class', 'core_dna_identity', 'core_dna_coverage_min',
                                                         'protein_identity_min')},
                                       estimated_span_a=span[(r['occurrence_a'], r['pair_scope'])]['estimated_span'],
                                       estimated_span_b=span[(r['occurrence_b'], r['pair_scope'])]['estimated_span']))
            for x in occs:
                o = s6[x]
                if not reads.get(o['sample']):
                    continue
                ps = partners[x]
                scopes = sorted({r['pair_scope'] for r in ps})
                val_rows.append(dict(
                    structure=key[0], core_key=key[1], product=c['products'].split(' | ')[0], occurrence_id=x,
                    sample=o['sample'], contig=o['contig'], start0=o['start0'], end0=o['end0'], strand=o['strand'],
                    s6_status=o['s6_status'],
                    estimated_span=';'.join('%s:%s' % (sc, span[(x, sc)]['estimated_span']) for sc in scopes),
                    partners_match_both_lt2500=len(ps),
                    same_assembly_partners=sum(r['pair_scope'] == 'same_assembly' for r in ps),
                    cross_assembly_partners=sum(r['pair_scope'] == 'cross_assembly' for r in ps),
                    partner_with_local_reads=sum(bool(reads.get(r['sample_b' if r['occurrence_a'] == x else 'sample_a']))
                                                 for r in ps),
                    local_runs=';'.join(reads[o['sample']])))
    out = dict(version=VERSION, batches=a.batches, evidence_source='assembly_only', production_ready=False,
               counts_are_pairs_occurrences_or_cores_not_sites_or_events=True, boundary_shards_missing=missing,
               sufficient_core_match=sorted(SUFFICIENT), assemblies_with_verified_local_reads=len(reads), structures={})
    for st in ('single_gene', 'two_gene'):
        cr_ = [c for c in core_rows if c['structure'] == st]
        vr = [v for v in val_rows if v['structure'] == st]
        out['structures'][st] = dict(
            pairs_by_state={sc: {s: pair_states[(st, sc, s)] for s in STATES} for sc in ('same_assembly', 'cross_assembly')},
            pairs_by_class_grouping_state={'%s:%s:%s:%s' % k[1:]: n for k, n in sorted(pair_states_by.items()) if k[0] == st},
            distinct_occurrences_in_a_match_both_lt2500_pair=sum(1 for (s, _), c in occ_states.items()
                                                                 if s == st and c['match_both_lt2500']),
            distinct_occurrences_only_pending=sum(1 for (s, _), c in occ_states.items() if s == st and not c['match_both_lt2500']
                                                  and (c['match_boundary_unknown_pending'] or c['partial_core_match_pending'])),
            cores_with_match_both_lt2500_pair=len(cr_),
            cores_by_hint=dict(Counter(c['product_hint'] for c in cr_)),
            cores_with_same_assembly_pair=sum(c['same_assembly'] > 0 for c in cr_),
            cores_only_cross_assembly=sum(c['same_assembly'] == 0 for c in cr_),
            cores_with_s6_supported_pair=sum(c['s6_supported_pairs'] > 0 for c in cr_),
            cores_where_k3_had_no_different_pair=sum(c['k3_different_pairs'] == '0' for c in cr_),
            cores_with_local_reads=sum(c['assemblies_with_local_reads'] > 0 for c in cr_),
            validation_batch1_occurrences=len(vr),
            validation_batch1_cores=len({v['core_key'] for v in vr}),
            validation_batch1_assemblies=len({v['sample'] for v in vr}))
    a.repo_out.mkdir(parents=True, exist_ok=True)
    core_rows.sort(key=lambda c: (c['structure'], c['identical_core_dna_pairs'] == 0, -c['assemblies_with_local_reads'],
                                  -c['pairs_match_both_lt2500'], c['core_key']))
    for name, rows in (('core_candidates', core_rows), ('match_both_lt2500_pairs', good_pairs),
                       ('validation_batch1', val_rows)):
        pc.write_tsv(a.repo_out / ('%s_%s.tsv.gz' % (name, tag)), rows, list(rows[0]) if rows else ['structure'])
    (a.repo_out / ('candidates_summary_%s.json' % tag)).write_text(json.dumps(out, indent=1) + '\n')
    print(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
