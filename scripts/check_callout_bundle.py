"""Check one assembly call-out plus v1 read evidence, without changing discovery.

Live mode reruns package coordinate/sequence/hash validation. --snapshot checks
archived small tables only and can never certify sequence or read correctness.
Neither mode implements the planned audited CIGAR/competing-locus classifier.
"""

import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path


def require(ok, message):
    if not ok:
        raise ValueError(message)


def load(path):
    return json.loads(Path(path).read_text())


def check_tables(summary, validation, families, bed, evidence, read_summary):
    """Fail on broken joins/counts; zero biological yield is not an error."""
    require(summary['subset'] == 'all', 'Scale contract requires --subset all')
    require(validation['status'] == 'pass', 'Package validation did not pass')
    require(validation['native_outputs_unchanged'] is True, 'Native hashes not validated')
    require(len({f['query'] for f in families}) == len(families), 'Duplicate family IDs')
    require(len(families) == summary['native_families'] == summary['selected_families'], 'Family count mismatch')
    require(sum(f['raw_occurrences'] for f in families) == summary['native_records'] == summary['selected_native_records'], 'Native inventory count mismatch')
    require(summary['exported_records'] + summary['invalid_records'] == summary['selected_native_records'], 'Export/invalid count mismatch')
    require(validation['exported_records'] == summary['exported_records'], 'Validation export count mismatch')
    require(validation['invalid_records_retained'] == summary['invalid_records'], 'Validation invalid count mismatch')
    require(sum(summary['boundary_sources'].values()) == summary['exported_records'], 'Boundary totals mismatch')
    require(len(bed) == summary['exported_records'], 'BED count mismatch')
    bm = {r[3]: r for r in bed}
    require(len(bm) == len(bed), 'Duplicate BED member IDs')
    em = {r['member']: r for r in evidence}
    require(len(em) == len(evidence), 'Duplicate evidence member IDs')
    require(set(bm) == set(em), 'Missing/extra read-evidence members')
    fm = {f['query']: f for f in families}
    expected_members = {f"{f['query']}_m{i}" for f in families for i in range(1, f['raw_occurrences'] + 1)}
    require(set(bm) <= expected_members, 'Unknown native member ID')
    for member, r in em.items():
        b = bm[member]
        require(r['family'] in fm, member + ': unknown family')
        f = fm[r['family']]
        require(member in {f"{f['query']}_m{i}" for i in range(1, f['raw_occurrences'] + 1)}, member + ': wrong family')
        require(r['method'] == f['method'], member + ': wrong method')
        require(r['tnp_supported'] in ('True', 'False'), member + ': invalid tnp boolean')
        require((r['tnp_supported'] == 'True') == bool(f['tnp_supported']), member + ': tnp label mismatch')
        require((r['contig'], int(r['start0']), int(r['end0'])) == (b[0], int(b[1]), int(b[2])), member + ': BED/evidence coordinates differ')
        require(0 <= int(r['start0']) < int(r['end0']), member + ': invalid interval')
        require(int(r['length']) == int(r['end0']) - int(r['start0']), member + ': wrong length')
        require(all(int(r[k]) >= 0 for k in ('full', 'left', 'right', 'empty', 'ambiguous', 'internal')), member + ': negative support')
    require(len(evidence) == read_summary['members'], 'Read summary count mismatch')
    require(Counter(r['class'] for r in evidence) == Counter({k:v for k,v in read_summary['member_classes'].items() if v}), 'Read class totals mismatch')
    require(Counter(r['boundary_source'] for r in evidence) == Counter(summary['boundary_sources']), 'Evidence boundary labels mismatch')
    return em


def summarize(evidence, boundaries=None, triplets=None):
    """Record counts; exact coordinate groups are diagnostic, not event counts."""
    by_family = defaultdict(list)
    for r in evidence:
        by_family[r['family']].append(r)
    rows = []
    for family, members in sorted(by_family.items()):
        long = [r for r in members if int(r['length']) >= 400]
        adequate = [r for r in long if int(r['full']) >= 2]
        eligible = None
        if boundaries is not None:
            eligible = sum(
                boundaries[r['member']]['flanks_complete']
                and r['boundary_source'] in ('deposited_mobile_element', 'provisional_homology_median')
                and set(str(triplets[r['member']].seq).upper()) <= set('ACGT')
                for r in adequate
            )
        rows.append(dict(family=family, method=members[0]['method'], tnp_supported=members[0]['tnp_supported'] == 'True',
                         exported_records=len(members), long_records=len(long),
                         exact_coordinate_groups=len({(r['contig'], int(r['start0']), int(r['end0'])) for r in members}),
                         long_v1_full_ge2_records=len(adequate),
                         long_v1_full_ge2_with_empty_signal_records=sum(int(r['empty']) > 0 for r in adequate),
                         long_v1_full_ge2_complete_unambiguous_triplet_records=eligible))
    return rows


def audit(sample_dir, reads_verify, snapshot=False):
    sample_dir = Path(sample_dir)
    callout = sample_dir / 'callout'
    results = callout if snapshot else callout / 'results'
    if snapshot:
        validation = load(results / 'validate.json')
    else:
        from multisite_insertion.validation import validate
        validation = validate(callout)
    paths = [results/'summary.json', results/'all_families_annotated.json', results/'candidates.bed',
             sample_dir/'evidence/member_read_support.tsv', sample_dir/'evidence/summary.json', Path(reads_verify)]
    summary, families = load(paths[0]), load(paths[1])
    with paths[2].open() as f:
        bed = [l.rstrip('\n').split('\t') for l in f if l.strip() and not l.startswith('#')]
    with paths[3].open() as f:
        evidence = list(csv.DictReader(f, delimiter='\t'))
    read_summary, verify = load(paths[4]), load(paths[5])
    check_tables(summary, validation, families, bed, evidence, read_summary)
    require(verify['status'] == 'verified', 'Reads download verification failed')
    for key in ('bytes_match', 'md5_match', 'base_count_match'):
        require(verify['checks'].get(key) is True, 'Reads verification missing/failed: ' + key)
    sample = read_summary['sample']
    require(verify['run']['sample_id'] == sample, 'Reads verification belongs to a different sample')
    require(all(f['organism'] == sample for f in families), 'Mixed/different samples in callout')
    boundaries, triplets = None, None
    if not snapshot:
        from multisite_insertion.io import fasta_dict
        staged = load(callout/'inputs/manifest.json')['samples']
        require(len(staged) == 1 and staged[0]['id'] == sample, 'Live check requires exactly this one sample')
        boundaries = {r['member']:r for r in load(results/'member_boundaries.json')}
        triplets = fasta_dict(callout/'sequences/triplet.fasta')
        require(set(boundaries) == {r['member'] for r in evidence}, 'Boundary/evidence IDs mismatch')
        for r in evidence:
            b = boundaries[r['member']]
            require((b['candidate_start0'], b['candidate_end0'], b['accession'], b['boundary_source']) ==
                    (int(r['start0']), int(r['end0']), r['contig'], r['boundary_source']), 'Boundary/evidence version mismatch')
    report = dict(schema_version=1, sample=sample,
                  status='snapshot_tables_consistent' if snapshot else 'live_exports_and_v1_tables_consistent',
                  scope='Table joins/counts and archived checks only' if snapshot else 'Package export validation plus v1 table joins/counts; not raw read remapping',
                  input_sha256={str(p.relative_to(sample_dir)) if p.is_relative_to(sample_dir) else str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
                  package_validation=validation, summary=summary,
                  family_rows=summarize(evidence, boundaries, triplets),
                  boundary_source_counts=dict(Counter(r['boundary_source'] for r in evidence)),
                  short_records=sum(int(r['length']) < 400 for r in evidence),
                  long_records=sum(int(r['length']) >= 400 for r in evidence),
                  read_evidence_version='v1', audited_read_check='not_performed',
                  independent_molecule_count='not_established_by_this_checker',
                  scale_release_approved=False,
                  limitations=['No raw read/PAF audit or independent molecule verification',
                               'Candidate boundaries are not certified element ends',
                               'Exact coordinate groups do not resolve overlapping calls or cross-strain contexts',
                               'Empty sites and tnp are not eligibility gates'])
    return report


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--sample-dir', type=Path, required=True)
    p.add_argument('--reads-verify', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--snapshot', action='store_true')
    a = p.parse_args()
    try:
        require(not a.out.exists(), 'Output already exists; choose a new report path')
        report = audit(a.sample_dir, a.reads_verify, a.snapshot)
        a.out.parent.mkdir(parents=True, exist_ok=True)
        with a.out.open('x') as f:
            json.dump(report, f, indent=2)
            f.write('\n')
        print(json.dumps({k:report[k] for k in ('sample','status','long_records','short_records','scale_release_approved')}))
    except (ValueError, OSError, KeyError, TypeError) as e:
        p.exit(1, f'Error: {e}\n')


if __name__ == '__main__':
    main()
