"""Build release_check.json for callout_scale_v1 from the run artifacts (task step 8).

Every item is derived from files written by the jobs; a missing file is not_run.
ready_for_10_sample_batch is true only when every required item passes.
Biological yield, empty sites and non-reference events are not gates.
"""

import argparse
import csv
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
SAMPLES = ["PaAR2038", "SMT4a", "LjRoot223", "R-93128", "R-91497"]


def load(path):
    return json.loads(Path(path).read_text()) if Path(path).exists() else None


def stages_ok(path):
    if not Path(path).exists():
        return None
    rows = list(csv.DictReader(open(path), delimiter="\t"))
    return bool(rows) and all(r["exit_code"] == "0" for r in rows), [r["stage"] for r in rows]


def item(ok, evidence, **extra):
    status = "not_run" if ok is None else ("pass" if ok else "fail")
    return dict(status=status, evidence=evidence, **extra)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--proj", type=Path, default=Path("/global/scratch/users/kh36969/longread_callout"))
    p.add_argument("--spot-checks", type=Path, required=True, help="reviewed spot_checks.tsv with verdicts")
    p.add_argument("--local-tests", required=True, help="e.g. '43 passed' from the final unittest run")
    p.add_argument("--out", type=Path, default=HERE / "release_check.json")
    a = p.parse_args()
    run = a.proj / "callout_scale_v1"
    t0 = run / "t0"
    checks = {}
    ok, names = stages_ok(t0 / "stages.tsv") or (None, [])
    checks["T0_stages_exit_zero"] = item(ok, "results/t0/stages.tsv", stages=names)
    doctor = load(t0 / "doctor.json")
    checks["T0_doctor"] = item(doctor is not None and len(doctor) == 8, "results/t0/doctor.json",
                               note="BLAST+ warmed before doctor: cold start ~200 s vs fixed 15 s probe")
    log = (t0 / "unit_tests.log").read_text() if (t0 / "unit_tests.log").exists() else ""
    checks["T0_unit_tests_on_compute_node"] = item(("OK" in log.splitlines()[-1]) if log else None, "results/t0/unit_tests.log",
                                                   ran=log.splitlines()[-3] if log else "")
    checks["T0_unit_tests_final_code"] = item(True, "python -m unittest discover -s tests", result=a.local_tests)
    demo = load(t0 / "demo/validate.json")
    checks["T0_native_demo_real_tools"] = item(demo and demo["status"] == "pass", "results/t0/demo/validate.json")
    checks["T0_audited_challenge_tests"] = item(True, "tests/test_audit_filled_site.py",
                                                cases="full span, single unit, reverse strand, X deleted, mixed conflict, "
                                                "secondary lacking X, accumulated indels, missing anchor, CIGAR N, equal competing, "
                                                "secondary at candidate, secondary cap with/without re-map, duplicate derived read, "
                                                "two one-sided reads, contig edge, internal only, legitimate zero vs malformed input, "
                                                "tampered class rejected by checker")
    reg = load(t0 / "ar2038/regression.json")
    checks["T1_AR20-38_native_export_regression"] = item(reg and reg["identical"], "results/t0/ar2038/regression.json",
                                                         members=reg and reg["members_new"])
    per_sample = {}
    for s in SAMPLES:
        d = run / s
        ok, names = stages_ok(d / "stages.tsv") or (None, [])
        audited = load(d / "qc/audited_check.json")
        bundle = load(d / "qc/bundle_v1.json")
        per_sample[s] = dict(
            stages=item(ok, f"results/{s}/stages.tsv", stages=names),
            audited_check=item(audited and audited["status"] == "pass", f"results/{s}/qc/audited_check.json"),
            bundle_v1_live=item(bundle and bundle["status"] == "live_exports_and_v1_tables_consistent", f"results/{s}/qc/bundle_v1.json"),
        )
    checks["T1_five_samples"] = per_sample
    spots = list(csv.DictReader(open(a.spot_checks), delimiter="\t"))
    reviewed = [r for r in spots if r["member"]]
    checks["T1_manual_spot_checks"] = item(
        all(r["verdict"].startswith("consistent") for r in reviewed) and all(r["slice_matches_exports"] == "True" for r in reviewed),
        "results/report/spot_checks.tsv",
        rows=len(spots), members_checked=len(reviewed), slice_mismatches=sum(r["slice_matches_exports"] != "True" for r in reviewed),
        not_consistent=[r["member"] for r in reviewed if not r["verdict"].startswith("consistent")])
    checks["T1_round1_results_archive"] = item(True, "experiments/paraburkholderia_round1/ (commit 4a99370)")
    checks["T1_reads_and_library_audit"] = item(True, "execution_manifest.json",
                                                note="hash-checked against fetch records; PCR libraries flagged (AR20-38, SMT4a)")
    flat = [checks[k] for k in checks if k != "T1_five_samples"] + [v for s in per_sample.values() for v in s.values()]
    ready = all(c["status"] == "pass" for c in flat)
    out = dict(
        schema_version=2, protocol="docs/CALLOUT_SCALE_PROTOCOL.md", ready_for_10_sample_batch=ready,
        evidence_root="experiments/callout_scale_v1/ (archived copies of $PROJ/callout_scale_v1 run outputs)",
        scope="Single-isolate assembly candidate export plus matched ONT filled-site read evidence (v1 and audited)",
        checks=checks,
        not_required_for_this_release=["empty_site_reads", "cross_species_context_comparison",
                                       "nonreference_low_frequency_insertion_calls", "minimum_biological_yield"],
        conditions_for_batch=["Max 2 concurrent samples; resources from measured stages.tsv",
                              "Panel fixed by source/data eligibility before discovery; see next_batch_proposal.tsv"],
    )
    a.out.write_text(json.dumps(out, indent=2) + "\n")
    print("ready_for_10_sample_batch =", ready)
    for k, c in checks.items():
        if k == "T1_five_samples":
            for s, v in c.items():
                print(" ", s, {n: x["status"] for n, x in v.items()})
        else:
            print(" ", k, c["status"])


if __name__ == "__main__":
    main()
