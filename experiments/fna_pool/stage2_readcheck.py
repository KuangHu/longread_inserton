"""Scale-up stage 2, batch 1: known-site long-read check of stage-1 candidates with reads already on local disk.

Validates known sites only (carrier present / partner absent at the called site, s3_readcheck v0.3 unchanged); it does
not search reads for new sites with the cores. Read evidence is a new layer next to the stage-1 assembly evidence,
which stays as it is (categories, lists and IDs of the stage-1 catalog are not changed). No reads are downloaded.

select      Frozen targets from the stage-1 catalog: main loci whose carrier has at least one usable run that is
            already verified on local disk (read_batch_manifest). One pool batch per run of the read check.
localcheck  For every sample of the tests (carriers and partners): its usable linked runs, each checked on disk
            (verify.json says verified, same BioSample and platform, file present, size and SHA256 recomputed equal to
            the record). Runs that pass go to runs_allowed.tsv (the only runs mapped); missing, failed or remote runs
            are listed with their reason and never count as evidence against a candidate.
layer       After s3_readcheck summarize: read_evidence.tsv (stage-1 schema, one row per test x role), per-category
            report counts, comparison with the baseline pilot tests, anomalies, and the next batch's deduplicated
            download list (metadata only).
"""

import argparse
import csv
import gzip
import hashlib
import json
import os
import subprocess
import sys
import time
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import s3_readcheck as rc3  # noqa: E402

VERSION = "stage2_readcheck_v1"
POOL = "scaleup_stage1_20261026"
CATALOG = HERE / POOL / "stage1_catalog"
BASELINE = HERE / "comamonadaceae_20261024" / "s5_readcheck_main"


def read(path):
    op = gzip.open if str(path).endswith(".gz") else open
    with op(str(path), "rt") as h:
        return list(csv.DictReader(h, delimiter="\t"))


def write(path, rows, cols=None):
    cols = cols or (list(rows[0]) if rows else ["id"])
    with open(str(path), "w") as h:
        w = csv.DictWriter(h, cols, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def split(x):
    return [v for v in x.split(";") if v] if x else []


def sha256(path):
    h = hashlib.sha256()
    with open(str(path), "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def cmd_select(a):
    a.work.mkdir(parents=True, exist_ok=True)
    core = {x["global_locus_id"]: x for x in read(CATALOG / "candidate_core_catalog.tsv.gz") if x["category"] == a.category}
    rows = []
    for m in read(CATALOG / "read_batch_manifest.tsv.gz"):
        if m["role"] == "carrier" and m["category"] == a.category and m["local_verified_runs"] and m["batch"] == a.batch:
            x = core[m["global_locus_id"]]
            rows.append(dict(locus_id=x["locus_id"], global_locus_id=x["global_locus_id"], batch=x["batch"], category=x["category"],
                             carrier=x["biosample"], assembly=x["assembly"], contig=x["contig"], core_start=x["core_start"],
                             core_end=x["core_end"], core_sha256=x["core_sha256"], carrier_local_runs=m["local_verified_runs"],
                             carrier_usable_runs=m["usable_runs"]))
    other = Counter(m["batch"] for m in read(CATALOG / "read_batch_manifest.tsv.gz")
                    if m["role"] == "carrier" and m["category"] == a.category and m["local_verified_runs"] and m["batch"] != a.batch)
    rows.sort(key=lambda r: r["global_locus_id"])
    write(a.work / "targets.tsv", rows)
    frozen = dict(version=VERSION, created=time.strftime("%Y-%m-%dT%H:%M:%S%z"), catalog_commit=a.catalog_commit, batch=a.batch,
                  category=a.category, rule="carrier has >= 1 usable run verified on local disk (stage-1 read_batch_manifest)",
                  targets=len(rows), same_rule_other_batches=dict(other),
                  inputs={p.name: sha256(p) for p in (CATALOG / "candidate_core_catalog.tsv.gz", CATALOG / "read_batch_manifest.tsv.gz")})
    (a.work / "targets.json").write_text(json.dumps(frozen, indent=1) + "\n")
    print(json.dumps(frozen, indent=1))


def cmd_localcheck(a):
    tests = read(a.work / "tests.tsv")
    roles = defaultdict(set)
    for t in tests:
        roles[t["carrier"]].add("carrier")
        roles[t["empty"]].add("partner")
    runs = {r["run_accession"]: r for r in read(a.pool_dir / "runs.tsv")}
    linked = defaultdict(set)
    for l in read(a.pool_dir / "assembly_run_links.tsv"):
        if l["biosample"] in roles:
            linked[l["biosample"]].add(l["run"])
    items = [(s, r) for s in sorted(roles) for r in sorted(linked[s])]

    def one(item):
        s, run = item
        r = runs.get(run, {})
        row = dict(sample=s, roles=";".join(sorted(roles[s])), run=run, platform=r.get("instrument_platform", ""),
                   read_type=r.get("read_type", ""), usable=r.get("usable", ""), unusable_reasons=r.get("unusable_reasons", ""),
                   fastq_bytes=r.get("fastq_bytes", ""), preset=rc3.PRESET.get(r.get("read_type", ""), ""))
        vp = a.reads / run / "verify.json"
        if r.get("usable") != "True":
            return dict(row, status="linked_not_usable")
        if not vp.exists():
            return dict(row, status="not_local")
        v = json.loads(vp.read_text())
        lf = Path(v.get("local_file", ""))
        row.update(local_file=str(lf), verify_status=v.get("status", ""), verify_sample=v["run"].get("sample_id", ""),
                   verify_platform=v["run"].get("platform", ""), bytes_record=v.get("bytes", ""), md5_record=v.get("md5", ""),
                   sha256_record=v.get("sha256", ""), repair=json.dumps(v["repair"]) if v.get("repair") else "")
        if v.get("status") != "verified":
            return dict(row, status="verify_record_not_verified")
        if v["run"].get("sample_id") != s or v["run"].get("platform") != r.get("instrument_platform"):
            return dict(row, status="verify_record_sample_or_platform_differs")
        if not lf.exists():
            return dict(row, status="local_file_missing")
        row["bytes_on_disk"] = lf.stat().st_size
        if row["bytes_on_disk"] != v["bytes"]:
            return dict(row, status="local_size_mismatch")
        row["sha256_now"] = sha256(lf)
        return dict(row, status="local_verified_ok" if row["sha256_now"] == v["sha256"] else "local_sha256_mismatch")

    t0 = time.time()
    with ThreadPoolExecutor(max_workers=a.jobs) as ex:
        rows = list(ex.map(one, items))
    write(a.work / "runs_checked.tsv", rows, ["sample", "roles", "run", "platform", "read_type", "preset", "usable", "unusable_reasons",
                                              "fastq_bytes", "status", "local_file", "verify_status", "verify_sample", "verify_platform",
                                              "bytes_record", "bytes_on_disk", "md5_record", "sha256_record", "sha256_now", "repair"])
    write(a.work / "runs_allowed.tsv", [r for r in rows if r["status"] == "local_verified_ok"], ["run", "sample", "read_type", "preset"])
    summ = dict(samples=len(roles), runs=len(rows), status=dict(Counter(r["status"] for r in rows)),
                bytes_rehashed=sum(int(r.get("bytes_on_disk") or 0) for r in rows if r.get("sha256_now")), seconds=round(time.time() - t0))
    (a.work / "runs_checked.json").write_text(json.dumps(summ, indent=1) + "\n")
    print(json.dumps(summ, indent=1))


def report_class(t):
    """Report category of one test from the v0.3 verdicts (no new thresholds)."""
    v = t["test_verdict"]
    if v != "not_tested":
        return v
    sides = [(t["C_verdict"], int(t["C_locus_ambiguous"])), (t["E_verdict"], int(t["E_locus_ambiguous"]))]
    if all(s == "no_reads" for s, _ in sides):
        return "untested_no_local_reads"
    if any(s == "insufficient" and la > 0 for s, la in sides):
        return "insufficient_with_repeat_placement_unclear"
    return "insufficient_coverage"


def cmd_layer(a):
    out = a.out
    tests = read(out / "readcheck_tests.tsv")
    loci = read(out / "readcheck_loci.tsv")
    targets = {r["locus_id"]: r for r in read(a.work / "targets.tsv")}
    runs_checked = read(a.work / "runs_checked.tsv")
    allowed = defaultdict(list)
    for r in runs_checked:
        if r["status"] == "local_verified_ok":
            allowed[r["sample"]].append(r["run"])
    used = {}
    for p in (a.work / "paf").glob("genome_*.json"):
        d = json.loads(p.read_text())
        used[d["sample"]] = d["runs"]
    created = time.strftime("%Y-%m-%dT%H:%M:%S%z")
    lev = {l["locus_id"]: l["evidence_level"] for l in loci}
    ev = []
    for t in tests:
        g = targets[t["locus_id"]]["global_locus_id"]
        for role, s, exp in (("C", t["carrier"], "present"), ("E", t["empty"], "absent")):
            rs = used.get(s, [])
            ev.append(dict(evidence_id="RE1_%s_%s" % (t["test_id"], role), global_locus_id=g, test_id=t["test_id"],
                           run=";".join(r["run"] for r in rs), sample_key="%s:%s" % (a.batch, s), role="carrier" if role == "C" else "partner",
                           caller_version=rc3.VERSION, count_unit=t["%s_count_unit" % role], present=t["%s_present" % role],
                           absent=t["%s_absent" % role], present_matched=t["%s_present_matched" % role],
                           absent_matched=t["%s_absent_matched" % role], ambiguous=t["%s_ambiguous" % role],
                           locus_ambiguous=t["%s_locus_ambiguous" % role], off_locus=t["%s_off_locus" % role],
                           verdict=t["%s_verdict" % role], evidence_level=lev[t["locus_id"]],
                           conflict=t["%s_verdict" % role] in ("mixed", "contradicted"), created=created,
                           expected_allele=exp, test_verdict=t["test_verdict"], report_class=report_class(t)))
    base_cols = open(str(CATALOG / "read_evidence.tsv")).readline().rstrip("\n").split("\t")
    write(out / "read_evidence.tsv", ev, base_cols + ["expected_allele", "test_verdict", "report_class"])
    # ---- counts
    rc = Counter(report_class(t) for t in tests)
    by_locus = defaultdict(list)
    for t in tests:
        by_locus[t["locus_id"]].append(report_class(t))
    samples_c = {t["carrier"] for t in tests}
    samples_e = {t["empty"] for t in tests}
    with_reads = {s for s, rs in used.items() if rs}
    runs_used = {r["run"] for rs in used.values() for r in rs}
    summ = dict(
        version=VERSION, caller=rc3.VERSION, evidence_source_added="long_read_known_site_check", batch=a.batch,
        targets=len(targets), loci_tested=len(loci), tests=len(tests),
        test_report_class=dict(rc), test_verdicts=dict(Counter(t["test_verdict"] for t in tests)),
        side_verdicts=dict(Counter("%s/%s" % (t["C_verdict"], t["E_verdict"]) for t in tests)),
        locus_evidence_levels=dict(Counter(l["evidence_level"] for l in loci)),
        loci_any_test_class={c: sum(1 for v in by_locus.values() if c in v) for c in sorted(rc)},
        samples=dict(carriers=len(samples_c), partners=len(samples_e - samples_c), both_roles=len(samples_c & samples_e),
                     carriers_with_local_reads=len(samples_c & with_reads), partners_with_local_reads=len((samples_e - samples_c) & with_reads)),
        runs=dict(checked=len(runs_checked), status=dict(Counter(r["status"] for r in runs_checked)), used=len(runs_used),
                  used_by_read_type=dict(Counter(r["read_type"] for r in runs_checked if r["run"] in runs_used))),
        molecules_excluded=dict(locus_ambiguous=sum(int(t["C_locus_ambiguous"]) + int(t["E_locus_ambiguous"]) for t in tests),
                                off_locus=sum(int(t["C_off_locus"]) + int(t["E_off_locus"]) for t in tests),
                                ambiguous=sum(int(t["C_ambiguous"]) + int(t["E_ambiguous"]) for t in tests),
                                cigar_fail=sum(int(t["C_cigar_fail"]) + int(t["E_cigar_fail"]) for t in tests)))
    # ---- anomalies
    anom = []
    for t in tests:
        k = report_class(t)
        if k in ("mixed", "contradicted", "insufficient_with_repeat_placement_unclear"):
            anom.append(dict(kind=k, test_id=t["test_id"], locus_id=t["locus_id"], global_locus_id=targets[t["locus_id"]]["global_locus_id"],
                             carrier=t["carrier"], partner=t["empty"], C_verdict=t["C_verdict"], E_verdict=t["E_verdict"],
                             C_counts="P%s/A%s matched P%s/A%s amb%s lamb%s off%s" % tuple(t["C_" + x] for x in (
                                 "present", "absent", "present_matched", "absent_matched", "ambiguous", "locus_ambiguous", "off_locus")),
                             E_counts="P%s/A%s matched P%s/A%s amb%s lamb%s off%s" % tuple(t["E_" + x] for x in (
                                 "present", "absent", "present_matched", "absent_matched", "ambiguous", "locus_ambiguous", "off_locus")),
                             difference_class=t["difference_class"], core_len=t["core_len"]))
    for r in runs_checked:
        if r["status"] not in ("local_verified_ok", "not_local", "linked_not_usable"):
            anom.append(dict(kind="run_" + r["status"], carrier=r["sample"], test_id="", locus_id="", global_locus_id="",
                             partner="", C_verdict="", E_verdict="", C_counts=r["run"], E_counts="", difference_class="", core_len=""))
    write(out / "anomalies.tsv", anom, ["kind", "test_id", "locus_id", "global_locus_id", "carrier", "partner", "C_verdict", "E_verdict",
                                        "C_counts", "E_counts", "difference_class", "core_len"])
    summ["anomalies"] = dict(Counter(x["kind"] for x in anom))
    # ---- comparison with the baseline pilot (exactly linked tests only)
    links = {t["test_id"]: t for t in read(CATALOG / "pilot_read_test_links.tsv") if t["link_scope"] != "not_linked"}
    pilot = {t["test_id"]: t for t in read(BASELINE / "readcheck_tests.tsv")}
    pruns = {}
    for p in BASELINE.glob("genome_*.json"):
        d = json.loads(p.read_text())
        pruns[d["sample"]] = sorted(r["run"] for r in d["runs"])
    new = {t["test_id"]: t for t in tests}
    geo = ["C_site", "E_site", "core_len", "other_side_gap", "aln_strand", "CP_len", "CA_len", "EA_len", "EP_len"]
    cnt = [r + x for r in "CE" for x in ("_present", "_absent", "_present_matched", "_absent_matched", "_ambiguous", "_locus_ambiguous",
                                         "_off_locus", "_cigar_fail")]
    comp = []
    for tid, l in sorted(links.items()):
        p, n = pilot[tid], new.get(tid)
        row = dict(test_id=tid, link_scope=l["link_scope"], pilot_locus_id=l["pilot_locus_id"], global_locus_id=l["global_locus_id"],
                   carrier=p["carrier"], partner=p["empty"], pilot_test_verdict=p["test_verdict"], pilot_C=p["C_verdict"], pilot_E=p["E_verdict"])
        if n is None:
            row.update(status="test_not_in_batch1", change="")
            comp.append(row)
            continue
        row.update(new_test_verdict=n["test_verdict"], new_C=n["C_verdict"], new_E=n["E_verdict"])
        why = []
        if any(p[k] != n[k] for k in geo):
            why.append("test_geometry_differs")
        for role, s in (("C", p["carrier"]), ("E", p["empty"])):
            if pruns.get(s, []) != sorted(r["run"] for r in used.get(s, [])):
                why.append("%s_runs_differ:%s->%s" % (role, ",".join(pruns.get(s, [])) or "none", ",".join(sorted(r["run"] for r in used.get(s, []))) or "none"))
        same_counts = all(p[k] == n[k] for k in cnt)
        if not same_counts:
            why.append("counts_differ")
        same = p["test_verdict"] == n["test_verdict"] and p["C_verdict"] == n["C_verdict"] and p["E_verdict"] == n["E_verdict"]
        row.update(status="same_verdicts" if same else "verdict_changed", counts_identical=same_counts, change=";".join(why) or "-")
        comp.append(row)
    write(out / "pilot_comparison.tsv", comp)
    summ["pilot_comparison"] = dict(linked_tests=len(comp), status=dict(Counter(c["status"] for c in comp)),
                                    by_scope={s: dict(Counter(c["status"] for c in comp if c["link_scope"] == s)) for s in sorted({c["link_scope"] for c in comp})},
                                    counts_identical=sum(1 for c in comp if c.get("counts_identical") is True),
                                    changes=dict(Counter(c["change"] for c in comp if c["status"] == "verdict_changed")))
    # ---- next batch download list (metadata only)
    rinfo = {r["run"]: r for r in read(CATALOG / "read_run_manifest.tsv.gz")}
    per = defaultdict(lambda: dict(carrier=[], partner=[]))
    for m in read(CATALOG / "read_batch_manifest.tsv.gz"):
        if m["category"] == "main":
            per[m["global_locus_id"]]["carrier" if m["role"] == "carrier" else "partner"].append(m)
    tested = {targets[t]["global_locus_id"] for t in targets}

    def usable(ms):
        return {r for m in ms for r in split(m["usable_runs"])}

    def local(ms):
        return {r for m in ms for r in split(m["local_verified_runs"])}

    tiers = defaultdict(set)
    loci_t = Counter()
    for g, d in per.items():
        if g in tested:
            continue
        cu, pu = usable(d["carrier"]), usable(d["partner"])
        if cu and pu:
            tier = "main_carrier_and_partner_usable"
        elif cu:
            tier = "main_carrier_only_usable"
        elif pu:
            tier = "main_partner_only_usable"
        else:
            continue
        loci_t[tier] += 1
        tiers[tier] |= (cu | pu) - local(d["carrier"] + d["partner"])
    nxt, seen = [], set()
    for tier in ("main_carrier_and_partner_usable", "main_carrier_only_usable", "main_partner_only_usable"):
        for r in sorted(tiers[tier] - seen):
            i = rinfo[r]
            nxt.append(dict(priority_tier=tier, run=r, batch=i["batch"], biosample=i["biosample"], platform=i["platform"], read_type=i["read_type"],
                            fastq_ftp=i["fastq_ftp"], fastq_bytes=i["fastq_bytes"], fastq_md5=i["fastq_md5"],
                            size_known=any(v.isdigit() for v in i["fastq_bytes"].split(";"))))
            seen.add(r)
    write(out / "next_batch_download_list.tsv", nxt)

    def nb(rows):
        return sum(sum(int(v) for v in r["fastq_bytes"].split(";") if v.isdigit()) for r in rows)

    summ["next_batch"] = {t: dict(loci=loci_t[t], new_runs_in_tier=sum(1 for r in nxt if r["priority_tier"] == t),
                                  bytes_known=nb([r for r in nxt if r["priority_tier"] == t]),
                                  runs_unknown_size=sum(1 for r in nxt if r["priority_tier"] == t and not r["size_known"]))
                          for t in ("main_carrier_and_partner_usable", "main_carrier_only_usable", "main_partner_only_usable")}
    summ["next_batch_note"] = "metadata only, nothing downloaded; each run counted once, in its highest tier"
    (out / "batch1_summary.json").write_text(json.dumps(summ, indent=1) + "\n")
    print(json.dumps(summ, indent=1))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("mode", choices=["select", "localcheck", "layer"])
    p.add_argument("--work", type=Path, required=True)
    p.add_argument("--out", type=Path)
    p.add_argument("--pool-dir", type=Path)
    p.add_argument("--reads", type=Path, default=Path("/global/scratch/users/kh36969/longread_callout/data/reads"))
    p.add_argument("--batch", default="comamonadaceae")
    p.add_argument("--category", default="main")
    p.add_argument("--catalog-commit", default="")
    p.add_argument("--jobs", type=int, default=8)
    a = p.parse_args()
    dict(select=cmd_select, localcheck=cmd_localcheck, layer=cmd_layer)[a.mode](a)


if __name__ == "__main__":
    main()
