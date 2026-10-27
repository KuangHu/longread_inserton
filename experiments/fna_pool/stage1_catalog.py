"""Scale-up stage 1, section G: cross-batch candidate catalog, unified protein families, stage-2 interface manifests.

Evidence source of everything here is assembly_only: assembly-derived candidate loci, long-read evidence pending. A
main locus is not an independent insertion event and a family is a preliminary protein family, not an insertion family.

Namespaces: batch = the family directory name under the pool (e.g. nitrobacteraceae). Locus, sample and local family
IDs from one batch are written as <batch>:<id> (global_locus_id, sample_key, local_family_id) so equal IDs from
different batches never overwrite each other. BioSample accessions are kept as is next to the namespaced key.

proteins  Proteins of complete CDS in the candidate lists (main, partner_mixed, direction_mixed, short_other_segment,
          replacement) of every batch, read from each carrier's FAA, deduplicated by the SHA256 of the full sequence
          (protein accessions != distinct sequences). Writes <work>/proteins.faa (named by sha256) and
          protein_instances.tsv.gz (one row per locus x coding unit, with the local hash and local family/homology IDs).
          Checks that each local 16-hex protein hash is the prefix of the global sha256.
cluster   MMseqs2 on the frozen global protein set at the two frozen S4 settings (s4_annotate.MMSEQS_SETTINGS, with
          createdb --shuffle 1 --createdb-mode 1 and --alignment-mode 3 exactly as easy-cluster runs them), then every
          member aligned to its final representative (mmseqs align on the cluster DB, alignment-mode 3, -e 10). This
          gives the same clusters and member values as the per-batch all-vs-all check (tested on Eggerthellaceae).
catalog   Global family / homology-group IDs bound to the catalog version and the representative's sequence hash
          (F1_<rep sha256[:12]>, H1_<...>); architecture = global family IDs in core order with orientation, canonical
          over both strands, ID A1_<sha256(architecture)[:12]>. Local-to-global family map per batch. Candidate core
          catalog (lists + annotation_pending) with U/core/D, coordinates, hashes, contig length/circularity, repeat
          and contig-end flags; core sequences in a FASTA on scratch. Coding anomalies (translation records other than
          match/pseudo, transl_table problems) per unit. Stage-2 manifests: assembly_run_links, read_run_manifest,
          read_batch_manifest (metadata only, nothing is downloaded), empty read_evidence / new_read_contexts with
          their schema. Pilot read evidence (Comamonadaceae baseline read check) linked by exact carrier/contig/core
          coordinates (overlapping but not identical stage-1 cores are listed separately, without transferring the
          evidence); it is not propagated to families or other loci.
report    stage1_report.json: per-batch and global counts gathered from the S1 / S2 / S3 / S4 / catalog outputs.
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
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from s3_pair_calls import read_fasta  # noqa: E402
from s4_annotate import CONTEXT, LISTS, MEMBER_MIN, MMSEQS_SETTINGS, PENDING, read, read_faa, sha256, write  # noqa: E402

VERSION = "stage1_catalog_v1"
CATALOG = "1"  # catalog version digit in family / homology / architecture IDs
EVIDENCE = "assembly_only"
READ_STATUS = "long_read_evidence_pending"
HERE = Path(__file__).resolve().parent
PROJ = Path("/global/scratch/users/kh36969/longread_callout")
POOL = "scaleup_stage1_20261026"
BASELINE = HERE / "comamonadaceae_20261024"
ALIGN = ["--alignment-mode", "3", "-e", "10", "-a"]
CODING_OK = {"match", "match_with_transl_except", "pseudo"}
PREFIX = {"family": "F", "homology_group": "H"}


def batches(a):
    return a.batches or sorted(p.name for p in (HERE / POOL).iterdir() if (p / "s2" / "catalog.tsv").exists())


def nbytes(x):
    return sum(int(v) for v in x.split(";") if v.isdigit())


def s4dir(a, b):
    return a.scratch / b / "s4_s3_full"


def sha(s):
    return hashlib.sha256(s.encode()).hexdigest()


def gl(b, x):
    return "%s:%s" % (b, x)


def split(x):
    return [v for v in x.split(";") if v] if x else []


def cmd_proteins(a):
    a.work.mkdir(parents=True, exist_ok=True)
    seqs, inst = {}, []
    for b in batches(a):
        s4 = s4dir(a, b)
        files = {}
        for p in read(a.scratch / b / "s3_full" / "pairs.tsv"):
            for k in "ab":
                files[p["sample_" + k]] = p["gff_" + k]
        lfam = {r["protein_hash"]: r for r in read(s4 / "protein_members.tsv")}
        faas = {}
        for r in read(s4 / "candidate_loci_families.tsv"):
            if r["carrier"] not in faas:
                faas[r["carrier"]] = read_faa(Path(files[r["carrier"]]).parent / "protein.faa")
            for i, (pid, strand, fid, hid) in enumerate(zip(split(r["protein_ids"]), r["strands"].split(";"),
                                                             r["family_ids"].split(";"), r["homology_group_ids"].split(";"))):
                s = faas[r["carrier"]][pid]
                h = sha(s)
                lh = h[:16]
                if lh not in lfam or lfam[lh]["family_id"] != fid:
                    sys.exit("local protein hash/family mismatch %s %s %s" % (b, r["locus_id"], pid))
                seqs[h] = s
                inst.append(dict(batch=b, global_locus_id=gl(b, r["locus_id"]), locus_id=r["locus_id"], sample_key=gl(b, r["carrier"]),
                                 biosample=r["carrier"], category=r["category"], unit_order=i + 1, strand=strand, protein_id=pid,
                                 protein_sha256=h, length=len(s), local_protein_hash=lh, local_family_id=gl(b, fid),
                                 local_homology_id=gl(b, hid)))
        local = json.loads((s4 / "families_summary.json").read_text())["distinct_proteins"]
        got = len({i["protein_sha256"] for i in inst if i["batch"] == b})
        if got != local:
            sys.exit("%s: %d distinct proteins, S4 families_summary says %d" % (b, got, local))
        print(b, len(faas), "carriers", got, "distinct proteins", flush=True)
    with open(str(a.work / "proteins.faa"), "w") as h:
        for k, s in sorted(seqs.items()):
            h.write(">%s\n%s\n" % (k, s))
    write(a.work / "protein_instances.tsv.gz", inst)
    print(len(inst), "protein instances;", len(seqs), "distinct sequences;", len({i["protein_id"] for i in inst}), "protein accessions")


def cmd_cluster(a):
    w = a.work / "mmseqs"
    w.mkdir(parents=True, exist_ok=True)
    log = open(str(w / "mmseqs.log"), "w")
    cmds = []

    def mm(*args):
        cmd = [a.mmseqs] + [str(x) for x in args]
        cmds.append(" ".join(cmd))
        subprocess.run(cmd, check=True, stdout=log, stderr=subprocess.STDOUT)

    mm("createdb", a.work / "proteins.faa", w / "db", "--shuffle", "1", "--createdb-mode", "1")
    for level, opts in MMSEQS_SETTINGS.items():
        mm("cluster", w / "db", w / (level + "_clu"), w / ("tmp_" + level), "--threads", a.threads, "--alignment-mode", "3", *opts)
        mm("createtsv", w / "db", w / "db", w / (level + "_clu"), w / (level + "_cluster.tsv"))
        mm("align", w / "db", w / "db", w / (level + "_clu"), w / (level + "_aln"), "--threads", a.threads, *ALIGN)
        mm("convertalis", w / "db", w / "db", w / (level + "_aln"), w / (level + "_members.m8"),
           "--format-output", "query,target,fident,qcov,tcov,alnlen")
    ver = subprocess.run([a.mmseqs, "version"], stdout=subprocess.PIPE, universal_newlines=True).stdout.strip()
    (w / "cluster_run.json").write_text(json.dumps(dict(version=VERSION, mmseqs=a.mmseqs, mmseqs_version=ver, settings=MMSEQS_SETTINGS,
                                                        align=ALIGN, commands=cmds, input_sha256=sha256(a.work / "proteins.faa")), indent=1) + "\n")


def canonical(ids, strands):
    fwd = ",".join("%s%s" % (i, s) for i, s in zip(ids, strands))
    rev = ",".join("%s%s" % (i, "+" if s == "-" else "-") for i, s in reversed(list(zip(ids, strands))))
    return min(fwd, rev)


REPAIR_COLS = ["level", "member", "original_rep", "original_identity", "original_member_cov", "original_rep_cov", "action",
               "new_rep", "new_identity", "new_member_cov", "new_rep_cov", "candidates_passing"]


def reassign_failed(a, level, reps, check):
    """Members below the frozen threshold against their final representative: drop that link, align each to every
    representative of the level (exhaustive, same alignment settings as the S4 check) and move it to the best
    representative that passes (identity, then the lower coverage); none passes -> its own singleton group.
    Thresholds unchanged; the original failing values are kept in the returned records."""
    mid, mc = MEMBER_MIN[level]
    failed = sorted(m for m, r in reps.items() if not check[(level, m)][4])
    if not failed:
        return {}
    w = a.work / "mmseqs" / ("repair_" + level)
    w.mkdir(parents=True, exist_ok=True)
    seqs, name = {}, None
    for line in open(str(a.work / "proteins.faa")):
        if line.startswith(">"):
            name = line[1:].strip()
        else:
            seqs[name] = line.strip()
    with open(str(w / "failed.faa"), "w") as h:
        for m in failed:
            h.write(">%s\n%s\n" % (m, seqs[m]))
    with open(str(w / "reps.faa"), "w") as h:
        for r in sorted(set(reps.values())):
            h.write(">%s\n%s\n" % (r, seqs[r]))
    subprocess.run([a.mmseqs, "easy-search", str(w / "failed.faa"), str(w / "reps.faa"), str(w / "hits.m8"), str(w / "tmp"),
                    "--threads", str(a.threads), "--exhaustive-search", "1", "--alignment-mode", "3", "-e", "10", "--max-seqs", "100000",
                    "--format-output", "query,target,fident,qcov,tcov,alnlen"], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    ok = defaultdict(list)
    for line in open(str(w / "hits.m8")):
        q, t, fid, qc, tc, _ = line.rstrip("\n").split("\t")
        fid, qc, tc = float(fid), float(qc), float(tc)
        if fid >= mid and qc >= mc and tc >= mc and t != reps[q]:
            ok[q].append((fid, min(qc, tc), t, qc, tc))
    out = {}
    for m in failed:
        r0, f0, m0, c0, _ = check[(level, m)]
        rec = dict(level=level, member=m, original_rep=r0, original_identity=round(f0, 3), original_member_cov=round(m0, 3),
                   original_rep_cov=round(c0, 3), candidates_passing=len(ok[m]))
        if ok[m]:
            fid, _, t, qc, tc = sorted(ok[m], key=lambda x: (-x[0], -x[1], x[2]))[0]
            rec.update(action="reassigned_to_passing_representative", new_rep=t, new_identity=fid, new_member_cov=qc, new_rep_cov=tc)
        else:
            rec.update(action="singleton_no_passing_representative", new_rep=m, new_identity=1.0, new_member_cov=1.0, new_rep_cov=1.0)
        out[m] = rec
    return out


def pilot_audit(a, core_rows, cat):
    """Baseline (comamonadaceae_20261024) read-check loci against the stage-1 Comamonadaceae batch. An exact link needs
    the same BioSample, the same assembly (accession.version or FNA SHA256), contig, core coordinates and core SHA256.
    Unmatched loci get only reasons established from the files; anything else is unresolved, with overlapping stage-1
    loci listed. Tests are linked separately: a carrier-partner test applies only where that partner is a partner of
    the stage-1 locus. No evidence moves to overlapping loci, other partners, families or architectures."""
    b = "comamonadaceae"
    base, bseed = {}, {r["locus_id"]: r for r in read(BASELINE / "s4_s3_full_v02" / "seed_catalog.tsv")}
    for c in LISTS:
        for r in read(BASELINE / "s4_s3_full_v02" / ("list_%s.tsv" % c)):
            base[r["locus_id"]] = r
    bcat = {r["sample_id"]: r for r in read(BASELINE / "s2" / "catalog.tsv")}
    rl = read(BASELINE / "s5_readcheck_main" / "readcheck_loci.tsv")
    carriers = {r["carrier"] for r in rl}
    cm = [x for x in core_rows if x["batch"] == b]
    exact = {(x["biosample"], x["contig"], int(x["core_start"]), int(x["core_end"])): x for x in cm}
    contigs = set()
    with gzip.open(str(a.scratch / b / "s2" / "contig_index.tsv.gz"), "rt") as h:
        for line in h:
            f = line.split("\t")
            if f[0] in carriers:
                contigs.add((f[0], f[2]))
    s4 = defaultdict(list)
    with gzip.open(str(s4dir(a, b) / "s4_loci.tsv.gz"), "rt") as h:
        for r in csv.DictReader(h, delimiter="\t"):
            if r["carrier"] in carriers:
                s4[(r["carrier"], r["contig"])].append(r)
    s3 = defaultdict(list)
    with open(str(a.scratch / b / "s3_full" / "merged" / "loci.tsv")) as h:
        for r in csv.DictReader(h, delimiter="\t"):
            if r["carrier"] in carriers:
                s3[(r["carrier"], r["contig"])].append(r)
    pairs = {frozenset((p["sample_a"], p["sample_b"])) for p in read(a.scratch / b / "s3_full" / "pairs.tsv")}
    tests = defaultdict(list)
    for t in read(BASELINE / "s5_readcheck_main" / "readcheck_tests.tsv"):
        tests[t["locus_id"]].append(t)
    links, tlinks = [], []
    for r in rl:
        bl, c = base[r["locus_id"]], r["carrier"]
        s, e, contig = int(bl["core_start"]), int(bl["core_end"]), bl["contig"]
        bsha = bseed[r["locus_id"]]["core_sha256"]
        sc = cat.get(gl(b, c))
        same_asm = bool(sc) and (sc["assembly"] == bcat[c]["assembly"] or sc["fna_sha256"] == bcat[c]["fna_sha256"])
        x = exact.get((c, contig, s, e))
        ok = x is not None and same_asm and x["core_sha256"] == bsha
        reason, ov4, ov3 = "-", [], []
        if not ok:
            ov4 = [y for y in s4.get((c, contig), []) if int(y["core_start"]) < e and s < int(y["core_end"])]
            ov3 = [y for y in s3.get((c, contig), []) if int(y["core_start"]) < e and s < int(y["core_end"])]
            same4 = [y for y in ov4 if (int(y["core_start"]), int(y["core_end"])) == (s, e)]
            same3 = [y for y in ov3 if (int(y["core_start"]), int(y["core_end"])) == (s, e)]
            if not sc:
                reason = "carrier_not_in_stage1_pool"
            elif not same_asm:
                reason = "carrier_assembly_differs:%s->%s" % (bcat[c]["assembly"], sc["assembly"])
            elif (c, contig) not in contigs:
                reason = "contig_not_in_stage1_assembly"
            elif x is not None:
                reason = "coordinates_match_core_sha256_differs"
            elif same4:
                reason = "identical_core_coordinates_in_stage1_s4_category:" + same4[0]["category"]
            elif same3:
                reason = "identical_core_coordinates_in_stage1_s3:status=%s,core_len=%s" % (same3[0]["status"], same3[0]["core_len"])
            elif ov3:
                reason = "unresolved_overlapping_stage1_loci"
            else:
                reason = "unresolved_no_overlapping_stage1_locus"
        partners = set(split(x["partners_anchored"]) + split(x["partners_not_anchored"])) if ok else set()
        pt = tests[r["locus_id"]]
        links.append(dict(
            pilot_locus_id=r["locus_id"], pilot_set="comamonadaceae_20261024/s5_readcheck_main", biosample=c,
            pilot_assembly=bcat[c]["assembly"], pilot_fna_sha256=bcat[c]["fna_sha256"], stage1_assembly=sc["assembly"] if sc else "",
            stage1_fna_sha256=sc["fna_sha256"] if sc else "", contig=contig, core_start=s, core_end=e, pilot_core_sha256=bsha,
            pilot_category=bl["category"], evidence_level=r["evidence_level"], test_verdicts=r["test_verdicts"],
            link_basis="exact_biosample_assembly_contig_coordinates_core" if ok else "not_linked",
            global_locus_id=x["global_locus_id"] if ok else "", stage1_category=x["category"] if ok else "",
            unmatched_reason=reason,
            overlapping_stage1_candidates=";".join(gl(b, y["locus_id"]) for y in ov4 if y["category"] in LISTS + [PENDING]),
            overlapping_stage1_s4_loci=";".join("%s:%s:%s-%s" % (gl(b, y["locus_id"]), y["category"], y["core_start"], y["core_end"]) for y in ov4),
            overlapping_stage1_s3_loci=";".join("%s:%s:%s-%s" % (gl(b, y["locus_id"]), y["status"], y["core_start"], y["core_end"]) for y in ov3),
            pilot_partners=";".join(sorted({t["empty"] for t in pt})),
            pilot_partners_in_stage1_locus=";".join(sorted({t["empty"] for t in pt} & partners)),
            pilot_pairs_in_stage1_pairs=sum(1 for t in pt if frozenset((c, t["empty"])) in pairs),
            scope="this carrier locus only; evidence stays with its pilot tests and is not propagated to overlapping loci, "
                  "other partners, families or architectures"))
        for t in pt:
            same_p = ok and t["empty"] in partners
            tlinks.append(dict(test_id=t["test_id"], pilot_locus_id=r["locus_id"], pair_id=t["pair_id"], carrier=c, partner=t["empty"],
                               test_verdict=t["test_verdict"], carrier_verdict=t["C_verdict"], partner_verdict=t["E_verdict"],
                               global_locus_id=x["global_locus_id"] if ok else "",
                               partner_in_stage1_locus=same_p, pair_in_stage1_pairs=frozenset((c, t["empty"])) in pairs,
                               link_scope="same_carrier_partner_test" if same_p else
                               "same_carrier_locus_partner_not_in_stage1_locus" if ok else "not_linked"))
    return links, tlinks


def read_coverage(core_rows, mrows, rinfo):
    """Per candidate category: loci by carrier / partner usable-read availability and carrier read state; distinct runs
    and estimated bytes still to download (each run once; runs without a size listed apart). Metadata only."""
    by = defaultdict(lambda: dict(carrier=[], partner=[]))
    for m in mrows:
        by[m["global_locus_id"]]["carrier" if m["role"] == "carrier" else "partner"].append(m)

    def state(ms):
        if any(split(m["local_verified_runs"]) for m in ms):
            return "local_verified"
        if any(split(m["usable_runs"]) for m in ms):
            return "remote_downloadable"
        if any(split(m["linked_runs"]) for m in ms):
            return "linked_not_usable"
        return "no_linked_runs"

    def runs_needed(ms):
        return {r for m in ms for r in split(m["usable_runs"]) if r not in split(m["local_verified_runs"])}

    out = {}
    for c in LISTS + [PENDING, "all"]:
        xs = [x["global_locus_id"] for x in core_rows if c == "all" or x["category"] == c]
        d = Counter()
        need_c, need_p, local = set(), set(), set()
        for g in xs:
            cs, ps = state(by[g]["carrier"]), state(by[g]["partner"])
            cu, pu = cs in ("local_verified", "remote_downloadable"), ps in ("local_verified", "remote_downloadable")
            d["loci"] += 1
            d["carrier_usable"] += cu
            d["carrier_and_partner_usable"] += cu and pu
            d["carrier_usable_no_partner_usable"] += cu and not pu
            d["partner_only_usable"] += pu and not cu
            d["no_usable_reads"] += not cu and not pu
            d["carrier_" + cs] += 1
            d["best_partner_" + ps] += 1
            need_c |= runs_needed(by[g]["carrier"])
            need_p |= runs_needed(by[g]["partner"])
            local |= {r for m in by[g]["carrier"] + by[g]["partner"] for r in split(m["local_verified_runs"])}
        allr = need_c | need_p
        unknown = sorted(r for r in allr if not any(v.isdigit() for v in rinfo[r]["fastq_bytes"].split(";")))
        out[c] = dict(d, distinct_runs_to_download=len(allr), carrier_runs_to_download=len(need_c),
                      partner_runs_to_download=len(need_p - need_c), est_bytes_known=sum(nbytes(rinfo[r]["fastq_bytes"]) for r in allr),
                      runs_unknown_size=len(unknown), runs_unknown_size_ids=unknown[:50], local_verified_runs_used=len(local))
    return out


def traceability(core_rows, inst, prow, maps, cat, check, assign, fa):
    """ID uniqueness across batches and back-references from every mapping; each entry is (ok, detail)."""
    t = {}
    ids = [x["global_locus_id"] for x in core_rows]
    t["global_locus_id_unique"] = (len(ids) == len(set(ids)), len(ids))
    bios = defaultdict(set)
    asms = defaultdict(set)
    for sk, r in cat.items():
        bios[r["sample_id"]].add(sk.split(":")[0])
        asms[r["assembly"]].add(sk)
    multi_b = sorted(k for k, v in bios.items() if len(v) > 1)
    multi_a = sorted(k for k, v in asms.items() if len(v) > 1)
    t["biosample_in_one_batch"] = (not multi_b, multi_b[:20])
    t["assembly_in_one_sample_key"] = (not multi_a, multi_a[:20])
    t["core_rows_sample_key_in_catalog"] = (all(x["sample_key"] in cat and cat[x["sample_key"]]["assembly"] == x["assembly"] for x in core_rows), "")
    pc = {p["protein_sha256"] for p in prow}
    t["instances_in_protein_catalog"] = (all(i["protein_sha256"] in pc for i in inst), len(inst))
    t["locus_protein_hashes_in_catalog"] = (all(h in pc for x in core_rows for h in split(x["protein_sha256s"])), "")
    t["instance_locus_in_core_catalog"] = (set(i["global_locus_id"] for i in inst) <= set(ids), "")
    mapped = {m["local_id"] for m in maps}
    t["local_ids_namespaced_and_mapped"] = (all(":" in i["local_family_id"] and i["local_family_id"] in mapped and i["local_homology_id"] in mapped
                                                for i in inst), len(mapped))
    for lv in MMSEQS_SETTINGS:
        fails = [m for m in assign[lv] if not check[(lv, m)][4]]
        t["members_pass_final_check_" + lv] = (not fails, fails[:20])
        t["representatives_are_members_" + lv] = (all(check[(lv, m)][0] in assign[lv] for m in assign[lv]), "")
    n = sum(1 for line in open(str(fa)) if line.startswith(">"))
    t["core_fasta_one_record_per_locus"] = (n == len(core_rows), n)
    return dict(ok=all(v[0] for v in t.values()), checks={k: dict(ok=v[0], detail=v[1]) for k, v in t.items()})


def circular_contigs(gff):
    if not gff or not Path(gff).exists():
        return None
    out = subprocess.run(["grep", "-F", "Is_circular=true", gff], stdout=subprocess.PIPE, universal_newlines=True).stdout
    return {l.split("\t")[0] for l in out.splitlines() if "\tregion\t" in l}


def cmd_catalog(a):
    t0 = time.time()
    out, w = a.out, a.work
    out.mkdir(parents=True, exist_ok=True)
    bs = batches(a)
    # ---- global clusters and member checks
    assign, check, fixes = {}, {}, []
    for level in MMSEQS_SETTINGS:
        reps = {}
        for line in open(str(w / "mmseqs" / (level + "_cluster.tsv"))):
            rep, mem = line.split()
            reps[mem] = rep
        names = {r: "%s%s_%s" % (PREFIX[level], CATALOG, r[:12]) for r in set(reps.values())}
        if len(set(names.values())) != len(names):
            sys.exit("representative hash prefix collision at " + level)
        aln = {}
        for line in open(str(w / "mmseqs" / (level + "_members.m8"))):
            q, t, fid, qc, tc, _ = line.rstrip("\n").split("\t")
            if (t, q) not in aln or float(fid) > aln[(t, q)][0]:
                aln[(t, q)] = (float(fid), float(tc), float(qc))  # query = representative, target = member
        mid, mc = MEMBER_MIN[level]
        for m, r in reps.items():
            fid, mcov, rcov = (1.0, 1.0, 1.0) if m == r else aln.get((m, r), (0.0, 0.0, 0.0))
            check[(level, m)] = (r, fid, mcov, rcov, fid >= mid and mcov >= mc and rcov >= mc)
        for m, rec in reassign_failed(a, level, reps, check).items():
            fixes.append(rec)
            reps[m] = rec["new_rep"]
            check[(level, m)] = (rec["new_rep"], rec["new_identity"], rec["new_member_cov"], rec["new_rep_cov"], True)
            names.setdefault(rec["new_rep"], "%s%s_%s" % (PREFIX[level], CATALOG, rec["new_rep"][:12]))
        if len(set(names.values())) != len(names):
            sys.exit("representative hash prefix collision at " + level)
        assign[level] = {m: names[r] for m, r in reps.items()}
    write(out / "member_check_repairs.tsv", fixes, REPAIR_COLS)
    seqlen = {}
    name = None
    for line in open(str(w / "proteins.faa")):
        if line.startswith(">"):
            name = line[1:].strip()
        else:
            seqlen[name] = len(line.strip())
    inst = read(w / "protein_instances.tsv.gz")
    for i in inst:
        i["family_id"] = assign["family"][i["protein_sha256"]]
        i["homology_id"] = assign["homology_group"][i["protein_sha256"]]
    write(out / "protein_instances.tsv.gz", inst)
    by_locus = defaultdict(list)
    for i in inst:
        by_locus[i["global_locus_id"]].append(i)
    # ---- per-batch metadata
    cat, contigs, runs, links, files, organism, bioproject = {}, {}, {}, [], {}, {}, {}
    for b in bs:
        for r in read(HERE / POOL / b / "s2" / "catalog.tsv"):
            cat[gl(b, r["sample_id"])] = r
        for p in read(a.scratch / b / "s3_full" / "pairs.tsv"):
            for k in "ab":
                files[gl(b, p["sample_" + k])] = dict(fna=p["fna_" + k], gff=p["gff_" + k])
        for r in read(HERE / POOL / b / "runs.tsv"):
            r["batch"] = b
            runs[r["run_accession"]] = r
        for r in read(HERE / POOL / b / "assembly_run_links.tsv"):
            links.append(dict(batch=b, **r))
    # ---- core catalog
    lists = {}
    for b in bs:
        s4 = s4dir(a, b)
        for r in read(s4 / "candidate_loci_families.tsv") + read(s4 / ("list_%s.tsv" % PENDING)):
            lists[gl(b, r["locus_id"])] = (b, r)
    seeds = {}
    for b in bs:
        for r in read(s4dir(a, b) / "seed_catalog.tsv"):
            seeds[gl(b, r["locus_id"])] = r
    cores = {}
    for b in bs:
        name = None
        for line in open(str(s4dir(a, b) / "seed_catalog.fa")):
            if line.startswith(">"):
                name = gl(b, line[1:].split("|")[0])
            else:
                cores[name] = line.strip()
    by_carrier = defaultdict(list)
    for g, (b, r) in lists.items():
        by_carrier[gl(b, r["carrier"])].append(g)
    cidx = {}
    for b in bs:
        with gzip.open(str(a.scratch / b / "s2" / "contig_index.tsv.gz"), "rt") as h:
            next(h)
            for line in h:
                s, asm, c, n, csha = line.rstrip("\n").split("\t")[:5]
                cidx[(gl(b, s), c)] = (int(n), csha)
    core_rows, problems = [], Counter()
    fa = w / "candidate_cores.fa"
    with open(str(fa), "w") as fh:
        for sk in sorted(by_carrier):
            gff = files[sk]["gff"]
            circ = circular_contigs(gff)
            seqs = None
            for g in sorted(by_carrier[sk]):
                b, r = lists[g]
                s, e = int(r["core_start"]), int(r["core_end"])
                n, csha = cidx[(sk, r["contig"])]
                if g in seeds:
                    core, U, D = cores[g], seeds[g]["U"], seeds[g]["D"]
                    if sha(core) != seeds[g]["core_sha256"]:
                        sys.exit("seed core hash mismatch " + g)
                else:
                    if seqs is None:
                        seqs = read_fasta(files[sk]["fna"])
                    q = seqs[r["contig"]]
                    if len(q) != n:
                        sys.exit("contig length differs from the contig index " + g)
                    core, U, D = q[s:e], q[max(0, s - CONTEXT):s], q[e:e + CONTEXT]
                if len(core) != e - s:
                    problems["core_length_mismatch"] += 1
                c = cat[sk]
                pi = sorted(by_locus.get(g, []), key=lambda x: int(x["unit_order"]))
                fids, hids, strands = [x["family_id"] for x in pi], [x["homology_id"] for x in pi], [x["strand"] for x in pi]
                arch = canonical(fids, strands) if pi else ""
                harch = canonical(hids, strands) if pi else ""
                flags = []
                if r["tandem_like"] == "True":
                    flags.append("tandem_like")
                if int(r["core_copies_carrier"] or 0) > 1:
                    flags.append("core_multicopy_in_carrier")
                if int(r["core_copies_empty"] or 0) > 0:
                    flags.append("core_copy_in_empty")
                if r["boundary_consistency"] != "exact":
                    flags.append("boundary_" + r["boundary_consistency"])
                ends = []
                if len(U) < CONTEXT:
                    ends.append("U_short_contig_start")
                if len(D) < CONTEXT:
                    ends.append("D_short_contig_end")
                partners = split(r["partners_anchored"]) + split(r["partners_not_anchored"])
                fh.write(">%s|%s|%s:%d-%d|%s|%s\n%s\n" % (g, c["assembly"], r["contig"], s, e, r["category"], sha(core), core))
                core_rows.append(dict(
                    catalog_version=VERSION, evidence_source=EVIDENCE, read_evidence_status=READ_STATUS, global_locus_id=g, batch=b,
                    locus_id=r["locus_id"], category=r["category"], sample_key=sk, biosample=r["carrier"], assembly=c["assembly"],
                    paired_accession=c["genbank_accession"], bioproject=c["bioproject"], organism_label=r["organism"],
                    assembly_level=c["assembly_level"], annotation=c["annotation"], annotation_status=c["annotation_status"],
                    fna_sha256=c["fna_sha256"], contig=r["contig"], contig_length=n, contig_sha256=csha,
                    contig_circular="unknown_no_gff" if circ is None else "circular" if r["contig"] in circ else "not_flagged_circular",
                    core_start=s, core_end=e, core_len=e - s, core_strand="carrier_forward", core_sha256=sha(core),
                    context=CONTEXT, U_len=len(U), D_len=len(D), context_end_flags=";".join(ends) or "-",
                    U_start=max(0, s - CONTEXT), D_end=min(n, e + CONTEXT), U_sha256=sha(U), D_sha256=sha(D),
                    repeat_flags=";".join(flags) or "-", difference_class=r["difference_class"], anchored_support=r["anchored_support"],
                    n_pairs=r["n_pairs"], partners_anchored=r["partners_anchored"], partners_not_anchored=r["partners_not_anchored"],
                    partner_assemblies=";".join(cat.get(gl(b, p), {}).get("assembly", "not_in_catalog") for p in partners),
                    cds_class=r["cds_class"], cds_units=r["cds_units"], protein_ids=r["protein_ids"], products=r["products"],
                    strands=r["strands"], protein_sha256s=";".join(x["protein_sha256"] for x in pi), family_ids=";".join(fids),
                    family_architecture=arch, architecture_id=("A%s_%s" % (CATALOG, sha(arch)[:12])) if arch else "",
                    homology_ids=";".join(hids), homology_architecture=harch,
                    local_family_ids=";".join(x["local_family_id"] for x in pi), translation=r["translation"],
                    coding_flags=r["coding_flags"], partner_breakdown=r["partner_breakdown"], U=U, D=D))
    write(out / "candidate_core_catalog.tsv.gz", core_rows)
    print("core catalog", len(core_rows), "%.0fs" % (time.time() - t0), flush=True)
    # ---- protein catalog, families, architectures, local-global map
    pinst = defaultdict(list)
    for i in inst:
        pinst[i["protein_sha256"]].append(i)
    prow = []
    repaired = {(f["level"], f["member"]): f["action"] for f in fixes}
    for h, ii in sorted(pinst.items()):
        row = dict(protein_sha256=h, length=seqlen[h], instances=len(ii), loci=len({i["global_locus_id"] for i in ii}),
                   batches=";".join(sorted({i["batch"] for i in ii})), protein_ids=";".join(sorted({i["protein_id"] for i in ii})))
        for level, nm in (("family", "family"), ("homology_group", "homology")):
            r, fid, mcov, rcov, ok = check[(level, h)]
            row.update({nm + "_id": assign[level][h], nm + "_rep": r, nm + "_identity": round(fid, 3), nm + "_member_cov": round(mcov, 3),
                        nm + "_rep_cov": round(rcov, 3), nm + "_member_pass": ok,
                        nm + "_repair": repaired.get((level, h), "-")})
        prow.append(row)
    write(out / "protein_catalog.tsv.gz", prow)
    loc = {x["global_locus_id"]: x for x in core_rows}
    fam = defaultdict(lambda: defaultdict(set))
    for i in inst:
        d = fam[i["family_id"]]
        x = loc[i["global_locus_id"]]
        d["proteins"].add(i["protein_sha256"])
        d["accessions"].add(i["protein_id"])
        d["instances"].add((i["global_locus_id"], i["unit_order"]))
        d["loci"].add(i["global_locus_id"])
        d["samples"].add(i["sample_key"])
        d["biosamples"].add(i["biosample"])
        d["batches"].add(i["batch"])
        d["organisms"].add(x["organism_label"])
        d["architectures"].add(x["architecture_id"])
        d["homology"].add(i["homology_id"])
        d["local"].add(i["local_family_id"])
        d["products"].add(x["products"].split(";")[int(i["unit_order"]) - 1] if x["products"] else "")
        d["loci_" + i["category"]].add(i["global_locus_id"])
    ftab = []
    for f, d in fam.items():
        rep = check[("family", next(iter(d["proteins"])))][0]
        row = dict(family_id=f, catalog_version=VERSION, rep_sha256=rep, rep_length=seqlen[rep], homology_ids=";".join(sorted(d["homology"])),
                   distinct_proteins=len(d["proteins"]), protein_accessions=len(d["accessions"]), instances=len(d["instances"]),
                   loci=len(d["loci"]), samples=len(d["samples"]), biosamples=len(d["biosamples"]), batches=len(d["batches"]),
                   batch_names=";".join(sorted(d["batches"])), organism_labels=len(d["organisms"]), architectures=len(d["architectures"]),
                   members_failing_check=sum(1 for h in d["proteins"] if not check[("family", h)][4]),
                   local_family_ids=";".join(sorted(d["local"])), products=" | ".join(sorted(d["products"]))[:300])
        for c in LISTS:
            row["loci_" + c] = len(d["loci_" + c])
        ftab.append(row)
    ftab.sort(key=lambda x: (-x["loci"], x["family_id"]))
    write(out / "families.tsv", ftab)
    arch = defaultdict(lambda: defaultdict(set))
    for x in core_rows:
        if x["architecture_id"]:
            d = arch[x["architecture_id"]]
            d["arch"].add(x["family_architecture"])
            d["harch"].add(x["homology_architecture"])
            d["loci"].add(x["global_locus_id"])
            d["samples"].add(x["sample_key"])
            d["batches"].add(x["batch"])
            d["loci_" + x["category"]].add(x["global_locus_id"])
    atab = []
    for k, d in arch.items():
        row = dict(architecture_id=k, family_architecture=next(iter(d["arch"])), units=next(iter(d["arch"])).count(",") + 1,
                   homology_architectures=";".join(sorted(d["harch"])), loci=len(d["loci"]), samples=len(d["samples"]),
                   batches=len(d["batches"]), batch_names=";".join(sorted(d["batches"])))
        for c in LISTS:
            row["loci_" + c] = len(d["loci_" + c])
        atab.append(row)
    atab.sort(key=lambda x: (-x["loci"], x["architecture_id"]))
    write(out / "architectures.tsv", atab)
    maps = []
    for level, lk, gk in (("family", "local_family_id", "family_id"), ("homology_group", "local_homology_id", "homology_id")):
        l2g, g2l = defaultdict(set), defaultdict(set)
        lprot = defaultdict(set)
        for i in inst:
            l2g[i[lk]].add(i[gk])
            g2l[i[gk]].add(i[lk])
            lprot[i[lk]].add(i["protein_sha256"])
        for l, gs in sorted(l2g.items()):
            others = {x for g in gs for x in g2l[g]} - {l}
            rel = []
            if len(gs) > 1:
                rel.append("split")
            if any(x.split(":")[0] == l.split(":")[0] for x in others):
                rel.append("merged_within_batch")
            if any(x.split(":")[0] != l.split(":")[0] for x in others):
                rel.append("merged_across_batches")
            maps.append(dict(level=level, batch=l.split(":")[0], local_id=l, distinct_proteins=len(lprot[l]), global_ids=";".join(sorted(gs)),
                             other_local_ids=";".join(sorted(others)), relation=";".join(rel) or "one_to_one"))
    write(out / "local_global_family_map.tsv", maps)
    # ---- coding anomalies
    anom = []
    for b in bs:
        with gzip.open(str(s4dir(a, b) / "s4_loci.tsv.gz"), "rt") as h:
            for r in csv.DictReader(h, delimiter="\t"):
                if not r["translation"]:
                    continue
                us, ps, ts, fl = r["unit_ids"].split(";"), r["protein_ids"].split(";"), r["translation"].split(";"), r["coding_flags"].split(";")
                for u, p, t, f in zip(us, ps, ts, fl):
                    if t not in CODING_OK or "transl_table" in f:
                        anom.append(dict(batch=b, global_locus_id=gl(b, r["locus_id"]), sample_key=gl(b, r["carrier"]), contig=r["contig"],
                                         core_start=r["core_start"], core_end=r["core_end"], category=r["category"], cds_class=r["cds_class"],
                                         unit_id=u, protein_id=p, translation=t, coding_flags=f, partial_units=r["partial_units"],
                                         pseudo_units=r["pseudo_units"], notes=r["notes"][:300]))
    write(out / "coding_anomalies.tsv", anom)
    # ---- stage-2 manifests (metadata only)
    write(out / "assembly_run_links.tsv.gz", links)
    verified = {}
    for p in (BASELINE / "s5_readcheck_main").glob("verify_*.json"):
        v = json.loads(p.read_text())
        verified[v["run"]["run"]] = v
    rrows = []
    for run, r in sorted(runs.items()):
        v = verified.get(run)
        lf = v["local_file"] if v else ""
        rrows.append(dict(batch=r["batch"], run=run, biosample=r["sample_accession"], study=r["study_accession"],
                          platform=r["instrument_platform"], instrument_model=r["instrument_model"], read_type=r["read_type"],
                          read_type_evidence=r["read_type_evidence"], usable=r["usable"], unusable_reasons=r["unusable_reasons"],
                          base_count=r["base_count"], fastq_ftp=r["fastq_ftp"], fastq_bytes=r["fastq_bytes"], fastq_md5=r["fastq_md5"],
                          local_status=("verified_local" if lf and os.path.exists(lf) else "verified_record_file_missing") if v else "not_downloaded",
                          local_file=lf, local_sha256=v["sha256"] if v else "", verify_record=("comamonadaceae_20261024/s5_readcheck_main/verify_%s.json" % run) if v else "",
                          repair=json.dumps(v["repair"]) if v and v.get("repair") else ""))
    write(out / "read_run_manifest.tsv.gz", rrows)
    rinfo = {r["run"]: r for r in rrows}
    runs_by = defaultdict(set)
    for l in links:
        runs_by[gl(l["batch"], l["biosample"])].add(l["run"])
    for r in rrows:
        runs_by[gl(r["batch"], r["biosample"])].add(r["run"])
    mrows = []
    for x in core_rows:
        b = x["batch"]
        for role, ss in (("carrier", [x["biosample"]]), ("partner_anchored", split(x["partners_anchored"])),
                         ("partner_not_anchored", split(x["partners_not_anchored"]))):
            for s in ss:
                sk = gl(b, s)
                rs = sorted(runs_by.get(sk, ()))
                known = [r for r in rs if r in rinfo]
                use = [r for r in known if rinfo[r]["usable"] == "True"]
                loc_ok = [r for r in use if rinfo[r]["local_status"] == "verified_local"]
                remote = [r for r in use if r not in loc_ok]
                est = sum(nbytes(rinfo[r]["fastq_bytes"]) for r in remote)
                status = ("no_linked_runs" if not rs else "runs_without_metadata" if not known else "no_usable_runs" if not use
                          else "usable_runs_local_verified" if not remote else "usable_runs_partly_local" if loc_ok else "usable_runs_remote_only")
                mrows.append(dict(global_locus_id=x["global_locus_id"], batch=b, category=x["category"], role=role, sample_key=sk, biosample=s,
                                  assembly=cat[sk]["assembly"], linked_runs=";".join(rs), read_types=";".join(sorted({rinfo[r]["read_type"] for r in known})),
                                  usable_runs=";".join(use), local_verified_runs=";".join(loc_ok), runs_missing_metadata=";".join(r for r in rs if r not in rinfo),
                                  status=status, est_download_bytes=est,
                                  remote_runs_unknown_size=";".join(r for r in remote if not any(v.isdigit() for v in rinfo[r]["fastq_bytes"].split(";")))))
    write(out / "read_batch_manifest.tsv.gz", mrows)
    schema = dict(
        read_evidence=["evidence_id", "global_locus_id", "test_id", "run", "sample_key", "role", "caller_version", "count_unit",
                       "present", "absent", "present_matched", "absent_matched", "ambiguous", "locus_ambiguous", "off_locus",
                       "verdict", "evidence_level", "conflict", "created"],
        new_read_contexts=["context_id", "query_global_locus_id", "query_core_sha256", "run", "sample_key", "core_match_identity",
                           "core_match_cov", "U_context", "D_context", "U_sha256", "D_sha256", "molecules", "count_unit",
                           "assembly_placement", "placement_status", "review_status", "created"])
    for k, cols in schema.items():
        write(out / ("%s.tsv" % k), [], cols)
    # ---- pilot read evidence audit (baseline Comamonadaceae read check)
    plinks, ptests = pilot_audit(a, core_rows, cat)
    write(out / "pilot_read_evidence_links.tsv", plinks)
    write(out / "pilot_read_test_links.tsv", ptests)
    # ---- summary
    def stats(xs, ii):
        fams = {i["family_id"] for i in ii}
        return dict(loci=len(xs), samples=len({x["sample_key"] for x in xs}), assemblies=len({x["assembly"] for x in xs}),
                    bioprojects=len({x["bioproject"] for x in xs}), organism_labels=len({x["organism_label"] for x in xs}),
                    categories=dict(Counter(x["category"] for x in xs)), protein_instances=len(ii),
                    protein_accessions=len({i["protein_id"] for i in ii}), distinct_proteins=len({i["protein_sha256"] for i in ii}),
                    families=len(fams), homology_groups=len({i["homology_id"] for i in ii}),
                    architectures=len({x["architecture_id"] for x in xs if x["architecture_id"]}),
                    distinct_core_sequences=len({x["core_sha256"] for x in xs}))
    summ = dict(version=VERSION, evidence_source=EVIDENCE, read_evidence_status=READ_STATUS, batches=bs,
                member_min=MEMBER_MIN, settings=MMSEQS_SETTINGS, cluster_run=json.loads((w / "mmseqs" / "cluster_run.json").read_text()))
    summ["global"] = stats(core_rows, inst)
    summ["global"]["member_check"] = {lv: dict(members=len(assign[lv]), passed=sum(1 for m in assign[lv] if check[(lv, m)][4]),
                                               fail=sorted(m for m in assign[lv] if not check[(lv, m)][4])) for lv in MMSEQS_SETTINGS}
    pb = defaultdict(set)
    for i in inst:
        pb[i["protein_sha256"]].add(i["batch"])
    fb = {f["family_id"]: f["batches"] for f in ftab}
    summ["global"]["proteins_in_more_than_one_batch"] = sum(1 for v in pb.values() if len(v) > 1)
    summ["global"]["families_in_more_than_one_batch"] = sum(1 for v in fb.values() if v > 1)
    summ["global"]["sum_of_local_families"] = len({i["local_family_id"] for i in inst})
    summ["global"]["sum_of_local_distinct_proteins"] = len({(i["batch"], i["protein_sha256"]) for i in inst})
    summ["global"]["local_family_relations"] = dict(Counter(m["relation"] for m in maps if m["level"] == "family"))
    summ["per_batch"] = {}
    for b in bs:
        xs = [x for x in core_rows if x["batch"] == b]
        ii = [i for i in inst if i["batch"] == b]
        st = stats(xs, ii)
        st["local_families"] = len({i["local_family_id"] for i in ii})
        st["coding_anomaly_records"] = dict(Counter(r["translation"] for r in anom if r["batch"] == b))
        mm = [m for m in mrows if m["batch"] == b and m["role"] == "carrier"]
        st["carrier_read_status_loci"] = dict(Counter(m["status"] for m in mm))
        summ["per_batch"][b] = st
    summ["coding_anomalies"] = dict(records=len(anom), by_translation=dict(Counter(r["translation"] for r in anom)),
                                    loci=len({r["global_locus_id"] for r in anom}),
                                    by_category=dict(Counter(r["category"] for r in anom)))
    mm = [m for m in mrows if m["role"] == "carrier"]
    summ["read_metadata"] = dict(
        carrier_loci_status=dict(Counter(m["status"] for m in mm)),
        runs=len(rrows), runs_usable=sum(1 for r in rrows if r["usable"] == "True"),
        runs_verified_local=sum(1 for r in rrows if r["local_status"] == "verified_local"),
        distinct_usable_remote_runs_for_candidate_loci=len({r for m in mrows for r in split(m["usable_runs"]) if r not in split(m["local_verified_runs"])}),
        est_download_bytes_distinct_runs=sum(nbytes(rinfo[r]["fastq_bytes"]) for r in
                                             {r for m in mrows for r in split(m["usable_runs"]) if r not in split(m["local_verified_runs"])}),
        manifest_rows=len(mrows), downloads_started=0)
    summ["pilot_read_evidence_links"] = dict(loci=len(plinks), link_basis=dict(Counter(p["link_basis"] for p in plinks)),
                                             unmatched_reason=dict(Counter(p["unmatched_reason"] for p in plinks if p["unmatched_reason"] != "-")),
                                             tests=len(ptests), test_link_scope=dict(Counter(t["link_scope"] for t in ptests)))
    summ["member_check_repairs"] = fixes
    summ["stage2_read_coverage"] = read_coverage(core_rows, mrows, rinfo)
    summ["traceability"] = traceability(core_rows, inst, prow, maps, cat, check, assign, fa)
    summ["problems"] = dict(problems)
    (out / "catalog_summary.json").write_text(json.dumps(summ, indent=1) + "\n")
    rebuild = "git checkout %s; sbatch experiments/fna_pool/stage1_catalog.sbatch (proteins -> cluster -> catalog), or the single mode named" % a.commit
    big = [(w / "proteins.faa", "proteins"), (w / "protein_instances.tsv.gz", "proteins"), (fa, "catalog")] + \
          [(w / "mmseqs" / (lv + x), "cluster") for lv in MMSEQS_SETTINGS for x in ("_cluster.tsv", "_members.m8")] + \
          [(w / "mmseqs" / ("repair_" + lv) / "hits.m8", "catalog") for lv in MMSEQS_SETTINGS if (w / "mmseqs" / ("repair_" + lv) / "hits.m8").exists()]
    inputs = {}
    for b in bs:
        f = {"s4/" + n: s4dir(a, b) / n for n in ("candidate_loci_families.tsv", "seed_catalog.fa", "seed_catalog.tsv", "list_annotation_pending.tsv",
                                                    "protein_members.tsv", "families_summary.json", "s4_summary.json", "s4_loci.tsv.gz")}
        f.update({"s3_full/pairs.tsv": a.scratch / b / "s3_full" / "pairs.tsv", "s2/contig_index.tsv.gz": a.scratch / b / "s2" / "contig_index.tsv.gz",
                  "repo/s2/catalog.tsv": HERE / POOL / b / "s2" / "catalog.tsv", "repo/runs.tsv": HERE / POOL / b / "runs.tsv",
                  "repo/assembly_run_links.tsv": HERE / POOL / b / "assembly_run_links.tsv"})
        if b == "comamonadaceae":
            f["s3_full/merged/loci.tsv"] = a.scratch / b / "s3_full" / "merged" / "loci.tsv"
        inputs[b] = {k: dict(path=str(v), sha256=sha256(v)) for k, v in f.items()}
    inputs["baseline_pilot"] = {str(p.relative_to(HERE)): sha256(p) for p in
                                [BASELINE / "s5_readcheck_main" / "readcheck_loci.tsv", BASELINE / "s5_readcheck_main" / "readcheck_tests.tsv",
                                 BASELINE / "s4_s3_full_v02" / "seed_catalog.tsv", BASELINE / "s2" / "catalog.tsv"] +
                                [BASELINE / "s4_s3_full_v02" / ("list_%s.tsv" % c) for c in LISTS]}
    cr = summ["cluster_run"]
    idx = dict(version=VERSION, created=time.strftime("%Y-%m-%dT%H:%M:%S%z"), commit=a.commit, evidence_source=EVIDENCE,
               parameters=dict(context=CONTEXT, lists=LISTS, pending=PENDING, mmseqs_settings=MMSEQS_SETTINGS, member_min=MEMBER_MIN,
                               member_align=ALIGN, id_scheme="F%s_/H%s_/A%s_ + sha256[:12]; <batch>:<id> for loci, samples, local families" % (CATALOG, CATALOG, CATALOG)),
               software=dict(python=sys.version.split()[0], mmseqs=cr["mmseqs"], mmseqs_version=cr["mmseqs_version"], mmseqs_sha256=sha256(cr["mmseqs"])),
               scratch_files={str(p): dict(sha256=sha256(p), bytes=p.stat().st_size, produced_by_mode=m, rebuild=rebuild) for p, m in big},
               inputs=inputs)
    (out / "frozen_manifest.json").write_text(json.dumps(idx, indent=1) + "\n")
    if (out / "scratch_files.json").exists():
        (out / "scratch_files.json").unlink()
    print(json.dumps(summ["global"], indent=1)[:4000])
    print("done %.0fs" % (time.time() - t0))


def cmd_report(a):
    """Per-batch and global stage-1 counts from S1 / S2 / S3 / S4 / catalog outputs (no recomputation)."""
    out = a.out
    tb = {r["family"].lower(): r for r in read(HERE / POOL / "target_batches.tsv")}
    cs = json.loads((out / "catalog_summary.json").read_text())
    overlap = read(HERE / POOL / "cross_batch_overlap.tsv")
    rep = dict(version=VERSION, evidence_source=EVIDENCE, read_evidence_status=READ_STATUS, per_batch={},
               cross_batch_overlap_rows=len([r for r in overlap if r.get("id")]))
    keys = ["pool_assemblies", "pool_biosamples", "pool_projects", "pool_annotated", "pool_annotation_coverage",
            "pool_samples_with_usable_long_reads", "reads_meta_eligible_wgs_fastq_samples"]
    for b in cs["batches"]:
        s1 = json.loads((HERE / POOL / b / "s1" / "s1_report.json").read_text())
        s2 = json.loads((HERE / POOL / b / "s2" / "s2_summary.json").read_text())
        s3 = json.loads((HERE / POOL / b / "s3_full" / "summary.json").read_text())
        s4 = json.loads((s4dir(a, b) / "s4_summary.json").read_text())
        fs = json.loads((s4dir(a, b) / "families_summary.json").read_text())
        excl = HERE / POOL / b / "s2" / "catalog_excluded.tsv"
        cb = cs["per_batch"][b]
        rep["per_batch"][b] = dict(
            batch_id=tb[b]["batch_id"], taxon_id=tb[b]["taxon_id"], **{k: tb[b][k] for k in keys},
            downloads=dict(rows=s1["rows"], status=s1["status"], failed=s1["failed"], annotation=s1["annotation"],
                           gca_gcf_pair_checks=s1["pair_checks"]),
            s2=dict(samples=s2["samples"], excluded_not_verified=len(read(excl)) if excl.exists() else 0, pair_status=s2["pair_status"],
                    pairs=s2["pairs"], species_clusters=s2["species_clusters"], clone_groups=s2["clone_groups"]),
            s3=dict(pairs=s3["pairs"], raw_calls=s3["calls"], loci=s3["loci"], status=s3["status"], anchored_lt2500=s3["anchored_lt2500"]),
            s4=dict(short_anchored_loci=s4["loci"], carriers=s4["carriers"], categories=s4["categories"], translation=s4["translation"],
                    local_distinct_proteins=fs["distinct_proteins"], local_families=fs["all"]["families"],
                    local_architectures=fs["all"]["family_architectures"]),
            catalog=cb)
    g = {}
    for k in ("pool_assemblies", "pool_biosamples", "pool_projects", "pool_annotated", "pool_samples_with_usable_long_reads"):
        g[k + "_sum"] = sum(int(v[k]) for v in rep["per_batch"].values())
    pb = rep["per_batch"].values()
    g["downloads"] = dict(Counter({k: c for v in pb for k, c in v["downloads"]["status"].items()}))
    g["s2_samples"] = sum(v["s2"]["samples"] for v in pb)
    g["s2_pairs"] = sum(v["s2"]["pairs"] for v in pb)
    g["s2_no_pair_selected"] = sum(v["s2"]["pair_status"].get("no_pair_selected_pending_review", 0) for v in pb)
    g["s3_raw_calls"] = sum(v["s3"]["raw_calls"] for v in pb)
    g["s3_loci"] = sum(v["s3"]["loci"] for v in pb)
    g["s3_anchored_lt2500"] = sum(v["s3"]["anchored_lt2500"] for v in pb)
    g["s4_short_anchored_loci"] = sum(v["s4"]["short_anchored_loci"] for v in pb)
    g["s4_categories"] = dict(Counter({k: c for v in pb for k, c in v["s4"]["categories"].items()}))
    g["s4_categories_close"] = sum(g["s4_categories"].values()) == g["s4_short_anchored_loci"]
    g["s4_translation"] = dict(Counter({k: c for v in pb for k, c in v["s4"]["translation"].items()}))
    g["sum_local_families"] = sum(v["s4"]["local_families"] for v in pb)
    g["sum_local_distinct_proteins"] = sum(v["s4"]["local_distinct_proteins"] for v in pb)
    g["catalog"] = cs["global"]
    g["coding_anomalies"] = cs["coding_anomalies"]
    g["read_metadata"] = cs["read_metadata"]
    g["pilot_read_evidence_links"] = cs["pilot_read_evidence_links"]
    g["member_check_repairs"] = cs.get("member_check_repairs")
    g["stage2_read_coverage"] = cs.get("stage2_read_coverage")
    g["traceability"] = cs.get("traceability")
    g["site_groups"] = "not computed in stage 1 (S6 merging runs after read evidence); counts here are loci, not site groups"
    rep["global"] = g
    (out / "stage1_report.json").write_text(json.dumps(rep, indent=1) + "\n")
    with open(str(out / "outputs.sha256"), "w") as h:
        for p in sorted(out.iterdir()):
            if p.name != "outputs.sha256":
                h.write("%s  %s\n" % (sha256(p), p.name))
    print(json.dumps(g, indent=1))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("mode", choices=["proteins", "cluster", "catalog", "report"])
    p.add_argument("--scratch", type=Path, default=PROJ / "fna_pool" / POOL)
    p.add_argument("--work", type=Path, default=PROJ / "fna_pool" / POOL / "stage1_catalog")
    p.add_argument("--out", type=Path, default=HERE / POOL / "stage1_catalog")
    p.add_argument("--mmseqs", default=str(PROJ / "env/s4/bin/mmseqs"))
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--commit", default="")
    p.add_argument("--batches", nargs="*", help="restrict to these batches (tests)")
    a = p.parse_args()
    dict(proteins=cmd_proteins, cluster=cmd_cluster, catalog=cmd_catalog, report=cmd_report)[a.mode](a)


if __name__ == "__main__":
    main()
