"""S4: GFF3 annotation of short anchored cores, candidate lists, preliminary protein families, DNA seed catalog.

Input: S3 loci/calls (merge output) and the pair table with FNA/GFF3 paths (protein.faa next to genomic.gff).
Coordinates are 0-based half-open internally; GFF3 is 1-based closed. Nothing here changes S3 coordinates: gene ends
only flag a locus for review, the DNA core boundary is never moved to a CDS boundary.

review    One row per anchored locus with core < SHORT bp (s4_loci.tsv.gz) plus annotation_files.tsv (GFF3/FAA sha256
          and annotation name per carrier; GFF IDs are only meaningful with that hash).
  Units   GFF3 CDS/RNA lines grouped by ID (multi-segment CDS such as programmed frameshifts are one unit). Each unit
          keeps strand, segments with phase, Parent, partial (partial=true / start_range / end_range), pseudo,
          exception, transl_except, protein_id, product, Note. A unit whose end exceeds the contig length crosses the
          origin of a circular contig and is handled as its two pieces.
  Coding  Every non-pseudo CDS unit inside a core is translated from the FNA (table 11; segments in transcript order,
          phase of the 5' segment; an alternative start codon reads as M) and compared with its protein in the FAA;
          start codon, terminal stop and internal stops are checked. "complete" = inside the core, not partial or
          pseudo, translation equals the FAA protein, valid start and stop, no internal stop. transl_except entries
          (e.g. selenocysteine) replace exactly the codon they name, located through the unit's segments; any other
          internal stop still fails. Exceptions such as ribosomal slippage are kept as flags.
  Bounds  Features crossing a core boundary: host_gene_interrupted (one feature spans the whole core),
          gene_mostly_inside (the core may be truncated: >= half of the feature is in the core), gene_mostly_outside
          (the core end overlaps a flank gene), rna_crossed (rRNA/tRNA/other RNA).
  Lists   category, in this order:
            main                 presence_absence, all_calls_anchored, 1-2 complete CDS, nothing crossing either
                                 boundary, boundary_consistency exact, not tandem-like
            partner_mixed        the same but some partner's calls are replacement/unresolved (partner_breakdown
                                 says which partner gives what)
            direction_mixed      the same but one mapping direction is not anchored
            short_other_segment  complex_small_other (1-49 bp of other sequence at the partner site), otherwise
                                 the main criteria (support kept as a column)
            replacement          replacement (>= 50 bp of other sequence), otherwise the main criteria
            boundary_intersection   a feature crosses a core boundary
            gene_integrity_anomaly  an inside CDS is pseudo/partial/translation or codon problem
            boundary_or_repeat_review  annotation passes but boundaries are not exact or the core is tandem-like
            other_structure      0 CDS or > 2 complete CDS
families  Proteins of complete CDS in main / partner_mixed / direction_mixed / short_other_segment / replacement,
          deduplicated by sequence, clustered with MMseqs2 easy-cluster at two frozen settings: family (min-seq-id
          0.9, coverage 0.8 of both sequences, cov-mode 0, cluster-mode 0 = set cover, every member aligned to the
          representative, --cluster-reassign 1) and homology_group (min-seq-id 0.5, same coverage) as candidate
          distant homology only. Cascaded clustering can leave a member below the threshold against its final
          representative, so every member is aligned to its representative (exhaustive all-vs-all MMseqs2 search,
          alignment-mode 3) and its identity, both coverages and pass/fail are kept (protein_members.tsv).
          Architecture = family ids in core order with orientation relative to the core, canonicalised over the two
          strands. Family table reports proteins, loci, carrier samples, species and architectures separately, per
          list. These are preliminary families, not final S6 families.
seeds     DNA seed catalog: the observed core DNA (carrier orientation) of every locus in the candidate lists, with
          U/D context, one record per locus (members, not consensus), sha256 per record.
"""

import argparse
import csv
import gzip
import hashlib
import json
import re
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path
from urllib.parse import unquote

sys.path.insert(0, str(Path(__file__).resolve().parent))
from s3_pair_calls import rc, read_fasta  # noqa: E402

VERSION = "s4_annotate_v0.2"
SHORT = 2500
CONTEXT = 300
LISTS = ["main", "partner_mixed", "direction_mixed", "short_other_segment", "replacement"]
MMSEQS_SETTINGS = {"family": ["--min-seq-id", "0.9", "-c", "0.8", "--cov-mode", "0", "--cluster-mode", "0", "--cluster-reassign", "1"],
                   "homology_group": ["--min-seq-id", "0.5", "-c", "0.8", "--cov-mode", "0", "--cluster-mode", "0", "--cluster-reassign", "1"]}
MEMBER_MIN = {"family": (0.9, 0.8), "homology_group": (0.5, 0.8)}  # identity, coverage of both sequences
MMSEQS_CHECK = ["--exhaustive-search", "1", "--alignment-mode", "3", "-e", "10", "--max-seqs", "100000",
                "--format-output", "query,target,fident,qcov,tcov,alnlen"]
BASES = "TCAG"
# NCBI genetic codes (amino acids, starts) in TCAG order; a CDS uses the transl_table of its GFF3 line. A CDS
# without transl_table, with conflicting values or with a table not listed here is reviewed, not translated.
GENETIC_CODES = {
    11: ("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
         "---M------**--*----M------------MMMM---------------M------------"),
    4: ("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "--MM------**-------M------------MMMM---------------M------------"),
}
TRIPLETS = [a + b + c for a in BASES for b in BASES for c in BASES]
CODES = {t: dict(codon=dict(zip(TRIPLETS, aa)), starts={x for x, m in zip(TRIPLETS, st) if m == "M"},
                 stops={x for x, m in zip(TRIPLETS, aa) if m == "*"}) for t, (aa, st) in GENETIC_CODES.items()}
CODON, STARTS, STOPS = CODES[11]["codon"], CODES[11]["starts"], CODES[11]["stops"]
RNA = {"rRNA", "tRNA", "tmRNA", "ncRNA", "RNase_P_RNA", "SRP_RNA", "antisense_RNA"}


def read(path):
    op = gzip.open if str(path).endswith(".gz") else open
    with op(str(path), "rt") as h:
        return list(csv.DictReader(h, delimiter="\t"))


def write(path, rows, cols=None):
    cols = cols or (list(rows[0]) if rows else ["locus_id"])
    op = gzip.open if str(path).endswith(".gz") else open
    with op(str(path), "wt") as h:
        w = csv.DictWriter(h, cols, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def sha256(path):
    h = hashlib.sha256()
    with open(str(path), "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def read_faa(path):
    out, name, buf = {}, None, []
    for line in open(str(path)):
        if line.startswith(">"):
            if name:
                out[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name:
        out[name] = "".join(buf)
    return out


def parse_gff(path, lengths):
    """{contig: [unit]}; unit segments are (start, end, phase) 0-based half-open, end may exceed the contig length
    for an origin-crossing feature (a segment of a multi-segment CDS can lie wholly past the end)."""
    units = defaultdict(dict)
    for line in open(str(path)):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or (f[2] != "CDS" and f[2] not in RNA):
            continue
        at = {k: unquote(v) for k, v in (kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)}
        key = at.get("ID") or "%s:%s:%s:%s" % (f[2], f[0], f[3], f[4])
        u = units[f[0]].get(key)
        if u is None:
            u = units[f[0]][key] = dict(id=key, type=f[2], strand=f[6], segments=[], parent=at.get("Parent", ""), partial=False,
                                        pseudo=False, exception=at.get("exception", ""), transl_except=at.get("transl_except", ""),
                                        protein=at.get("protein_id", ""), product=at.get("product", ""), note=at.get("Note", ""),
                                        transl_tables=set())
        u["segments"].append((int(f[3]) - 1, int(f[4]), int(f[7]) if f[7] in "012" else 0))
        u["transl_tables"].add(at.get("transl_table", ""))
        u["partial"] |= at.get("partial") == "true" or "start_range" in at or "end_range" in at
        u["pseudo"] |= at.get("pseudo") == "true"
    out = {}
    for contig, us in units.items():
        n = lengths.get(contig, 0)
        for u in us.values():
            u["segments"].sort()
            u["crosses_origin"] = bool(n) and any(e > n for _, e, _ in u["segments"])
            iv = []
            for s, e, _ in u["segments"]:
                if n and s >= n:
                    iv.append((s - n, e - n))
                elif n and e > n:
                    iv += [(s, n)] + ([(0, e - n)] if e - n > 0 else [])
                else:
                    iv.append((s, e))
            # one interval per unit unless it crosses the origin (then its pieces)
            u["intervals"] = iv if u["crosses_origin"] else [(min(s for s, _ in iv), max(e for _, e in iv))]
        out[contig] = list(us.values())
    return out


def unit_nt(u, seq):
    """CDS nucleotides in transcript order (segments joined; minus strand reverse-complemented) and 5' phase."""
    n = len(seq)
    parts = []
    for s, e, _ in u["segments"]:
        parts.append(seq[s - n:e - n] if s >= n else seq[s:n] + seq[0:e - n] if e > n else seq[s:e])
    if u["strand"] == "-":
        return "".join(rc(p) for p in reversed(parts)), u["segments"][-1][2]
    return "".join(parts), u["segments"][0][2]


AA3 = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
       "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
       "Sec": "U", "Pyl": "O", "TERM": "*", "OTHER": "X"}
TRANSL_EXCEPT = re.compile(r"\(pos:(complement\()?(\d+)\.\.(\d+)\)?,aa:(\w+)\)")


def transl_exceptions(u, n):
    """{codon index: amino acid} from transl_except: the codon whose 5' base is at the stated position (the high end
    for complement), located through the unit's own segments in transcript order; unplaceable entries map to None."""
    pos = []
    for s, e, _ in u["segments"]:
        pos += [x % n for x in range(s, e)]
    if u["strand"] == "-":
        pos.reverse()
    phase = u["segments"][-1][2] if u["strand"] == "-" else u["segments"][0][2]
    out = {}
    for comp, a, b, aa in TRANSL_EXCEPT.findall(u["transl_except"]):
        five = int(b) - 1 if comp else int(a) - 1
        k = pos.index(five) - phase if five in pos else -1
        out[k // 3 if k >= 0 and k % 3 == 0 else None] = AA3.get(aa, AA3.get(aa.capitalize(), "X"))
    return out


def translate(nt, table=11):
    codon = CODES[table]["codon"]
    return "".join(codon.get(nt[i:i + 3], "X") for i in range(0, len(nt) - 2, 3))


def unit_table(u):
    """(table, flag): the unit's single supported transl_table, else (None, reason)."""
    tt = u.get("transl_tables", {"11"})
    if tt == {""}:
        return None, "transl_table_missing"
    if len(tt) > 1:
        return None, "transl_table_conflict=" + "/".join(sorted(x or "none" for x in tt))
    t = next(iter(tt))
    if not t.isdigit() or int(t) not in CODES:
        return None, "unsupported_transl_table=" + t
    return int(t), ""


def coding_check(u, seq, faa):
    """Flags for one CDS unit: translation vs FAA, start, stop, internal stops, under the unit's transl_table."""
    table, why = unit_table(u)
    if table is None:
        return ("pseudo" if u["pseudo"] else "transl_table_review"), [why]
    code = CODES[table]
    nt, phase = unit_nt(u, seq)
    nt = nt[phase:].upper()
    aa = translate(nt, table)
    flags = []
    if not u["partial"]:
        if nt[:3] not in code["starts"]:
            flags.append("no_start")
        if nt[-3:] not in code["stops"]:
            flags.append("no_stop")
        if len(nt) % 3:
            flags.append("length_not_codon_multiple")
        if aa[:1] and nt[:3] in code["starts"]:
            aa = "M" + aa[1:]
    exc = transl_exceptions(u, len(seq)) if u["transl_except"] else {}
    if None in exc:
        flags.append("transl_except_unplaced")
    raw = aa
    for k, x in exc.items():
        if k is not None and k < len(aa):
            aa = aa[:k] + x + aa[k + 1:]
    body = aa[:-1] if aa.endswith("*") else aa
    if "*" in body:
        flags.append("internal_stop")
    if u["pseudo"]:
        return "pseudo", flags
    if not u["protein"]:
        return "no_protein", flags + ["no_protein"]
    p = faa.get(u["protein"])
    if p is None:
        return "protein_not_in_faa", flags + ["protein_not_in_faa"]
    if body == p:
        return ("match_with_transl_except" if aa != raw else "match"), flags
    return "mismatch", flags + ["translation_mismatch"]


def boundary_relations(us, s, e):
    rel = []
    for u in us:
        for a, b in u["intervals"]:
            if not (a < s < b or a < e < b):
                continue
            if u["type"] != "CDS":
                rel.append(("rna_crossed", u))
            elif a < s and b > e:
                rel.append(("host_gene_interrupted", u))
            else:
                inside = max(0, min(b, e) - max(a, s))
                rel.append(("gene_mostly_inside" if inside >= (b - a) / 2.0 else "gene_mostly_outside", u))
            break
    return rel


def categorise(l, cds_class, rel):
    structure_ok = l["boundary_consistency"] == "exact" and l["tandem_like"] == "False"
    if rel:
        return "boundary_intersection"
    if cds_class in ("contains_pseudo", "contains_partial", "integrity_anomaly"):
        return "gene_integrity_anomaly"
    if cds_class != "1-2_complete_CDS":
        return "other_structure"
    if not structure_ok:
        return "boundary_or_repeat_review"
    if l["difference_class"] == "replacement":
        return "replacement"
    if l["difference_class"] == "complex_small_other":
        return "short_other_segment"
    return {"all_calls_anchored": "main", "partner_mixed": "partner_mixed"}.get(l["anchored_support"], "direction_mixed")


def partner_breakdown(calls):
    by = defaultdict(lambda: dict(calls=0, status=Counter(), difference_class=Counter(), other_side_gap=[]))
    for c in calls:
        b = by[c["empty"]]
        b["calls"] += 1
        b["status"][c["status"]] += 1
        b["difference_class"][c.get("difference_class", "")] += 1
        b["other_side_gap"].append(int(c["other_side_gap"]))
    return json.dumps({p: dict(calls=b["calls"], status=dict(b["status"]), difference_class=dict(b["difference_class"]),
                               other_side_gap=sorted(set(b["other_side_gap"]))) for p, b in sorted(by.items())}, sort_keys=True)


def sample_files(pairs):
    out = {}
    for p in read(pairs):
        for k in "ab":
            out[p["sample_" + k]] = dict(fna=p["fna_" + k], gff=p["gff_" + k], organism=p["organism_" + k])
    return out


def cmd_review(a):
    files = sample_files(a.pairs)
    loci = [l for l in read(a.loci) if l["status"] == "two_flank_anchored" and int(l["core_len"]) < SHORT]
    want = {l["locus_id"] for l in loci}
    calls = defaultdict(list)
    for c in read(a.calls):
        if c["locus_id"] in want:
            calls[c["locus_id"]].append(c)
    by_carrier = defaultdict(list)
    for l in loci:
        by_carrier[l["carrier"]].append(l)
    rows, ann = [], []
    for carrier, ls in sorted(by_carrier.items()):
        f = files[carrier]
        seqs = read_fasta(f["fna"])
        gff = Path(f["gff"])
        faa_path = gff.parent / "protein.faa"
        faa = read_faa(faa_path)
        units = parse_gff(gff, {k: len(v) for k, v in seqs.items()})
        v = json.load(open(str(gff.parent / "verify.json")))
        ann.append(dict(sample=carrier, organism=f["organism"], accession=v.get("accession", ""), annotation=v.get("annotation", ""),
                        gff=str(gff), gff_sha256=sha256(gff), faa_sha256=sha256(faa_path), loci=len(ls)))
        for l in ls:
            s, e = int(l["core_start"]), int(l["core_end"])
            seq, us = seqs[l["contig"]], units.get(l["contig"], [])
            rel = boundary_relations(us, s, e)
            inside = [u for u in us if u["type"] == "CDS" and not u["crosses_origin"] and s <= u["intervals"][0][0] and u["intervals"][0][1] <= e]
            inside.sort(key=lambda u: u["intervals"][0])
            checks = [coding_check(u, seq, faa) for u in inside]
            ok = [(not u["pseudo"] and not u["partial"] and st in ("match", "match_with_transl_except") and not fl)
                  for u, (st, fl) in zip(inside, checks)]
            if any(u["pseudo"] for u in inside):
                cds_class = "contains_pseudo"
            elif any(u["partial"] for u in inside):
                cds_class = "contains_partial"
            elif not all(ok):
                cds_class = "integrity_anomaly"
            else:
                cds_class = "0_CDS" if not inside else "1-2_complete_CDS" if len(inside) <= 2 else ">2_complete_CDS"
            orient = ["+" if u["strand"] == "+" else "-" for u in inside]
            rows.append(dict(
                locus_id=l["locus_id"], carrier=carrier, organism=f["organism"], contig=l["contig"], core_start=s, core_end=e,
                core_len=e - s, difference_class=l["difference_class"], anchored_support=l["anchored_support"],
                partners_anchored=l["partners_anchored"], partners_not_anchored=l["partners_not_anchored"],
                boundary_consistency=l["boundary_consistency"], tandem_like=l["tandem_like"], n_pairs=l["n_pairs"],
                core_copies_carrier=l["core_copies_carrier"], core_copies_empty=l["core_copies_empty"],
                cds_class=cds_class, category=categorise(l, cds_class, rel),
                boundary_relations=";".join(sorted({r for r, _ in rel})),
                crossing_features=";".join("%s:%s:%s" % (r, u["type"], u["product"]) for r, u in rel),
                cds_units=len(inside), unit_ids=";".join(u["id"] for u in inside), protein_ids=";".join(u["protein"] for u in inside),
                products=";".join(u["product"] for u in inside), strands=";".join(orient),
                segments=";".join(str(len(u["segments"])) for u in inside), parents=";".join(u["parent"] for u in inside),
                translation=";".join(st for st, _ in checks), coding_flags=";".join(",".join(fl) or "-" for _, fl in checks),
                exceptions=";".join(u["exception"] or "-" for u in inside), partial_units=sum(u["partial"] for u in inside),
                pseudo_units=sum(u["pseudo"] for u in inside), notes=";".join(u["note"] or "-" for u in inside),
                partner_breakdown=partner_breakdown(calls[l["locus_id"]])))
        print(carrier, len(ls), file=sys.stderr)
    a.out.mkdir(parents=True, exist_ok=True)
    write(a.out / "s4_loci.tsv.gz", rows)
    write(a.out / "annotation_files.tsv", ann)
    for cat in LISTS:
        write(a.out / ("list_%s.tsv" % cat), [r for r in rows if r["category"] == cat], list(rows[0]))
    summ = dict(version=VERSION, short=SHORT, loci=len(rows), carriers=len(ann), categories=dict(Counter(r["category"] for r in rows)),
                cds_class=dict(Counter(r["cds_class"] for r in rows)),
                boundary_relations=dict(Counter(x for r in rows for x in r["boundary_relations"].split(";") if x)),
                translation=dict(Counter(x for r in rows for x in r["translation"].split(";") if x)),
                coding_flags=dict(Counter(f for r in rows for x in r["coding_flags"].split(";") for f in x.split(",") if f not in ("", "-"))),
                lists={c: dict(loci=sum(1 for r in rows if r["category"] == c),
                               carriers=len({r["carrier"] for r in rows if r["category"] == c}),
                               cds_units=sum(int(r["cds_units"]) for r in rows if r["category"] == c),
                               protein_ids=len({p for r in rows if r["category"] == c for p in r["protein_ids"].split(";") if p}))
                       for c in LISTS})
    (a.out / "s4_summary.json").write_text(json.dumps(summ, indent=1) + "\n")
    print(json.dumps(summ, indent=1))


def cmd_families(a):
    rows = [r for c in LISTS for r in read(a.s4 / ("list_%s.tsv" % c))]
    files = sample_files(a.pairs)
    faas, seqs, members, loc_prots = {}, {}, defaultdict(list), {}  # seqs/members keyed by protein sequence hash
    for r in rows:
        if r["carrier"] not in faas:
            faas[r["carrier"]] = read_faa(Path(files[r["carrier"]]["gff"]).parent / "protein.faa")
        hs = []
        for pid in r["protein_ids"].split(";"):
            s = faas[r["carrier"]][pid]
            h = hashlib.sha256(s.encode()).hexdigest()[:16]
            seqs[h] = s
            members[h].append((r["locus_id"], r["carrier"], pid))
            hs.append(h)
        loc_prots[r["locus_id"]] = hs
    a.work.mkdir(parents=True, exist_ok=True)
    fa = a.work / "proteins.faa"
    with open(str(fa), "w") as h:
        for k, s in sorted(seqs.items()):
            h.write(">%s\n%s\n" % (k, s))
    assign, reps_by = {}, {}
    for level, opts in MMSEQS_SETTINGS.items():
        prefix = a.work / level
        subprocess.run([a.mmseqs, "easy-cluster", str(fa), str(prefix), str(a.work / ("tmp_" + level)), "--threads", str(a.threads)] + opts,
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        reps = {}
        for line in open(str(prefix) + "_cluster.tsv"):
            rep, mem = line.split()
            reps[mem] = rep
        names = {rep: "%s_%04d" % ("F" if level == "family" else "H", i + 1) for i, rep in enumerate(sorted(set(reps.values()),
                 key=lambda x: (-sum(1 for v in reps.values() if v == x), x)))}
        assign[level] = {m: names[r] for m, r in reps.items()}
        reps_by[level] = reps
    aln = {}
    res = a.work / "allvsall.m8"
    subprocess.run([a.mmseqs, "easy-search", str(fa), str(fa), str(res), str(a.work / "tmp_check"), "--threads", str(a.threads)]
                   + MMSEQS_CHECK, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    for line in open(str(res)):
        q, t, fid, qc, tc, _ = line.split("\t")
        if (q, t) not in aln or float(fid) > aln[(q, t)][0]:
            aln[(q, t)] = (float(fid), float(qc), float(tc))
    check = {}
    for level, (mid, mc) in MEMBER_MIN.items():
        for m, rep_ in reps_by[level].items():
            fid, qc, tc = (1.0, 1.0, 1.0) if m == rep_ else aln.get((m, rep_), (0.0, 0.0, 0.0))
            check[(level, m)] = (rep_, fid, qc, tc, fid >= mid and qc >= mc and tc >= mc)
    ver = subprocess.run([a.mmseqs, "version"], stdout=subprocess.PIPE, universal_newlines=True).stdout.strip()
    species = {r["carrier"]: r["organism"] for r in rows}
    out_loci, fam = [], defaultdict(lambda: defaultdict(set))
    for r in rows:
        hs = loc_prots[r["locus_id"]]
        strands = r["strands"].split(";")
        for level in MMSEQS_SETTINGS:
            ids = [assign[level][h] for h in hs]
            fwd = ",".join("%s%s" % (i, s) for i, s in zip(ids, strands))
            rev = ",".join("%s%s" % (i, "+" if s == "-" else "-") for i, s in reversed(list(zip(ids, strands))))
            r[level + "_ids"] = ";".join(ids)
            r[level + "_architecture"] = min(fwd, rev)
        out_loci.append(r)
        for h, f, prod in zip(hs, r["family_ids"].split(";"), r["products"].split(";")):
            d = fam[f]
            d["proteins"].add(h)
            d["loci"].add(r["locus_id"])
            d["carriers"].add(r["carrier"])
            d["species"].add(species[r["carrier"]])
            d["architectures"].add(r["family_architecture"])
            d["lists"].add(r["category"])
            d["homology_group"].add(assign["homology_group"][h])
            d["products"].add(prod)
            for c in LISTS:
                if r["category"] == c:
                    d["loci_" + c].add(r["locus_id"])
                    d["carriers_" + c].add(r["carrier"])
    ftab = []
    for f, d in sorted(fam.items()):
        row = dict(family_id=f, homology_groups=";".join(sorted(d["homology_group"])), proteins=len(d["proteins"]),
                   members_failing_check=sum(1 for h in d["proteins"] if not check[("family", h)][4]), loci=len(d["loci"]),
                   carriers=len(d["carriers"]), species=len(d["species"]), architectures=len(d["architectures"]),
                   products=" | ".join(sorted(d["products"]))[:300])
        for c in LISTS:
            row["loci_" + c] = len(d["loci_" + c])
            row["carriers_" + c] = len(d["carriers_" + c])
        ftab.append(row)
    ftab.sort(key=lambda x: (-x["loci"], x["family_id"]))
    write(a.out / "families.tsv", ftab)
    write(a.out / "candidate_loci_families.tsv", out_loci)
    pm = []
    for h, m in sorted(members.items()):
        row = dict(protein_hash=h, length=len(seqs[h]), loci=len(m), members=";".join("%s:%s:%s" % x for x in m))
        for level, name in (("family", "family"), ("homology_group", "homology")):
            rep_, fid, qc, tc, ok = check[(level, h)]
            row.update({name + "_id": assign[level][h], name + "_rep": rep_, name + "_identity": round(fid, 3),
                        name + "_member_cov": round(qc, 3), name + "_rep_cov": round(tc, 3), name + "_member_pass": ok})
        pm.append(row)
    write(a.out / "protein_members.tsv", pm)
    summ = dict(version=VERSION, mmseqs=ver, settings=MMSEQS_SETTINGS, member_min=MEMBER_MIN, check_search=MMSEQS_CHECK,
                distinct_proteins=len(seqs),
                member_check={lv: dict(members=sum(1 for r in pm), passed=sum(1 for r in pm if r[n + "_member_pass"]),
                                       fail=sorted(r["protein_hash"] for r in pm if not r[n + "_member_pass"]))
                              for lv, n in (("family", "family"), ("homology_group", "homology"))})
    for c in LISTS + ["all"]:
        rs = [r for r in out_loci if c == "all" or r["category"] == c]
        summ[c] = dict(loci=len(rs), carriers=len({r["carrier"] for r in rs}), species=len({species[r["carrier"]] for r in rs}),
                       cds_units=sum(int(r["cds_units"]) for r in rs), distinct_proteins=len({h for r in rs for h in loc_prots[r["locus_id"]]}),
                       families=len({f for r in rs for f in r["family_ids"].split(";")}),
                       homology_groups=len({f for r in rs for f in r["homology_group_ids"].split(";")}),
                       family_architectures=len({r["family_architecture"] for r in rs}))
    (a.out / "families_summary.json").write_text(json.dumps(summ, indent=1) + "\n")
    print(json.dumps(summ, indent=1))


def cmd_seeds(a):
    files = sample_files(a.pairs)
    rows = sorted((r for c in LISTS for r in read(a.s4 / ("list_%s.tsv" % c))), key=lambda r: (r["carrier"], r["locus_id"]))
    fam = {r["locus_id"]: r for r in read(a.out / "candidate_loci_families.tsv")} if (a.out / "candidate_loci_families.tsv").exists() else {}
    out, cat = a.work / "seed_catalog.fa", []
    seqs, cur = None, None
    with open(str(out), "w") as h:
        for r in rows:
            if r["carrier"] != cur:
                seqs, cur = read_fasta(files[r["carrier"]]["fna"]), r["carrier"]
            seq = seqs[r["contig"]]
            s, e = int(r["core_start"]), int(r["core_end"])
            core = seq[s:e]
            fid = fam.get(r["locus_id"], {}).get("family_ids", "")
            h.write(">%s|%s|%s:%d-%d|%s|%s\n%s\n" % (r["locus_id"], r["carrier"], r["contig"], s, e, r["category"], fid, core))
            cat.append(dict(locus_id=r["locus_id"], carrier=r["carrier"], contig=r["contig"], core_start=s, core_end=e, category=r["category"],
                            family_ids=fid, core_sha256=hashlib.sha256(core.encode()).hexdigest(),
                            U=seq[max(0, s - CONTEXT):s], D=seq[e:e + CONTEXT]))
    write(a.out / "seed_catalog.tsv", cat)
    print(len(cat), "seed records;", len({c["core_sha256"] for c in cat}), "distinct core sequences; fasta", out, sha256(out))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("mode", choices=["review", "families", "seeds"])
    p.add_argument("--loci", type=Path)
    p.add_argument("--calls", type=Path)
    p.add_argument("--pairs", type=Path, required=True)
    p.add_argument("--s4", type=Path, help="review output directory (families/seeds)")
    p.add_argument("--work", type=Path)
    p.add_argument("--mmseqs", default="mmseqs")
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    dict(review=cmd_review, families=cmd_families, seeds=cmd_seeds)[a.mode](a)


if __name__ == "__main__":
    main()
