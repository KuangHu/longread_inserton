"""Read length and base-quality statistics of one verified FASTQ (no filtering).

Per read: length and mean Phred computed from mean error probability (not mean QV).
Writes JSON with counts, bases, length quantiles/N50, base-QV quantiles, and read
fractions at mean Q >= 10/20/30. Compares records and bases with verify.json.
"""

import argparse
import gzip
import json
import math
import statistics
from pathlib import Path

ERR = [10 ** (-q / 10) for q in range(94)]


def quantiles(hist, total, ps):
    out, acc, it = {}, 0, iter(ps)
    p = next(it)
    for q in range(len(hist)):
        acc += hist[q]
        while p is not None and acc >= p * total:
            out["p%02d" % round(p * 100)] = q
            p = next(it, None)
    return out


def main():
    a = argparse.ArgumentParser(description=__doc__)
    a.add_argument("--verify", type=Path, required=True)
    a.add_argument("--out", type=Path, required=True)
    args = a.parse_args()
    v = json.loads(args.verify.read_text())
    if v["status"] != "verified":
        raise SystemExit("reads not verified: " + v["status"])
    lengths, read_q, base_hist = [], [], [0] * 94
    with gzip.open(v["local_file"], "rb") as h:
        while True:
            head = h.readline()
            if not head:
                break
            seq, _, qual = h.readline().rstrip(b"\r\n"), h.readline(), h.readline().rstrip(b"\r\n")
            if len(qual) != len(seq):
                raise SystemExit("sequence/quality length mismatch near " + head[:60].decode(errors="replace"))
            counts = [0] * 94
            for c in qual:
                counts[min(c - 33, 93)] += 1
            err = sum(n * ERR[q] for q, n in enumerate(counts) if n)
            for q, n in enumerate(counts):
                base_hist[q] += n
            lengths.append(len(seq))
            read_q.append(-10 * math.log10(max(err / len(seq), 1e-10)) if seq else 0.0)
    total = sum(lengths)
    ordered = sorted(lengths, reverse=True)
    acc = 0
    for n50 in ordered:
        acc += n50
        if acc * 2 >= total:
            break
    rq = sorted(read_q)
    res = dict(
        file=v["local_file"], sha256=v["sha256"], records=len(lengths), bases=total,
        matches_verify=len(lengths) == v["summary"]["records"] and total == v["summary"]["bases"],
        length=dict(min=ordered[-1], median=statistics.median(ordered), mean=round(total / len(ordered), 1),
                    n50=n50, max=ordered[0], p05=ordered[-1 - len(ordered) // 20]),
        base_qv=dict(mean=round(sum(q * n for q, n in enumerate(base_hist)) / total, 2),
                     **quantiles(base_hist, total, [0.05, 0.25, 0.5, 0.75, 0.95]),
                     fraction_qv93=round(base_hist[93] / total, 4), placeholder_qv0=base_hist[0]),
        read_mean_q=dict(min=round(rq[0], 2), median=round(statistics.median(rq), 2), max=round(rq[-1], 2),
                         fraction_ge_q10=round(sum(q >= 10 for q in rq) / len(rq), 4),
                         fraction_ge_q20=round(sum(q >= 20 for q in rq) / len(rq), 4),
                         fraction_ge_q30=round(sum(q >= 30 for q in rq) / len(rq), 4)),
        definition="read mean Q = -10 log10(mean base error probability); no reads filtered")
    args.out.write_text(json.dumps(res, indent=2) + "\n")
    print(json.dumps(res, indent=2))
    if not res["matches_verify"]:
        raise SystemExit("record/base totals differ from verify.json")


if __name__ == "__main__":
    main()
