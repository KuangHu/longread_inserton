"""Create artificial DNA for plumbing tests; not a biological benchmark."""
import argparse
import json
from pathlib import Path
import random


def rc(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def create(out):
    out.mkdir(parents=True, exist_ok=True)
    rng = random.Random(72)
    dna = lambda n: "".join(rng.choices("ACGT", k=n))
    references, truths, molecules = {}, [], []
    for kind in ("INS", "DEL"):
        left, right, insert = dna(15000), dna(15000), dna(1200)
        insert = rng.choice([x for x in "ACGT" if x != right[0]])+insert[1:-1]+rng.choice([x for x in "ACGT" if x != left[-1]])
        chrom = "synthetic_"+kind
        references[chrom] = left+(insert if kind == "DEL" else "")+right
        truths.append({"contig": chrom, "svtype": kind, "start0": len(left),
                       "end0": len(left)+(len(insert) if kind == "DEL" else 0),
                       "insertion": insert, "empty_site": left[-1000:]+right[:1000],
                       "empty_site_insertion_pos0": 1000})
        for state, hap in (("occupied", left+insert+right), ("empty", left+right)):
            # Uniform background coverage with many independent break-spanning reads.
            for j in range(100):
                start = rng.randrange(0, len(hap)-8000)
                sequence = hap[start:start+8000]
                if j % 2:
                    sequence = rc(sequence)
                molecules.append((f"{kind}_{state}_{j}", sequence))
    with (out/"assembly.fasta").open("w") as fh:
        for name, sequence in references.items():
            fh.write(f">{name}\n{sequence}\n")
    (out/"truth.json").write_text(json.dumps(truths, indent=2)+"\n")
    for platform, sub, indel in (("hifi", 0.001, 0.0001), ("ont", 0.015, 0.003)):
        with (out/f"{platform}.fastq").open("w") as fh:
            for name, seq in molecules:
                result = []
                for base in seq:
                    p = rng.random()
                    if p < indel/2:
                        continue
                    if p < indel:
                        result.append(rng.choice("ACGT"))
                    result.append(rng.choice([x for x in "ACGT" if x != base]) if rng.random() < sub else base)
                sequence = "".join(result)
                fh.write(f"@{name}\n{sequence}\n+\n{'I'*len(sequence)}\n")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("out", type=Path)
    create(p.parse_args().out)
