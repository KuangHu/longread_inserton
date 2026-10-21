"""Check recovered sequence and coordinates against the artificial truth."""
import json
from pathlib import Path
import sys

from site_recovery.core import identity

truth = {r["contig"]: r for r in json.loads(Path(sys.argv[1]).read_text())}
records = [json.loads(line) for line in Path(sys.argv[2]).read_text().splitlines()]
accepted = [r for r in records if r["status"] == "read_supported_pair"]
assert len(accepted) == 2, [(r["source_call"], r["status"], r["reasons"]) for r in records]
for rec in accepted:
    t = truth[rec["reference"]["contig"]]
    assert rec["reference"]["start0"] == t["start0"]
    assert rec["reference"]["end0"] == t["end0"]
    assert rec["empty_site_insertion_pos0"] == t["empty_site_insertion_pos0"]
    assert rec["sequences"]["empty_site"] == t["empty_site"]
    score = identity(rec["sequences"]["insertion"], t["insertion"])
    assert score >= 0.995, score
    print(json.dumps({"contig": t["contig"], "platform": rec["platform"],
                      "insertion_identity": score, "n_occupied": rec["support"]["n_occupied"],
                      "n_empty": rec["support"]["n_empty"], "coordinates_correct": True}))
