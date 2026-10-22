"""Input validation, immutable input staging, and command provenance."""

import hashlib
import json
import re
import shutil
import subprocess
import time
from pathlib import Path

from Bio import SeqIO


def read_json(path):
    return json.loads(Path(path).read_text())


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + "\n")


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1048576), b""):
            h.update(block)
    return h.hexdigest()


def load_manifest(path):
    path = Path(path).resolve()
    data = read_json(path)
    if (
        data.get("schema_version") != 1
        or not isinstance(data.get("samples"), list)
        or not data["samples"]
    ):
        raise ValueError(
            "Manifest requires schema_version=1 and a nonempty samples list"
        )
    seen = set()
    samples = []
    for value in data["samples"]:
        sample = dict(value)
        name = sample["id"]
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]*", name) or name in seen:
            raise ValueError("Sample IDs must be unique filename-safe labels: " + name)
        seen.add(name)
        sample["fasta"] = str((path.parent / sample["fasta"]).resolve())
        records = fasta_dict(sample["fasta"])
        if not records or any(len(r) == 0 for r in records.values()):
            raise ValueError("Empty FASTA or sequence: " + sample["fasta"])
        if any(
            not re.fullmatch(r"[ACGTRYKMSWBDHVNacgtrykmswbdhvn]+", str(r.seq))
            for r in records.values()
        ):
            raise ValueError(
                "FASTA must contain ungapped IUPAC DNA: " + sample["fasta"]
            )
        sample["annotations"] = [
            str((path.parent / f).resolve()) for f in sample.get("annotations", [])
        ]
        if sample.get("annotation_format", "embl") not in ("embl", "genbank"):
            raise ValueError("annotation_format must be embl or genbank")
        for f in sample["annotations"]:
            if not Path(f).is_file():
                raise ValueError("Missing annotation: " + f)
        sample["native"] = {
            k: str((path.parent / v).resolve())
            for k, v in sample.get("native", {}).items()
        }
        if set(sample["native"]) - {"PILER", "RS"}:
            raise ValueError("native keys must be PILER or RS")
        samples.append(sample)
    return samples


def stage_inputs(samples, out):
    out = Path(out).resolve()
    for sample in samples:
        for source in sample.get("native", {}).values():
            if out == Path(source).resolve() or Path(source).resolve() in out.parents:
                raise ValueError("Output cannot be inside an imported native directory")
    if out.exists() and any(out.iterdir()):
        raise ValueError(
            "Output directory must be new or empty; existing runs are never overwritten"
        )
    out.mkdir(parents=True, exist_ok=True)
    staged = []
    for original in samples:
        s = dict(original)
        folder = out / "inputs" / s["id"]
        folder.mkdir(parents=True)
        target = folder / "genome.fasta"
        shutil.copyfile(s["fasta"], target)
        s["source_fasta"] = s["fasta"]
        s["fasta"] = str(target.relative_to(out))
        s["fasta_sha256"] = digest(target)
        s["source_annotations"] = s["annotations"]
        s["annotations"] = []
        s["annotation_sha256"] = {}
        for i, source in enumerate(original["annotations"]):
            target = folder / f'annotation_{i}.{s.get("annotation_format", "embl")}'
            shutil.copyfile(source, target)
            rel = str(target.relative_to(out))
            s["annotations"].append(rel)
            s["annotation_sha256"][rel] = digest(target)
        staged.append(s)
    write_json(out / "inputs/manifest.json", {"schema_version": 1, "samples": staged})
    return staged


def samples_in_run(out):
    samples = read_json(Path(out) / "inputs/manifest.json")["samples"]
    for s in samples:
        for rel, expected in [
            (s["fasta"], s["fasta_sha256"]),
            *s["annotation_sha256"].items(),
        ]:
            if digest(Path(out) / rel) != expected:
                raise ValueError("Staged input changed: " + rel)
    return samples


class Runner:
    def __init__(self, out, bins=None):
        self.out = Path(out).resolve()
        self.bins = bins or {}
        self.path = self.out / "provenance.json"
        self.provenance = (
            read_json(self.path)
            if self.path.exists()
            else {"commands": [], "executables": {}}
        )

    def executable(self, name):
        path = shutil.which(self.bins.get(name, name))
        if not path:
            raise ValueError(f"Missing executable {name}; see docs/INSTALL.md")
        path = str(Path(path).resolve())
        self.provenance["executables"][name] = {"path": path, "sha256": digest(path)}
        return path

    def run(self, label, name, args, cwd=None):
        argv = [self.executable(name), *map(str, args)]
        work = Path(cwd or self.out).resolve()
        log = self.out / "logs" / (label + ".log")
        log.parent.mkdir(parents=True, exist_ok=True)
        entry = {"label": label, "argv": argv, "cwd": str(work), "status": "running"}
        self.provenance["commands"].append(entry)
        write_json(self.path, self.provenance)
        start = time.monotonic()
        try:
            with log.open("w") as handle:
                result = subprocess.run(
                    argv, cwd=work, stdout=handle, stderr=subprocess.STDOUT, check=True
                )
            entry.update(status="complete", returncode=result.returncode)
        except subprocess.CalledProcessError as error:
            entry.update(status="failed", returncode=error.returncode)
            raise RuntimeError(f"{label} failed; see {log}") from error
        finally:
            entry["wall_seconds"] = round(time.monotonic() - start, 3)
            write_json(self.path, self.provenance)
        print(label + " complete", flush=True)


def fasta_dict(path):
    with Path(path).open() as handle:
        return SeqIO.to_dict(SeqIO.parse(handle, "fasta"))


def completed_blast(path, queries):
    with Path(path).open() as handle:
        count = sum(line.startswith("# Query: ") for line in handle)
    if count != queries:
        raise ValueError(f"Incomplete BLAST output: {path}: {count}/{queries} queries")
