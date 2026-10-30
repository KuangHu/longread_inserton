"""Download / compute separation: scripts/reads_ready.py, scripts/fetch_long_reads.py and the compute entry points."""

import contextlib
import hashlib
import io
import json
import os
import re
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import fetch_long_reads as flr  # noqa: E402
import reads_ready as rr  # noqa: E402

# Compute entry points: none may run a downloader (comment and echo lines are documentation, not commands).
COMPUTE_ENTRIES = [
    "experiments/callout_scale_v1/batch10/batch_sample.sbatch",
    "experiments/callout_scale_v1/batch50/batch_sample.sbatch",
    "experiments/pacbio_hifi_pilot/pilot.sbatch",
    "experiments/read_callout/01_inputs.sbatch",
    "experiments/read_callout/02_callout.sbatch",
    "experiments/read_callout/03_reads.sbatch",
    "experiments/fna_pool/s3_readcheck.sbatch",
    "experiments/fna_pool/s5_readcheck.sbatch",
    "experiments/fna_pool/stage2_readcheck.sbatch",
    "experiments/core_context/ucd_readcheck.sbatch",
]
DOWNLOADERS = re.compile(r"fetch_long_reads\.py|reads_ready\.py prepare|s3_readcheck\.py fetch|datasets download|"
                         r"prepare_assembly\.py|verify_assembly\.py|download_fna\.py|kingfisher|curl |wget |prefetch")


def write_run(reads, run, data=b"@r1\nACGT\n+\nIIII\n", status="verified", sample="S1"):
    d = reads / run
    d.mkdir(parents=True)
    f = d / (run + ".fastq.gz")
    f.write_bytes(data)
    v = dict(run=dict(run=run, sample_id=sample), local_file=str(f), bytes=len(data), md5=hashlib.md5(data).hexdigest(),
             sha256=hashlib.sha256(data).hexdigest(), status=status, checks={})
    (d / "verify.json").write_text(json.dumps(v))
    return f


class ReadyManifestTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.reads = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def test_missing_manifest_means_nothing_is_ready(self):
        write_run(self.reads, "R1")
        self.assertTrue(rr.not_ready(["R1"], self.reads)["R1"].startswith("ready_manifest_missing"))

    def test_ready_run_gives_its_file(self):
        f = write_run(self.reads, "R1")
        rr.index(self.reads)
        self.assertEqual(rr.not_ready(["R1"], self.reads, rehash=True), {})
        self.assertEqual(rr.require(["R1"], self.reads), {"R1": str(f)})

    def test_each_kind_of_unready_input(self):
        write_run(self.reads, "OK")
        write_run(self.reads, "MD5ONLY", status="md5_verified")
        gone = write_run(self.reads, "GONE")
        cut = write_run(self.reads, "CUT")
        changed = write_run(self.reads, "CHANGED")
        flipped = write_run(self.reads, "FLIPPED")
        rr.index(self.reads)
        gone.unlink()
        cut.write_bytes(b"@r1\nAC")
        v = self.reads / "CHANGED" / "verify.json"
        v.write_text(json.dumps(dict(json.loads(v.read_text()), status="mismatch")))
        flipped.write_bytes(b"@r1\nACGA\n+\nIIII\n")  # same size, other content
        bad = rr.not_ready(["OK", "MD5ONLY", "GONE", "CUT", "CHANGED", "FLIPPED", "NEVER"], self.reads, rehash=True)
        self.assertNotIn("OK", bad)
        self.assertEqual(bad["MD5ONLY"], "status_md5_verified")
        self.assertTrue(bad["GONE"].startswith("file_missing"))
        self.assertTrue(bad["CUT"].startswith("size_mismatch"))
        self.assertEqual(bad["CHANGED"], "verify_json_differs_from_manifest")
        self.assertEqual(bad["FLIPPED"], "sha256_mismatch")
        self.assertEqual(bad["NEVER"], "not_in_ready_manifest")
        self.assertNotIn("FLIPPED", rr.not_ready(["FLIPPED"], self.reads))  # size-only check without --rehash

    def test_require_stops_with_exit_3(self):
        write_run(self.reads, "R1")
        rr.index(self.reads)
        with mock.patch("sys.stderr"):
            with self.assertRaises(SystemExit) as cm:
                rr.require(["R1", "R2"], self.reads)
        self.assertEqual(cm.exception.code, rr.NOT_READY)
        self.assertEqual(cm.exception.reasons, {"R2": "not_in_ready_manifest"})

    def test_check_cli_reports_input_not_ready_and_never_loads_the_downloader(self):
        write_run(self.reads, "R1")
        rr.index(self.reads)
        code = ("import sys, runpy; sys.argv = ['reads_ready.py', 'check', '--reads', %r, '--runs', 'R1', 'R2'];\n"
                "try:\n    runpy.run_path(%r, run_name='__main__')\n"
                "except SystemExit as e:\n    print('exit', e.code, 'fetch_long_reads' in sys.modules, 'urllib.request' in sys.modules)\n"
                % (str(self.reads), str(ROOT / "scripts" / "reads_ready.py")))
        p = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True,
                           env=dict(os.environ, SLURM_JOB_ID="123"))
        self.assertIn("input_not_ready R2 not_in_ready_manifest", p.stderr)
        self.assertEqual(p.stdout.strip(), "exit 3 False False")


class DownloadGuardTests(unittest.TestCase):
    def test_downloads_refused_inside_slurm(self):
        with tempfile.TemporaryDirectory() as t, mock.patch.dict(os.environ, {"SLURM_JOB_ID": "42"}), \
                mock.patch("subprocess.run") as run:
            with self.assertRaises(SystemExit):
                flr.download("https://example.invalid/x.fastq.gz", Path(t) / "x.fastq.gz", 10)
            with self.assertRaises(SystemExit):
                rr.cmd_prepare(mock.Mock(runs=["R1"], runs_tsv=None, reads=Path(t)))
            run.assert_not_called()

    def test_compute_entries_never_call_a_downloader(self):
        for rel in COMPUTE_ENTRIES:
            for n, line in enumerate((ROOT / rel).read_text().splitlines(), 1):
                cmd = line.strip()
                if not cmd or cmd.startswith("#") or cmd.startswith("echo") or "input_not_ready" in cmd:
                    continue
                self.assertIsNone(DOWNLOADERS.search(cmd), "%s:%d runs a downloader: %s" % (rel, n, cmd))

    def test_compute_entries_check_inputs(self):
        for rel in COMPUTE_ENTRIES:
            text = (ROOT / rel).read_text()
            if rel.endswith(("02_callout.sbatch", "stage2_readcheck.sbatch", "ucd_readcheck.sbatch")):
                continue  # no read input / inputs checked in Python (verified_runs, run_info + require)
            self.assertRegex(text, r"reads_ready\.py check|input_not_ready", rel)


class StorageTests(unittest.TestCase):
    def test_batch_refused_without_space_and_partial_counts_as_done(self):
        with tempfile.TemporaryDirectory() as t:
            reads = Path(t)
            write_run(reads, "DONE")
            (reads / "R2").mkdir()
            (reads / "R2" / "R2.fastq.gz.unverified.part").write_bytes(b"x" * 40)
            runs = [dict(run="DONE", fastq_bytes=10 ** 12), dict(run="R2", fastq_bytes=100)]
            free = os.statvfs(t).f_bavail * os.statvfs(t).f_frsize
            out = io.StringIO()
            with contextlib.redirect_stdout(out):
                rr.check_space(runs, reads, 0)
            self.assertIn("0.00 GB still to download", out.getvalue())
            with mock.patch("sys.stdout"), mock.patch("sys.stderr"):
                with self.assertRaises(SystemExit) as cm:
                    rr.check_space(runs, reads, free / 1e9 + 1)
            self.assertEqual(cm.exception.code, rr.NO_SPACE)


class TempFileTests(unittest.TestCase):
    """A download keeps a temporary name until every check passes."""

    def run_fetch(self, payload, md5, content=True):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        out = Path(tmp.name)
        run = dict(run="R1", sample_id="S1", url="https://example.invalid/R1.fastq.gz", fastq_bytes=len(payload),
                   fastq_md5=md5, platform="OXFORD_NANOPORE", archive_base_count=4, archive_read_count=1)
        seen = []

        def fake_download(url, target, expected_bytes, attempts=8):
            seen.append(target.name)
            target.write_bytes(payload)

        def fake_content(result, target, run, checks):
            seen.append(Path(target).name)
            checks.update(spot_count_match=True, base_count_match=True)

        with mock.patch.object(flr, "download", fake_download), mock.patch.object(flr, "content_checks", fake_content):
            result = flr.fetch(run, out, content=content)
        return out / "R1", result, seen

    def test_final_name_only_after_all_checks(self):
        payload = b"@r1\nACGT\n+\nIIII\n"
        d, result, seen = self.run_fetch(payload, hashlib.md5(payload).hexdigest())
        self.assertEqual(seen, ["R1.fastq.gz.unverified", "R1.fastq.gz.unverified"])
        self.assertEqual(result["status"], "verified")
        self.assertEqual(sorted(p.name for p in d.iterdir()), ["R1.fastq.gz", "verify.json"])
        self.assertEqual(result["local_file"], str(d / "R1.fastq.gz"))

    def test_mismatch_never_gets_the_final_name(self):
        payload = b"@r1\nACGT\n+\nIIII\n"
        d, result, _ = self.run_fetch(payload, "0" * 32)
        self.assertEqual(result["status"], "mismatch")
        self.assertFalse((d / "R1.fastq.gz").exists())
        rr.index(d.parent)
        self.assertEqual(rr.not_ready(["R1"], d.parent), {"R1": "status_mismatch"})

    def test_md5_only_stays_unverified_until_content_check(self):
        payload = b"@r1\nACGT\n+\nIIII\n"
        d, result, _ = self.run_fetch(payload, hashlib.md5(payload).hexdigest(), content=False)
        self.assertEqual(result["status"], "md5_verified")
        self.assertTrue((d / "R1.fastq.gz.unverified").exists())
        with mock.patch.object(flr, "content_checks",
                               lambda r, t, run, c: c.update(spot_count_match=True, base_count_match=True)):
            done = flr.check_content(d)
        self.assertEqual(done["status"], "verified")
        self.assertEqual(done["local_file"], str(d / "R1.fastq.gz"))
        self.assertFalse((d / "R1.fastq.gz.unverified").exists())


if __name__ == "__main__":
    unittest.main()
