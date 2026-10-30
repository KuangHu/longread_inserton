# Read downloads for the core-context check, batch by batch

Plan made on 2026-10-30 by `plan_batches.py` (run from the repo root, no network). Downloads run only on lrc-xfer.
Compute only reads runs in the ready manifest (see the repository README).

## What is planned

- **Starting point:**
  - 151 assemblies have verified reads.
  - 186 more assemblies have usable single-file runs in the pool that have not been downloaded.
  - Of the 6,661 cores in `match_both_lt2500_pairs_b1-2-3.tsv.gz`, 795 currently have a pair that can be tier A or B.
- **Selection:** 82 of the 186 assemblies would add a core or an A/B pair, for 42.7 GB in total. The other 104
  (118.4 GB) add neither and are not planned.
- **Order:** greedy, by the number of new cores that would get an A/B pair.
  - These counts are an upper bound: a pair becomes testable, not supported.
  - Support is known only after the read check.

| batch | assemblies | runs | GB | new cores (A/B possible) | new A/B pairs |
|---|---|---|---|---|---|
| 1 | 20 | 20 | 9.38 | 219 | 1,530 |
| 2 | 15 | 19 | 9.60 | 61 | 841 |
| 3 | 17 | 19 | 9.77 | 31 | 319 |
| 4 | 19 | 20 | 9.71 | 9 | 1,433 |
| 5 | 11 | 13 | 4.26 | 0 | 375 |

- **Batch tables:** `plan.tsv` has one row per run, with batch, rank, sample, platform, bytes and gain.
- **Input to prepare:** `batch<n>_runs.tsv` are rows of the pool `runs.tsv`, in the format `reads_ready.py prepare --runs-tsv` reads.
- **Re-planning:** if an earlier batch leaves runs unverified, rerun `plan_batches.py` before the next batch.
  Ready runs drop out of the plan.

## Each batch

1. **Storage check (any node, no download):**
   ```sh
   python scripts/reads_ready.py prepare --runs-tsv experiments/core_context/scale_comamonadaceae_v0.1/download_batches/batch1_runs.tsv --dry-run
   ```
   - It prints the reads directory size and the GB this batch still needs (partial files already on disk count
     as done).
   - It also prints the filesystem free space before and after the batch.
   - It refuses the batch (exit 4) if less than `--min-free-gb` (default 500) would remain.
2. **Download (lrc-xfer, OTP login, outside SLURM):** use the same command without `--dry-run`.
   - The storage check runs again first, and usage is printed again at the end.
   - Exit 1 lists runs that are not verified. Rerunning the same command resumes their `.part` files.
3. **Check the result:**
   - The last line of prepare gives verified / requested runs.
   - Unverified runs are not used by compute.
   - Record the storage lines in the batch log (`batch<n>_prepare.log`).
4. **Compute (lr6):** start a new cumulative round on the ready runs (`ucd_readcheck.py select`, then the
   `genome` → `check` → `status` → `summarize` stages of `../../ucd_readcheck.sbatch`).
   - The formal summary is written only when the round is complete by id.
   - Until then, use `--allow-incomplete` for a progress snapshot.
5. **Next batch:** go on to the next batch only after checking storage again.

## Storage at planning time

- **Reads directory** `/global/scratch/users/kh36969/longread_callout/data/reads`: 226.5 GB in all.
  - 185.1 GB is verified files.
  - 37 GB is rejected copies kept for the record (`md5_mismatch_attempt*`).
  - The rest is other files.
- **`/global/scratch`:** about 2,078 TB free. `checkquota` shows no per-user scratch quota.
- **All five batches:** add 42.7 GB.
