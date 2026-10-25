# FNA 数据池：S0 清单与关联、S1 批量下载

计划见 [PLAN_CROSS_ASSEMBLY_GFF3_SHORT_INSERTIONS.md](../../docs/PLAN_CROSS_ASSEMBLY_GFF3_SHORT_INSERTIONS.md) 的 S0–S1。
本轮决定（2026-10-24）：

- FNA 来源：NCBI Assembly，优先 RefSeq，GenBank 补充；
- 首个 ≥100 样本数据池：Comamonadaceae（taxid 80864）；
- 内容：所有与 ONT/PacBio 通过 exact BioSample 关联的 assembly，加上科内全部 RefSeq Complete Genome（assembly-only 背景）；
- 批量下载在 lrc-xfer 数据传输节点上运行。

脚本只用 Python 标准库（3.6 可运行）和 NCBI datasets CLI。原始响应和下载文件留在
`/global/scratch/users/kh36969/longread_callout/fna_pool/<pool>/`，git 只存表格、校验值和小结。

## S0：[build_pool.py](build_pool.py)

```sh
W=/global/scratch/users/kh36969/longread_callout/fna_pool/comamonadaceae_20261024
/usr/bin/python3 build_pool.py fetch --taxon 80864 --work $W
/usr/bin/python3 build_pool.py build --taxon 80864 --work $W --out comamonadaceae_20261024
```

- `fetch`：NCBI Datasets genome summary（RefSeq、GenBank 当前版本）；ENA read_run（ONT + PacBio）；
  PacBio run 的 SRA XML（原始提交文件名，作为 read type 证据）。每个原始文件的命令/URL、行数和 SHA256 在 `raw/provenance.json`，并抄入 `summary.json`。
- `build`：
  - 可用 run：GENOMIC、WGS、非 PCR、ENA FASTQ 有 MD5 和 bytes、平均读长 ≥500 bp（排除标成长读段的短读段数据，如 PacBio Onso、150–302 bp 的 R1/R2 文件）。
  - read type：ONT 按平台；PacBio 只按原始提交文件名判断 HiFi/CCS 或 CLR/subreads，ENA 的 `_subreads.fastq.gz` 是 SRA 命名，不作证据；判断不了标 `pacbio_unresolved`。
  - 关联：run 的 sample_accession 与 assembly BioSample 精确相等 → `biosample_matched`。是否是组装所用 reads（`assembly_source_run`）本轮**未核查**。
  - 每个 BioSample 一个 assembly：当前版本 → 组装级别 → 最新发布 → accession；有 RefSeq 配对且 NCBI 未标“不同”时下载 GCF，同时下载 GCA 的 FNA 做序列一致性核对。
  - 不做来源（临床/宿主）筛、物种排除或克隆合并；这些在分 batch 时处理。

### Comamonadaceae 结果（[summary.json](comamonadaceae_20261024/summary.json)）

| 项 | 数 |
|---|---:|
| NCBI 当前 assembly | RefSeq 2,400；GenBank 6,819 |
| ENA run | 422（ONT 267，PacBio 155） |
| 可用 run | 373：ONT 243，HiFi/CCS 21，CLR/subreads 17，PacBio 类型未定 92 |
| 可用且同 BioSample 有 assembly 的 run | 258 |
| 数据池 BioSample（= assembly） | 499：linked 236，RefSeq Complete 背景 263 |
| linked 样本按 read type | ONT 163，HiFi 21，CLR 12，PacBio 未定 40（一个样本可有多种） |
| linked assembly 级别 | Complete 155，Chromosome 5，Contig 76 |
| 序列来源 | RefSeq 485，GenBank 14 |
| BioProject | 295 |
| NCBI atypical 标记 | 5 |
| 下载清单 | 984 行：主包 499（FNA+GFF3+FAA+seq report），GCA 核对 485（FNA+seq report） |

表格：`pool.tsv`（每样本一行）、`download_manifest.tsv`、`assembly_run_links.tsv`、`runs.tsv`（含不可用原因）、`assemblies.tsv.gz`（全科目录）。
linked 的 Contig 级 assembly 有 76 个，跨 FNA 调用时边界可靠性需单独评估。

## S1：[download_fna.py](download_fna.py)

每块（默认 50 个 accession）：dehydrated 包 → `datasets rehydrate` → 逐 accession 核查 → 通过后才移入 `data/genomes/<accession>/`。
每次尝试追加到 `data/state.jsonl`；重跑只处理未通过的，失败最多重试 3 次；同一时间只允许一个实例（`data/.lock`）。

核查：包内 MD5 全部一致；FNA 记录与 NCBI sequence report（accession、长度）一致且全部保留；
assembly report 的 accession、BioSample 与清单一致；GFF3 `##sequence-region` 与 FNA 一致；
GFF3 每个 CDS 的 protein_id 都在 FAA，FAA 每个 id 都在 GFF3。没有注释记为 `annotation_missing`，不算失败。
`report` 再逐对比较 GCF 与 GCA 的记录序列 SHA256，写 `pair_checks.tsv`；不一致的 pair 不去重。

小试验（4 个 BioSample、6 个 accession，登录节点，系统 Python 3.6.8）：6/6 通过，2 对 GCF/GCA 序列一致
（含一对版本号不同：GCF_964252695.1 / GCA_964252695.2），49 MB；重跑不重复下载；不存在的 accession 记为 failed 并保留 NCBI 错误信息。

### 在 lrc-xfer 上运行

lrc-xfer 需要 OTP 交互登录，只能由用户启动。它与集群共享 home 和 `/global/scratch`。

```sh
ssh lrc-xfer.lbl.gov
tmux new -s fna        # 断线后用 tmux attach -t fna 回来
cd ~/middle_insertion_callout/longread_inserton/experiments/fna_pool
P=/global/scratch/users/kh36969/longread_callout/fna_pool/comamonadaceae_20261024
/usr/bin/python3 download_fna.py run --manifest comamonadaceae_20261024/download_manifest.tsv \
    --data $P/data 2>&1 | tee -a $P/download.log
# Ctrl-b d 脱离；中断后重跑同一命令即可续传
```

估计：约 7 GB（主包约 10 MB/个，GCA FNA 约 4 MB/个），按小试验速度约半小时。完成后在登录节点运行：

```sh
/usr/bin/python3 download_fna.py report --manifest comamonadaceae_20261024/download_manifest.tsv \
    --data $P/data --out comamonadaceae_20261024/s1
```

reads（S5）不在本步下载。linked 样本的可用 FASTQ 共约 213 GB（ONT 151，PacBio CLR 25，HiFi 12，类型未定 25），等候选和优先级明确后再分批下载。
