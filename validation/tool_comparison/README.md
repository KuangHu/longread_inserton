# 四工具实测的复现入口

本目录的代码和说明发布到GitHub；`results/`目前仅保留本地，等待明确的数据上传授权。下文的结果文件路径用于说明本地记录位置。

目标是比较完整 insertion 与边界证据；empty-site 支持单独记录，不是 INS 的接受门槛。
这是一个真实 CCS 数据的有限 pilot，不是已知真值 benchmark，也不是 Nanopore benchmark。
结论见 [实测报告](../../docs/TOOL_COMPARISON.md)。

## 输入

复用原 pilot 的 `hifi_marine_pilot/`：完整 SRR28395891 FASTQ、418 MAG 竞争比对 BAM、参考和冻结的 10 MAG 选择。
下载及对齐步骤见 [原 pilot](../../docs/REAL_HIFI_PILOT.md)。不要用单独对 10 MAG 重比对的 MAPQ 冒充原竞争比对结果。

以下命令从工作目录运行，仓库目录名是 `longread_inserton/`，工作输出是 `tool_comparison/`：

```bash
python longread_inserton/validation/tool_comparison/prepare.py hifi_marine_pilot tool_comparison/inputs
python longread_inserton/validation/tool_comparison/alias_contigs.py tool_comparison/inputs
```

`prepare.py` 保留原始 reads，选取 primary MAPQ ≥20 的 read 名用于 read-only 工具；BAM 保留目标 contig 上的原始 alignment 记录和 MAPQ。参考只重命名为 `chr1...chr235` 以兼容 DeBreak 的命名假设，`alias_map.json` 可逆。
因此工具仍有不同的内部 read 选择策略，不能宣称所有内部 evidence 完全相同。
focal 子集保留与 C11 的 [80000,110000) 重叠的完整 reads。该区域按既知 1.95 kb 案例选择，仅用于诊断。

## 依赖和工具版本

精确上游 commits、软件包版本见 `results/provenance.json`。上游源码未更改；DIVE 的两个旧接口别名由 wrapper 提供兼容。

| 工具 | 本次源码commit |
|---|---|
| TLDR | `46d9718591505f66dad4a76c8419549c9b054fed` |
| DeBreak | `21d9b4294089ccd9df0feec976239551f4359096` |
| Palidis4 | `d85894ebae269a8dbf37df6897cdd8d6e14cf2b3` |
| 旧Palidis | `f1105ec9cfbdc7fe5cd7ee6cf16f1c30afee3981` |
| DIVE / biodive | `180964c10b2c8cb1a05111741f09b3a28d50d82b` |
| pal-MEM | `fe8f5fb87ba9dc8d907f3cc0f2e33e9a8cd18701` |
| bsalign | `6cff5337b7be6588bd79a5bf388b1e64841d503d` |

本机使用 micromamba 环境，主要包为 Python 3.11、pysam、numpy、scipy、biopython 1.81、matplotlib、statsmodels、scikit-learn、editdistance、bx-python、tqdm、psutil、samtools、minimap2、mafft、exonerate、wtdbg、bowtie2、cd-hit、blast、prodigal、libboost-headers。
TLDR 与 biodive 从冻结的 Git commit 进行 `pip install --no-deps -e` 安装。bsalign 用 `make -j2` 编译；pal-MEM 在其目录用以下命令编译：

```bash
g++ -std=c++11 -O3 -pthread -I /path/to/env/include pal-mem.cpp -o pal-mem
```

将环境、bsalign 所在目录加入 PATH。运行 `compare.py` 的 Python 还需安装本仓库及其 edlib/pyspoa 依赖。

## 实际运行命令

每次尝试（包括失败）均在 `results/runs/<name>/run.json` 保存完整参数、cwd、退出码、时间和内存上限，日志在同目录。
工具源码放在 `tool_comparison/tools/`。按 `provenance.json` 的 remote/commit checkout 后，复现相应 `command` 数组。

`run_logged.py` 的内存限制是每进程虚拟地址空间限制，不是进程组总 RSS；`max_rss_kib` 是系统 child-resource 指标，不能拿它做严格并行总内存 benchmark。不同尝试曾并行运行，耗时也仅供本机诊断。

主要比较参数：

| 流程 | 参数/范围 |
|---|---|
| TLDR | `-e none -m 3 --min_te_len 500 --max_te_len 50000`；最后全范围配置额外 `-p 1 --max_cluster_size 20`，不启用detail输出，4 GiB地址空间上限；50-read/detail配置中断，不能计作完成 |
| DeBreak | `--min_size 500 --max_size 50000 --min_support 3 --min_quality 20 -t 2 --no_genotype --poa --rescue_large_ins` |
| DeBreak 深度诊断 | 同上，加 `--maxcov 1000000000`，仅诊断全局深度过滤影响，不是推荐生产阈值 |
| DIVE | 原默认统计参数；wrapper 固定随机种子 42；先427-read focal，再60415-read targets，3 GiB 上限 |
| Palidis IR 模块 | pal-MEM 原生 `-fu` unpaired 输入，`-t 2 -l 25 -m 50 -k 15`；随后原生 clip_reads、Bowtie2；不等同完整 Palidis4 |

Palidis 模块后续步骤：

```bash
python tool_comparison/tools/palidis4/bin/clip_reads.py \
  --read_fasta tool_comparison/runs/palmem/targets_IR.fasta \
  --tab_file tool_comparison/runs/palmem/targets_IR.tab \
  --output_prefix tool_comparison/runs/palmem/targets
bowtie2-build --threads 2 -q tool_comparison/inputs/aliases.fasta tool_comparison/inputs/bowtie_targets
bowtie2 --very-sensitive-local -x tool_comparison/inputs/bowtie_targets \
  -U tool_comparison/runs/palmem/targets_irs.fasta \
  -S tool_comparison/runs/palidis_ir_mapping/irs.sam -p 2 -f
python longread_inserton/validation/tool_comparison/summarize_ir_module.py tool_comparison
```

`summarize_ir_module.py` 是本项目的探索性后处理，不是 Palidis 官方完整调用器：只保留原生 Bowtie2 报告的 primary、无剪切、无 gap、零 mismatch 的 IR 命中，再配对反向互补末端，区间长 500–3000 bp。不枚举所有多重命中，可能漏检；不能把结果视作 Palidis 的 recall。

最后汇总：

```bash
python longread_inserton/validation/tool_comparison/compare.py \
  tool_comparison hifi_marine_pilot tool_comparison/comparison
python longread_inserton/validation/tool_comparison/check_diagnostics.py tool_comparison hifi_marine_pilot
python longread_inserton/validation/tool_comparison/collect_results.py \
  tool_comparison longread_inserton/validation/tool_comparison/results
python longread_inserton/validation/tool_comparison/update_report.py \
  tool_comparison longread_inserton/docs/TOOL_COMPARISON.md
```

## 指标的含义

- 原生 candidate、PASS、PRECISE 与本项目 read 核验分开。
- 500–50000 bp 是本次筛选范围，不是 IS 分类标准。
- INS read 核验：primary、MAPQ ≥20；距候选 ≤100 bp；长度差 ≤25%；100 bp 双侧 anchor 各 ≥95% identity；read 内完整 insertion 与候选 ≥95% 全长 identity。允许边界略有偏移并报告范围，不强制精确坐标相同。
- `audit_three_occupied_reads` 对 INS 要求至少3个满足上述条件的 read 名；不要求 empty reads。对 DEL 只说明至少3条 reads 跨越参考中的 occupied 区间，**不能独自验证删除边界**；另看 `audit_three_occupied_and_gap_reads`。
- 核验未处理 split-only 支持，也未做局部重组装或 strain phasing。未通过不能直接判作假阳性；通过也不是单碱基边界真值或 IS 分类。
- DEL 是 reads 相对参考删除的区间，可能提供 reference-carried 插入候选，不能和非参考 INS 数量直接混为一谈。
- `candidate_sequences.fasta` 包括未通过核验的原生候选，header 标记证据状态。DeBreak/Sniffles flanks 是参考上下文，TLDR flanks 是 occupied consensus；都不自动变成真实 empty allele。
- `reference_like_or_empty_gap_reads` 是局部支持计数，没有完成真实 empty allele 的全长拼接、SNP polishing 或跨 MAG 唯一性证明。
- 交集按同 contig/同类型、断点差≤100 bp、长度比≥0.75，仅表示 caller agreement；不是 precision/recall。

仓库不保存原始 FASTQ/BAM，也不修改原 pipeline 的正式 accepted 结果。

TLDR额外focal诊断使用同一个C11区间：`focal.chroms`内容为`chr186`，`focal.bed`内容为`chr186\t80000\t110000`（实际tab分隔）。该次未启用`--detail_output`；完整参数见`results/runs/tldr_focal/run.json`。它只检查既知位点，未加入全范围caller交集计数。
