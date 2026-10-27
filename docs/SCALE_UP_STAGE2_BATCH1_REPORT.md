# 扩量第二阶段第一批：本地 reads 的已知位点检查

输入：第一阶段目录 `23f195b`（[SCALE_UP_STAGE1_REPORT.md](SCALE_UP_STAGE1_REPORT.md)）。
本批**只检查已知位点**：检查 carrier 的 present 和 partner 的 absent 是否得到该样本自己的 long reads 支持。
本批**没有**用 core 在 reads 中搜索新位点；普通 core 命中不算新增 site。**没有下载新的 reads。**
第一阶段的候选、类别和 ID 都没有改动。read 证据作为新的一层保存，assembly 证据保持原样（`assembly_only` 不被覆盖）。

这些 reads 通常就是当初组装该 assembly 的数据，所以本检查是 assembly 一致性检查，用来发现组装错误、坍缩重复或样本内等位混合，
不是独立测序验证。

结果目录：[experiments/fna_pool/scaleup_stage1_20261026/stage2_readcheck_batch1/](../experiments/fna_pool/scaleup_stage1_20261026/stage2_readcheck_batch1/)；
代码：[stage2_readcheck.py](../experiments/fna_pool/stage2_readcheck.py)、[stage2_readcheck.sbatch](../experiments/fna_pool/stage2_readcheck.sbatch)、
[s3_readcheck.py](../experiments/fna_pool/s3_readcheck.py)（v0.3）。

## 1. 冻结目标

- 规则：main 清单中，carrier 至少有一个 usable run 已在本地核验过（第一阶段 `read_batch_manifest`）。
- 共 **204 个 loci**，与第一阶段报告的 204 一致，全部属于 Comamonadaceae（只有这个 batch 有本地 reads）。
  目标清单和输入哈希见 `targets.tsv` / `targets.json`。
- 每个 locus 的全部 calls 都形成 test：共 **320 个 carrier–partner tests**，涉及 208 个样本。
  partner 有本地 reads 时一并检查；没有时保留单侧结果（`no_reads`），不当作不支持。

## 2. 本地 reads 核查（`runs_checked.tsv`）

逐一检查 208 个样本的全部关联 run：verify 记录状态、BioSample、平台、文件是否存在、文件大小，并**重新计算 SHA256**（共 113 GB）。

| 状态 | runs |
|---|---:|
| local_verified_ok（实际使用） | 130 |
| not_local（partner：SRR40356297、SRR40356277、SRR40356274，ONT） | 3 |
| linked_not_usable（carrier SAMN11265848 的 SRR31313491，PacBio） | 1 |
| 记录失效、文件缺失、大小或哈希不符 | 0 |

- 只有通过核查的 run 被比对（`runs_allowed.tsv`）。
- 不在本地或不可用的 run 单独列出，不计入任何判定。
- SRR21808314 是第一阶段经多数投票修复的文件，本次核查通过（verify 中保留 repair 记录）。
- 读段类型与 minimap2 预设沿用 v0.3：ONT 97 个 run 用 map-ont；PacBio HiFi/CCS 3 个用 map-hifi；CLR 8 个和类型未定的 PacBio 22 个用 map-pb。

## 3. 方法（v0.3，没有调整阈值）

判定规则沿用 v0.3，未作修改：
- 读段必须完整跨越两侧连接点并各留 300 bp；
- CIGAR 中的插入缺失检查；
- 全基因组竞争：两个 anchor 必须定位到本位点；
- 分子去重：ZMW、SRA spot、ONT duplex；
- 按长度匹配计数；
- `MIN_SPAN=3`，`MAX_OTHER_FRACTION=0.1`。

本批只给 v0.3 增加了两项输出，计数方法不变：
1. 每个 分子 × test × 角色 × run 一行支持记录（`support_records.tsv.gz`），内含 read 名和逐 read 判定；
2. 可选的 run 白名单。

**回归：** 用修改后的代码对旧 pilot 重新 summarize，`readcheck_tests.tsv`、`readcheck_loci.tsv`、`readcheck_summary.json` 与已提交版本逐字节一致（job 27318897）。
S3/S4 单元测试分别通过 24/24、15/15。

## 4. 结果

### 4.1 Tests（320 个）

| 类别 | tests |
|---|---:|
| 双侧一致（carrier present 且 partner absent） | 150 |
| 仅 carrier 一致 | 167 |
| 仅 partner 一致 | 1 |
| 覆盖不足（有 reads 的一侧跨越分子少于 3 个） | 2 |
| 重复定位不清导致不足 | 0 |
| mixed | 0 |
| contradicted | 0 |
| 未测试（两侧都没有本地 reads） | 0 |

两侧判定组合：

| carrier / partner | tests |
|---|---:|
| consistent / consistent | 150 |
| consistent / no_reads | 165 |
| consistent / insufficient | 2 |
| insufficient / consistent | 1 |
| insufficient / no_reads | 2 |

- 167 个仅 carrier 一致的 test 中，165 个是因为 partner 没有本地 reads，2 个是 partner 有 reads 但跨越分子不足。
- 被排除的分子–test 组合：
  - CIGAR 不合格 486；
  - 重复定位不清（locus_ambiguous）129；
  - 定位到别处（off_locus）68；
  - 同一分子内各条记录不一致 0。

  这些都在判定之前排除，不计作支持或反对。

### 4.2 Loci（204 个）

| 证据等级 | loci |
|---|---:|
| E2（所有 test 两侧一致） | 85 |
| E2s（部分 test 两侧一致） | 26 |
| E1c（仅 carrier 一致） | 90 |
| E1p（仅 partner 一致） | 1 |
| E0（未能判定） | 2 |
| X（mixed / contradicted） | 0 |

E1c 只说明 carrier 的 assembly 与它自己的 reads 一致；没有 partner reads，就不能说明 partner 确实缺少该片段。

### 4.3 样本与 runs

- carrier 样本 97 个，全部有本地 reads。
- 只作 partner 的样本 111 个，其中 23 个有本地 reads。另有 41 个样本在不同 test 中既是 carrier 又是 partner。
- 实际使用 130 个 run：ONT 97，PacBio 未定类型 22，CLR 8，HiFi/CCS 3。
- 计数单位：无法按分子归组的 run 以 read record 计，各 test 都标出了 `count_unit`。

### 4.4 异常清单（`anomalies.tsv`）

| 类型 | 个数 | 内容 |
|---|---:|---|
| 有 reads 但跨越分子不足 | 5 | 涉及 SAMEA115676864，在 3 个 test 中作 carrier、1 个中作 partner，跨越分子为 0；另有 SAMEA115676832 作 partner 的 1 个 test |
| run 不在本地 | 3 | 均为 partner |
| run 已关联但不可用 | 1 | 一个 carrier run |

没有 mixed、contradicted，也没有记录失效。

## 5. 与旧 pilot 的比较（`pilot_comparison.tsv`）

对象是第一阶段精确关联的 382 个 pilot tests（同一 carrier–partner test 268 个，只同一 carrier 位点 114 个）：

- **146 个在本批重新检查：** 两侧判定和 test 判定全部相同，所有计数逐项相同。
  输入 run、test 几何和参数都没有变化，所以结果一致。
- **236 个不在本批，各有原因：**
  - 192 个：该 locus 的 carrier 没有本地 reads，不属于本批目标；pilot 中这些 test 的 carrier 也是 no_reads。
  - 42 个：locus 在本批，但旧 partner 不是第一阶段 locus 的 partner。旧 partner 的证据不转给新 partner。
  - 2 个：第一阶段类别为 partner_mixed，不在 main 目标中（pilot 中为 E2）。
- 旧证据只用于比较，本批结论都来自实际重新运行。

## 6. Read evidence 层

`read_evidence.tsv` 沿用第一阶段定义的列，每个 test × 角色一行，共 640 行。每行记录：
- 所用 runs、count_unit、各类计数、判定、locus 证据等级、是否冲突；
- 所属 test 和 locus 的全局 ID（`comamonadaceae:L_…`）。

每个判定都能通过 `support_records.tsv.gz` 追到具体的 run、分子和 read 记录。第一阶段的 `candidate_core_catalog` 和类别都没有改写。

## 7. 资源（`resources.json`）

| 步骤 | 作业 | 墙钟时间 | CPU 时间 | 最大内存 |
|---|---|---:|---:|---:|
| 回归 | 27318897 | 6.2 min | 0.10 h | 3.7 GB |
| 准备：选择、等位、本地核查（重算 113 GB 的 SHA256） | 27318898 | 2.6 min | 0.14 h | 12.2 GB |
| 比对（208 个 array 任务） | 27319067 | 合计 109 min，单任务最长 8.4 min | 3.92 h | 8.1 GB |
| 汇总 | 27319068 | 4.7 min | 0.05 h | 2.4 GB |

- 代码 commit 记录在每个作业日志中（dadb281）。layer 步骤在汇总之后补充了比较原因和异常类别，并重新运行（见本次提交）。
- 比对产生的 PAF 在 scratch：`/global/scratch/users/kh36969/longread_callout/fna_pool/scaleup_stage1_20261026/stage2_readcheck_batch1/paf/`。
  可用 `stage2_readcheck.sbatch map`（array）重建；仓库中保存每个样本的 run 与参数记录（`paf_json/`）。

## 8. 下一批下载清单（只有元数据，未下载）

范围：`next_batch_download_list.tsv`，覆盖除本批 204 个之外的 main loci；每个 run 只计一次，计入它所属的最高优先级。

| 优先级 | loci | 新增 runs | 预计字节 |
|---|---:|---:|---:|
| 1. carrier 和 partner 都有可用 reads | 973 | 600 | 500.4 GB |
| 2. 仅 carrier 有可用 reads | 1,982 | 736 | 585.7 GB |
| 3. 仅 partner 有可用 reads | 1,621 | 170 | 135.5 GB |
| 合计 | 4,576 | 1,506 | 1.22 TB |

大小未知的 run 为 0。建议下一批先下载优先级 1。

## 9. 下一项任务（尚未开始）

用 core 引导的 long-read 新位点搜索：用 core 实例在 reads 中检索，产生新的上下文（`new_read_contexts`），
再经 assembly 定位和复核。这是独立模块，与已知位点检查分开；v0.3 不能描述为已经实现新增位点搜索。
