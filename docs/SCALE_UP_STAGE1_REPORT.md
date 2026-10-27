# 扩量第一阶段报告：assembly 候选目录

> **2026-10-27 主线纠偏：** 发现目标是“相似 core、不同上下文”，不要求空位点。Gap 搜索仅作辅助。
> 当前主线见 [纠偏计划](PLAN_CORE_CONTEXT_DISCOVERY.md) 和 [可运行的 core-context 原型](../experiments/core_context/README.md)。
> 下文已完成的 fna_pool Stage 1/2 是 **gap-derived 辅助目录与已知位点验证**，不是新主线完成。
> 新入口已通过 H16/H850 真实 smoke（旧 R33 六位置找回），但 GFF3 完整性、独立蛋白入口和精细上下文分组仍待接入。


计划：[PLAN_SCALE_UP_STAGE1_ASSEMBLY_DISCOVERY.md](PLAN_SCALE_UP_STAGE1_ASSEMBLY_DISCOVERY.md)（d1c663a）。
状态：**`assembly_catalog_ready_for_long_reads`**。本阶段没有新的 long-read 下载、比对或验证任务。

**证据等级：** 本报告中的全部候选都是 *assembly-derived candidate loci*，证据来源为 `assembly_only`，
long-read evidence pending。它们来自 assembly 两两比较和 GFF3 注释复核，没有用 reads 检查过。
main loci 不是独立插入事件，也不是新的 insertion 家族；蛋白家族是初步蛋白家族。
本阶段没有做 S6 位点合并，所以没有“位点组数”（site groups），下面的计数都是 loci。

计数单位必须分开：loci（候选位点）、样本（BioSample；本轮一个 BioSample 只保留一个 assembly）、
蛋白（实例 / accession / 不同序列）、初步蛋白家族、架构。organism 字符串只计作 organism label 数，不代表物种数。

## 1. 范围与冻结

- 10 个 P1 科（wave 1），每科一个 batch，取全部 isolate 级 assembly；MAG/SAG/atypical 不纳入。
  批次定义与来源见 [target_batches.tsv](../experiments/fna_pool/scaleup_stage1_20261026/target_batches.tsv)，
  运行前快照见 [baseline.json](../experiments/fna_pool/scaleup_stage1_20261026/baseline.json)。
- 跨 batch 重复：没有 BioSample 或 accession 同时出现在两个 batch
  （[cross_batch_overlap.tsv](../experiments/fna_pool/scaleup_stage1_20261026/cross_batch_overlap.tsv) 为空）。
  catalog 的可追溯性检查也再次确认了这一点。
- 与旧 Comamonadaceae 基线池相比，485 个样本的 assembly 不变，新增 1,220 个，去掉 14 个
  （9 个 origin_metagenome，5 个 ncbi_atypical；见 [baseline_pool_delta.tsv](../experiments/fna_pool/scaleup_stage1_20261026/baseline_pool_delta.tsv)）。
  没有因为缺 reads 而排除任何 assembly。
- 各步参数沿用已冻结的 S2/S3/S4 设置；第一阶段 catalog 的参数、软件版本、代码 commit、输入与输出哈希见
  [frozen_manifest.json](../experiments/fna_pool/scaleup_stage1_20261026/stage1_catalog/frozen_manifest.json)
  和 [outputs.sha256](../experiments/fna_pool/scaleup_stage1_20261026/stage1_catalog/outputs.sha256)。

## 2. 逐 batch 计数

来源：[stage1_report.json](../experiments/fna_pool/scaleup_stage1_20261026/stage1_catalog/stage1_report.json)，
它汇总各 batch 的 S1/S2/S3/S4 输出和 catalog，没有重新计算。

### 2.1 Assembly、下载、比较、callout

| batch | assemblies | 项目 | 有注释 / 缺注释 | 下载失败 | S2 样本 | 未选到 pair | pairs | raw calls | S3 loci | 短 anchored loci |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Actinomycetaceae | 1,139 | 322 | 1,132 / 7 | 0 | 1,139 | 23 | 3,611 | 1,291,940 | 450,224 | 52,667 |
| Comamonadaceae | 1,705 | 894 | 1,682 / 23 | 0 | 1,705 | 14 | 4,399 | 3,579,886 | 1,250,815 | 125,528 |
| Eggerthellaceae | 329 | 128 | 315 / 14 | 0 | 329 | 4 | 1,004 | 674,748 | 209,687 | 16,341 |
| Erwiniaceae | 2,500 | 950 | 2,476 / 23 | 1 | 2,499 | 13 | 8,710 | 3,058,673 | 964,029 | 105,616 |
| Flavobacteriaceae | 3,704 | 2,054 | 3,682 / 22 | 0 | 3,704 | 165 | 7,531 | 4,292,683 | 1,695,818 | 107,571 |
| Lachnospiraceae | 3,254 | 669 | 2,968 / 284 | 2 | 3,252 | 141 | 8,425 | 4,420,424 | 1,429,768 | 261,573 |
| Metamycoplasmataceae | 1,738 | 396 | 1,561 / 177 | 0 | 1,738 | 17 | 5,418 | 427,601 | 140,480 | 11,845 |
| Microbacteriaceae | 3,444 | 1,550 | 3,266 / 178 | 0 | 3,444 | 20 | 8,259 | 4,771,431 | 1,809,841 | 209,674 |
| Micrococcaceae | 2,975 | 1,364 | 2,958 / 16 | 1 | 2,974 | 11 | 8,031 | 4,480,277 | 1,490,247 | 171,615 |
| Nitrobacteraceae | 1,975 | 1,115 | 1,971 / 4 | 0 | 1,975 | 4 | 6,262 | 6,890,874 | 2,224,868 | 295,417 |
| **合计** | **22,763** | 9,442* | 22,011 / 748 | **4** | **22,759** | 412 | **61,650** | **33,888,537** | **11,665,777** | **1,357,847** |

\* 项目数是逐 batch 相加；同一个 BioProject 可能跨 batch，所以这个合计不是去重数。

- 每个 batch 中 assemblies 与 BioSamples 数相同（一个 BioSample 一个 assembly），合计都是 22,763。
- 下载：43,735 个下载条目中 43,731 个核验通过，4 个连续 3 次失败，已排除，没有替换：
  - Lachnospiraceae GCF_902652075.1、GCF_902652265.1：package MD5 不符；
  - Micrococcaceae GCA_977014115.1：注释文件不完整；
  - Erwiniaceae GCA_003568535.1：FNA 总长与 assembly report 不一致。

  它们没有进入 S2，S2 样本为 22,759。
- 可检查的 GCA/GCF 配对序列全部一致。另有 2 对因下载失败没有检查。
- 缺注释（748 个 carrier）的位点进入 annotation_pending，不进入候选清单（§3）。
- 有 412 个样本没有选到比较对象，状态为 `no_pair_selected_pending_review`；它们不计作零命中。

### 2.2 候选类别、蛋白、家族

| batch | main | partner_mixed | direction_mixed | short_other_segment | replacement | annotation_pending | 不同蛋白序列 | 局部家族 | 全局家族† | 局部架构 | 翻译不符记录 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Actinomycetaceae | 373 | 11 | 3 | 574 | 303 | 114 | 993 | 682 | 682 | 543 | 7 |
| Comamonadaceae | 1,010 | 35 | 2 | 2,697 | 502 | 853 | 3,673 | 2,460 | 2,460 | 1,979 | 24 |
| Eggerthellaceae | 117 | 5 | 0 | 169 | 66 | 855 | 290 | 231 | 231 | 174 | 2 |
| Erwiniaceae | 1,072 | 28 | 8 | 1,097 | 505 | 481 | 1,983 | 1,144 | 1,144 | 874 | 11 |
| Flavobacteriaceae | 2,312 | 75 | 18 | 3,186 | 1,080 | 246 | 5,859 | 4,847 | 4,847 | 3,936 | 10 |
| Lachnospiraceae | 3,306 | 109 | 14 | 2,368 | 537 | 22,840 | 4,764 | 3,180 | 3,179 | 2,578 | 52 |
| Metamycoplasmataceae | 510 | 12 | 4 | 105 | 22 | 314 | 409 | 226 | 226 | 190 | 4 |
| Microbacteriaceae | 1,737 | 34 | 9 | 6,753 | 1,148 | 9,068 | 8,806 | 6,441 | 6,442 | 5,169 | 21 |
| Micrococcaceae | 1,626 | 64 | 22 | 3,664 | 651 | 582 | 5,094 | 3,552 | 3,552 | 3,001 | 39 |
| Nitrobacteraceae | 4,351 | 216 | 49 | 4,866 | 1,460 | 885 | 6,024 | 3,890 | 3,890 | 3,277 | 14 |
| **合计 / 全局** | **16,414** | **589** | **129** | **25,479** | **6,274** | **36,238** | 37,895 | 26,653（局部相加） | **26,643**（全局去重） | — | **184** |

† 某 batch 的 loci 用到的全局家族数。各 batch 的全局家族数不能直接相加，因为有 10 个全局家族跨两个 batch。

## 3. S4 类别闭合

本阶段共复核 **1,357,847 个短 anchored loci**（core < 2,500 bp），类别互斥：

| 类别 | loci |
|---|---:|
| main | 16,414 |
| partner_mixed + direction_mixed | 589 + 129 = 718 |
| short_other_segment | 25,479 |
| replacement | 6,274 |
| annotation_pending | 36,238 |
| boundary_intersection | 655,138 |
| other_structure（0 个或 >2 个完整 CDS） | 520,025 |
| boundary_or_repeat_review | 88,836 |
| gene_integrity_anomaly | 8,725 |
| **合计** | **1,357,847** |

- 合计与短 anchored loci 总数一致（`s4_categories_close: true`）。
- 边界诊断标记可以叠加，但类别本身互斥。
- 五类候选共 48,885 loci，加上 annotation_pending 36,238，共 85,123 loci。
  这些数直接从最终的 `candidate_core_catalog.tsv.gz` 统计，core FASTA 也是 85,123 条。
- annotation_pending 约占全部短 anchored loci 的 2.7%，在 Lachnospiraceae 内约 8.7%（22,840 / 261,573）。
  它们是待审查位点，补齐注释后不一定成为 main；它们不参与蛋白聚类。

## 4. 编码检查与翻译表

- CDS 按 GFF3 中的 `transl_table` 翻译。表 4 与表 11 已实现；缺失、冲突或不受支持的表进入审查，不会一律按表 11 处理。
  本轮没有出现 transl_table 问题记录。
- 翻译检查共有 361,291 条 match、16,015 条 pseudo、**184 条 mismatch**。
  - 这 184 条是 coding-unit 检查记录，分布在 180 个 loci，不等于 184 个独立基因。
  - 对应 loci 的类别：boundary_intersection 88，gene_integrity_anomaly 96。没有一条进入 main 或其他候选清单。
  - 逐条明细见 [coding_anomalies.tsv](../experiments/fna_pool/scaleup_stage1_20261026/stage1_catalog/coding_anomalies.tsv)。
- 因此不能沿用旧基线“全部翻译一致”的结论。

## 5. 蛋白目录与初步家族（跨 batch 统一）

| 计数 | 值 |
|---|---:|
| 蛋白实例（locus × 编码单元） | 60,747 |
| protein accession | 39,065 |
| 不同蛋白序列（完整序列 SHA256） | 37,895 |
| 初步家族（90% identity，双侧 80% coverage） | 26,643 |
| 同源分组（50%，双侧 80%；仅作候选同源） | 20,973 |
| 家族架构（按 core 顺序和方向，两链取规范形式） | 21,714 |

- **聚类方式：** 10 个 batch 的全部蛋白合并后一次聚类，使用冻结的 S4 设置（MMseqs2 17.b804f）。
  然后每个成员与其最终代表比对，保存 identity、双向 coverage 和是否通过。
  在 Eggerthellaceae 上，这种做法与原 S4 的全对全检查得到完全相同的簇和数值。
- **成员验收：** 90% 家族 37,895 / 37,895 通过。
  50% 同源分组中，WP_349626926.1（Nitrobacteraceae，111 aa）对原代表的 identity 为 0.619，coverage 为 0.748 / 0.706，不合格。
  按同一阈值在全部代表中寻找合格归属，结果为 0 个，因此它单独成组（H1_77004f244708）。
  原失败记录保存在 [member_check_repairs.tsv](../experiments/fna_pool/scaleup_stage1_20261026/stage1_catalog/member_check_repairs.tsv)，
  `protein_catalog` 的 `homology_repair` 列也有标记。修复后两个层级全部通过；阈值没有改，90% 家族没有重跑。
- **局部到全局：** 各 batch 局部家族相加为 26,653，全局为 26,643。
  - 一对一 26,632；跨 batch 合并 18；同 batch 内合并 2（Lachnospiraceae F_0064 与 F_1932）；拆分并跨 batch 合并 1（Microbacteriaceae F_0686）。
  - 没有蛋白序列跨 batch 出现。共有 10 个全局家族跨两个 batch，主要是 Micrococcaceae 与 Microbacteriaceae 或 Actinomycetaceae。
  - 映射见 [local_global_family_map.tsv](../experiments/fna_pool/scaleup_stage1_20261026/stage1_catalog/local_global_family_map.tsv)。
- **ID 规则：** 全局 ID 绑定目录版本和代表序列哈希（`F1_` / `H1_` / `A1_` + sha256 前 12 位）。
  locus、样本和局部家族一律写成 `<batch>:<id>`，避免跨 batch 重名覆盖。
  以后新增 batch 或代表发生变化时，要提供新旧 ID 映射。

## 6. Core 与侧翼目录

[candidate_core_catalog.tsv.gz](../experiments/fna_pool/scaleup_stage1_20261026/stage1_catalog/candidate_core_catalog.tsv.gz)
每个 locus 一行，共 85,123 行（五类候选 + annotation_pending），保存：
- 实际观察到的 core 坐标（0-based 半开区间，carrier 正链）、长度、SHA256；
- 两侧各 300 bp 的 U/D 序列及其哈希，以及用于从完整 assembly 扩展窗口的坐标；
- contig 长度、contig SHA256、环状标记（没有 GFF 时为 `unknown_no_gff`）；
- 重复和近边缘标记、partner 及其 assembly、逐 partner 争议（partner_breakdown）、蛋白哈希、全局家族与架构、原始翻译检查；
- `evidence_source=assembly_only`。

core 序列本身（共 63,501 种不同序列）每个 locus 一条，存于 scratch FASTA；保存路径、SHA256 和重建方法见 frozen_manifest。
目录不用代表序列或共识序列替代各个 DNA 实例。

## 7. 旧 read 证据关联审计

对象是 Comamonadaceae 基线 read check 的 464 个 loci 和 743 个 test。明细见
[pilot_read_evidence_links.tsv](../experiments/fna_pool/scaleup_stage1_20261026/stage1_catalog/pilot_read_evidence_links.tsv)
和 [pilot_read_test_links.tsv](../experiments/fna_pool/scaleup_stage1_20261026/stage1_catalog/pilot_read_test_links.tsv)。

- **精确关联 231 个**：BioSample 相同；assembly accession.version 相同或 FNA SHA256 相同；contig、core 坐标和 core SHA256 都相同。
  这 231 个在第一阶段中是 main 229 个、partner_mixed 2 个。
  第一次 catalog 运行报告 0 个匹配，原因是整数和字符串坐标比较的 bug，已修复（4c381f7）。
  修复后重跑（27315389）和最终运行都是 231。
- **未精确匹配 233 个**，只写文件能证实的原因：

  | 原因 | loci |
  |---|---:|
  | carrier 不在第一阶段池（被去掉的 14 个样本之一） | 11 |
  | 相同 core 坐标出现在第一阶段 S4，但类别为 boundary_or_repeat_review | 85 |
  | 相同 core 坐标出现在第一阶段 S3，状态为 unresolved（core 1,096 / 876 bp） | 2 |
  | **unresolved**：同一 contig 上有重叠但不相同的第一阶段 loci（已列出） | 88 |
  | **unresolved**：没有重叠的第一阶段 locus | 47 |

  88 个重叠情况没有被认定为“边界变化”，原因仍未查明。
- **test 层面：**
  - 268 个属于同一 carrier–partner test：精确关联，且该 partner 也是第一阶段 locus 的 partner。
  - 114 个只关联到同一 carrier 位点：旧 partner 不在新 locus 的 partner 中，这些 test 的 partner 侧结论不能转给新 partner。
  - 361 个无关联。
- 旧证据只停留在对应的 locus 和 test 上，不传播到重叠位点、其他 partner、家族或架构；新候选默认为 assembly_only。

## 8. 第二阶段接口（只有元数据，没有执行）

| 文件 | 内容 |
|---|---|
| `assembly_run_links.tsv.gz` | assembly–run 关联、平台与类型、关联依据 |
| `read_run_manifest.tsv.gz` | 4,414 个 run：可用性、大小、MD5、本地核验状态（162 个已核验在本地；SRR21808314 带有 repair 记录） |
| `read_batch_manifest.tsv.gz` | 223,309 行 locus × carrier/partner：关联 runs、可用 runs、本地已核验 runs、状态、远端预计字节、大小未知的 runs |
| `read_evidence.tsv`、`new_read_contexts.tsv` | 只有表头的空表，列定义在代码中；以后按 locus/test/run 追加，不覆盖 assembly 发现证据 |

各类候选的 reads 覆盖（“可用”= 至少一个 usable run，本地或远端均可）：

| 类别 | loci | carrier 有可用 reads | 其中 partner 也有 | 仅 partner 有 | 都没有 | carrier 本地已核验 | carrier 仅远端 | carrier 已关联但不可用 | carrier 无关联 run |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| main | 16,414 | 3,159 | 1,088 | 1,621 | 11,634 | 204 | 2,955 | 221 | 13,034 |
| partner_mixed | 589 | 153 | 76 | 110 | 326 | 10 | 143 | 11 | 425 |
| direction_mixed | 129 | 23 | 12 | 24 | 82 | 1 | 22 | 0 | 106 |
| short_other_segment | 25,479 | 3,240 | 851 | 2,792 | 19,447 | 404 | 2,836 | 265 | 21,974 |
| replacement | 6,274 | 709 | 143 | 601 | 4,964 | 52 | 657 | 58 | 5,507 |
| annotation_pending | 36,238 | 321 | 0 | 6,527 | 29,390 | 0 | 321 | 506 | 35,411 |
| 全部 | 85,123 | 7,605 | 2,170 | 11,675 | 65,843 | 671 | 6,934 | 1,061 | 76,457 |

- annotation_pending 中 carrier 有可用 reads 的 321 个 loci，其 partner 都没有关联 run，所以“partner 也有”为 0。
- **main 清单的下载估计：** 1,506 个去重 run（carrier 1,231 个，另有只作为 partner 用到的 275 个），
  约 1.22 TB（1,221,551,833,200 bytes）。同一 run 被多个 loci 使用时只计一次；大小未知的 run 为 0；本地已核验、无需下载的 run 有 146 个。
- **全部候选：** 1,908 个去重 run，约 1.56 TB。
- 本阶段开始的下载数为 0。

## 9. 可追溯性检查（`catalog_summary.json → traceability`，全部通过）

- `global_locus_id` 唯一（85,123）。每个 BioSample 只在一个 batch 中，每个 assembly 只对应一个 sample_key。
- 每个 locus 的 sample_key 和 assembly 都能回到 S2 catalog。
- 每个蛋白实例和 locus 蛋白哈希都在 protein_catalog 中；每个实例的 locus 都在 core 目录中。
- 局部家族和局部同源 ID 都带 batch 命名空间，并出现在 local_global_family_map 中。
- 两个层级全部成员通过最终代表检查，代表本身也是成员。core FASTA 每个 locus 一条（85,123）。
- **scratch 文件：** frozen_manifest 记录路径、字节数、SHA256、产生它的模式和重建命令。这些文件包括：
  - proteins.faa、protein_instances；
  - candidate_cores.fa；
  - 两个层级的 cluster 表与成员比对；
  - 修复搜索结果。
- **输入：** 记录每个 batch 用到的 S2/S3/S4 文件的路径和 SHA256，以及基线 read check 文件的哈希。

## 10. 文件名与计划名称的对应

| 计划中的名称 | 实际文件 |
|---|---|
| assembly_catalog | `<batch>/s2/catalog.tsv`（排除项见 `catalog_excluded.tsv`） |
| download_manifest / annotation_files | `<batch>/download_manifest.tsv`、`<batch>/s1/files.tsv.gz`、`<batch>/s1/accessions.tsv.gz`；S4 的 `annotation_files.tsv` |
| assembly_aliases（GCA/GCF） | `<batch>/s1/pair_checks.tsv.gz` |
| 距离与 pair 目录 | `<batch>/s2/pairs.tsv`、`groups.tsv`、`skani_run.json` |
| calls 与 loci | `<batch>/s3_full/`（仓库内为 gzip 表与 sha256；完整表在 scratch） |
| S4 清单 | `<batch>/s4_s3_full/list_*.tsv`、`s4_loci.tsv.gz` |
| 蛋白家族与架构目录 | `stage1_catalog/protein_catalog.tsv.gz`、`families.tsv`、`architectures.tsv`、`local_global_family_map.tsv`、`member_check_repairs.tsv` |
| candidate_core_catalog | `stage1_catalog/candidate_core_catalog.tsv.gz` |
| 第二步 manifests | `stage1_catalog/assembly_run_links.tsv.gz`、`read_run_manifest.tsv.gz`、`read_batch_manifest.tsv.gz`、`read_evidence.tsv`、`new_read_contexts.tsv` |

## 11. 验收

| # | 条件 | 结果 |
|---|---|---|
| 1 | batch 与 accession 范围冻结，没有因 reads 缺失而排除 assembly | 通过 |
| 2 | 冻结条目都有状态；FNA 全部核验；缺失注释有去向 | 通过（4 个失败下载排除并记录原因；748 个缺注释 carrier → annotation_pending） |
| 3 | 计划内 pair 都成功或明确排除 | 通过（61,650 对全部完成；412 个样本 no_pair_selected_pending_review） |
| 4 | 回归与针对性测试 | 通过（S4 在基线上回归一致；跨 origin 的 frameshift 有测试；member check 新方法与 S4 全对全一致） |
| 5 | locus 可追溯到 assembly、contig、calls、注释和 core 序列 | 通过 |
| 6 | 类别分开且合计闭合 | 通过（1,357,847） |
| 7 | 家族成员对最终代表达到阈值；各计数单位分开统计 | 通过（修复 1 个 50% 层级成员后） |
| 8 | 跨 batch 映射去重，局部与全局计数的关系可解释 | 通过 |
| 9 | 交付第二阶段接口，没有 long-read 下载或搜索 | 通过（下载 0） |

## 12. 限制

- 全部候选只有 assembly 证据。main 不代表经过 reads 支持，也不代表独立插入事件。
- 76,457 / 85,123 个 loci 的 carrier 没有关联的 long-read run。第二阶段对这些 loci 只能使用 partner 的 reads 或另找数据。
- annotation_pending（36,238）需要先补注释或另行审查，才能参与家族和 main 判断；不能假定它们会成为 main。
- 233 个旧 read check loci 没有精确对应，其中 135 个原因未查明（unresolved）。
- 初步蛋白家族只基于完整 CDS 的蛋白，不覆盖非编码或部分编码的 core。
  50% 同源分组只是候选同源，不能当作家族。
- 环状标记只来自 GFF 的 `Is_circular`。没有标记不等于线性；没有 GFF 的 contig 状态未知。
- 项目数逐 batch 相加时没有去重；物种数没有给出，organism label 数不能当作物种数。
- 没有 site groups：S6 合并要在 read 证据之后进行。
