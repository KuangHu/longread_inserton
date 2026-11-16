# Paraburkholderia round 1：四株单株 call-out、v1 read 证据与同种株间比较（归档）

运行日期：2026-10-23，LBL Lawrencium（account pc_rubinlab，partition lr6）。
面板在发现前冻结于 commit ee22986（[samples.tsv](samples.tsv)）。
本目录是 round-1 已完成运行的归档；扩容验收按
[CALLOUT_SCALE_PROTOCOL.md](../../docs/CALLOUT_SCALE_PROTOCOL.md) 在
[callout_scale_v1](../callout_scale_v1/) 另行执行，不覆盖这里的结果。

## 单株产量（records，不是独立 insertion 或事件数）

| 株 | 原生家族（tnp-supported） | 导出 records | 精确坐标组 | ≥400 bp records | 其中 v1 full≥2 | 分子深度 | boundary_source |
|---|---|---:|---:|---:|---:|---:|---|
| SMT4a | 62 (26) | 333 | 268 | 211 | 210 | 61.8× | provisional 316，censored 14，unresolved 3 |
| LjRoot223 | 44 (20) | 287 | 246 | 150 | 150 | 14.1× | provisional 276，censored 4，unresolved 7 |
| R-93128 | 63 (17) | 491 | 449 | 203 | 173 | 86.7× | provisional 446，censored 29，unresolved 16 |
| R-91497 | 43 (9) | 491 | 442 | 131 | 131 | 198.1× | provisional 472，censored 11，unresolved 8 |

“精确坐标组”只合并完全相同坐标的 records；跨工具非完全重合记录仍未解析，不是 locus 数。
四株 package validate 全部 pass；`scripts/check_callout_bundle.py --snapshot` 对本目录四株均为
`snapshot_tables_consistent`（仅小表一致性，不复核序列或 reads）。

## 同种株间比较（辅助观察，不是单株产量）

[compare_strains.py](compare_strains.py)（[04_compare.sbatch](04_compare.sbatch)）把两株候选 X 做严格家族对应
（单 HSP 双向覆盖 ≥90%、identity ≥95%，组内两两满足），株内 ≥90% 互重叠合并为 locus，
再用 nucmer --maxmatch 把每个 locus 的 500 bp 侧翼放到另一株；未唯一定位记 unresolved。
合成两基因组真值测试（共享 2、对侧空 2、对侧异序列 1）全部正确。

| 组 | 两株均有的家族 X≥400 bp（tnp） | locus 合计 | 共享位置 | 对侧无 X 且本株 v1 支持 | unresolved |
|---|---:|---:|---:|---:|---:|
| terricola：SMT4a＋LjRoot223 | 24 (7) | 96 | 48 | 12 | 36 |
| fungorum：R-93128＋R-91497 | 15 (1) | 43 | 34 | 2 | 7 |

**12 / 2 是跨株比较量，不能替代单株产量，也未经对侧 reads 检查空位点**
（“对侧无 X”目前只来自对侧 assembly）。<400 bp 家族两组分别 52 / 40 个，未产生任何新位置。
terricola 新位置主要来自 IS5（1.34 kb，6 个；0.86 kb，3 个）及 IS256、ISNCY、DDE 各 1 个。
unresolved 主要是侧翼在对侧无唯一定位（未比对上或多处定位）。

## 已知问题与运行记录

- **下载损坏**：LjRoot223 与 SMT4a 首次下载长度正确但 MD5 不符；与 curl 重新下载的正确文件比较，
  各只有一段约 198 bp 连续字节损坏（偏移 62.45 MB / 187.35 MB，非块对齐）。连同此前 PaAR2038，
  原因归于脚本的 urllib 下载路径。`scripts/fetch_long_reads.py` 已改用 curl（断点续传、重试），
  MD5 不符时保留坏文件并重下最多 3 次；改后在 LjRoot223 上端到端验证通过。
  坏文件 verify 记录见 [download_audit/](results/download_audit/)。
- **短 read 过滤**：v1 改名步骤丢弃 <500 bp 片段（LjRoot223：18,146 条中 8,462 条，约 1% 碱基）；
  进入比对的 reads 99% 可比对。新规范不再额外删除短 reads，重算须用新目录。
- **SMT4a 7 个 full_span_and_empty_allele**：对应 1 个 1,360 bp tnp-supported 位点（4 条重叠 records，
  2/61 分子缺 X）及 3 个 81–102 bp RS 片段（小缺失噪声）。SMT4a 为 PCR 文库，低频等位型不作结论。
- 节点 n0348.lr6 / n0039.lr6 此前有文件系统异常，本轮作业均排除。
- 资源：callout 6–16 min；reads 1–8 min（16 CPU）；比较 3–5 min（4 CPU），见各 `resources.txt`。

## 文件

`results/<株>/`：`callout/`（summary、validate、candidates.bed、families、family_summary、
member_boundaries.json.gz、resources）、`evidence/`（v1 member_read_support、深度曲线、非参考表、summary）、
`assembly_check.json`（GCA/GCF 逐条序列 SHA256 一致、文件 hash）、`reads_verify.json`、`read_unit.json`、
`rename_stats.json`、minimap2 命令与版本。`results/compare_<组>/`：loci、position_comparisons、
family_increment、family_pair_matches、summary、resources。
raw reads、PAF、FASTQ 中间文件与 nucmer/blast 工作目录留在集群 `$PROJ`，不入库。
