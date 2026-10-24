# 10 株受控批次（protocol T2）

运行：2026-10-23，LBL Lawrencium（pc_rubinlab / lr6），作业 [batch_sample.sbatch](batch_sample.sbatch)，
每株一个任务、最多 2 个并发；样本冻结于 [samples.tsv](samples.tsv)（来源见 `../next_batch_inputs.tsv`、`../NEXT_BATCH_10_AUDIT.md`）。
流程：reads 下载校验（ENA MD5/bytes/counts）→ GCA/GCF 选择（序列逐条一致才用 GCF 及其注释）→
`multisite-insertion run --methods PILER RS --subset all --threads 2` → validate → minimap2 map-ont（固定命令）
→ 审核版 filled-site 检查（`nearby_competition_v2`）→ 独立 checker。发现与边界算法未改动。

## 判定

**10/10 株全部终态 complete，所有阶段退出码 0；package validate 与审核 checker 全部 pass。**
唯一失败为节点文件系统故障（见下），定位清楚，无共享流程缺陷。按 protocol T2 的技术条件，可进入 50 株规模；
但 50 株面板尚未按来源/数据资格选取，需另行冻结后才能提交。

## 每株结果（records，不是独立 insertion 数）

“审核支持”= audited_full_span（含带冲突）；“有 read 支持的候选结构”，不是完全确认的 insertion。

| 株 | 原生家族（tnp） | 导出 records | 审核支持 | ≥400 bp 审核支持 / ≥400 bp | edge | 竞争歧义 | 单侧 | 其它 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| CT11-68 | 20 (2) | 142 | 139 | 63 / 66 | 0 | 0 | 3 | 0 |
| DD10 | 51 (24) | 293 | 288 | 190 / 195 | 5 | 0 | 0 | 0 |
| FW306-02-B | 59 (5) | 714 | 713 | 209 / 210 | 0 | 0 | 0 | 1 |
| LjRoot267 | 185 (59) | 1,361 | 1,351 | 1,024 / 1,032 | 9 | 0 | 0 | 1 |
| R-92578 | 98 (50) | 541 | 536 | 377 / 382 | 5 | 0 | 0 | 0 |
| R-92583 | 67 (13) | 470 | 434 | 174 / 208 | 22 | 7 | 7 | 0 |
| R-92561 | 68 (14) | 463 | 429 | 176 / 207 | 18 | 8 | 8 | 0 |
| R-91495 | 62 (15) | 488 | 455 | 168 / 197 | 12 | 21 | 0 | 0 |
| R-92973 | 62 (15) | 488 | 455 | 168 / 197 | 12 | 21 | 0 | 0 |
| Bs869T2 | 24 (0) | 242 | 242 | 92 / 92 | 0 | 0 | 0 | 0 |

- 全部 10 株使用 GCF 序列与注释（与 GCA 逐条序列一致）；无 read 达 secondary 上限。
- R-91495 与 R-92973 各项计数完全相同，已核实为不同输入：assembly SHA256 不同（8,958,119 vs 8,958,099 bp）、
  reads 来自不同 run、成员 reads 总数不同（198,990 vs 147,444）。二者属 100% ANI fungorum 克隆，
  候选集合与由组装结构决定的分类相同。其余克隆成员 R-92583/R-92561 接近（470/463 records）。
- fungorum 克隆四株 edge/歧义/单侧较多，与其 draft 组装的 contig 断点相符。
- Bs869T2 有 GCF 注释但 0 个 tnp-supported 家族；其 242 条 records 全部审核支持。tnp 不是资格条件。

## 资源（[resources.tsv](results/resources.tsv)，[batch_summary.tsv](results/batch_summary.tsv)）

- 每株总 wall 16–99 min；主要耗时为 ENA 下载（0.14–1.69 GB，4–79 min，约 0.3–1 MB/s）。
- callout 3–43 min（LjRoot267 最长），map ≤1.5 min，classify ≤1.2 min；峰值 RSS 远低于 64 GB 预算。
- 每株输出 0.06–0.83 GB，reads 0.14–1.69 GB；10 株合计约 12 GB（输出 2.6、reads 8.6、assemblies 0.7 GB）。
- BLAST+ 首次加载偶有数分钟（warm_tools 0–7.5 min），与此前冷启动现象一致。

## 失败与重跑

- DD10 第 1 次：节点 n0258.lr6 上 `makeblastdb` 无法加载共享库（Lustre `Error 108`），callout 失败；
  reads 与 assembly 已校验。作业改为对已准备的 assembly 先核对 hash 再复用（`edfe05b`），旧目录保留于集群 `batch10_failed/`。
- DD10 第 2 次：节点 n0083.lr6 scratch 探针 60 s 超时，按设计以 75 退出（未写任何结果）。
- DD10 第 3 次：完成。累计排除的故障节点：n0348、n0039、n0292、n0082、n0031、n0258、n0083（均未被调度系统 drain）。
- 记录：[failed_attempts/](results/failed_attempts/)。

## 归档范围

`results/<株>/`：stages、source commit、callout summary/BED/families/family_summary、审核 summary 与 member_audit、
validate/audited_check、映射命令与版本、assembly_check、reads verify。
`member_boundaries.json` 与逐 read `read_evidence.tsv` 体积较大（单株压缩后可达 18 MB），留在集群，
路径、SHA256 与字节数见 [omitted_large_files.tsv](results/omitted_large_files.tsv)；raw reads 与 PAF 同样留在集群。
