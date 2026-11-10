# H16/H850 GFF3 注释试验

日期：2026-10-24（America/Los_Angeles）。这是冻结候选坐标的注释审计，不是新的 insertion caller。
完整后续计划见 [跨 assembly 短 insertion 计划](../../docs/PLAN_CROSS_ASSEMBLY_GFF3_SHORT_INSERTIONS.md)。

## 实测结果

| 项目 | 数量 |
|---|---:|
| 输入 assembly presence/absence 候选 | 70 |
| core 长度严格小于 2500 bp | 44 |
| 通过本轮严格注释条件 | 14，全部完整包含 1 个 CDS |
| 需复核 | 30 |
| 与注释特征边界相交 | 25 |
| 完整包含 2 个 CDS | 4，全部仍需复核 |
| 14 个合格 CDS 的基因组翻译与 NCBI FAA 一致 | 14/14 |

原因可重叠，不应相加。严格条件是完整包含 1–2 个 CDS、相同数量的 gene、无其它相交边界、
无 RNA 重叠、无已标记的 partial/pseudogene/translation exception，且有匹配 protein_id 的 FAA。
这只是主要分析子集：保留全部 44 行和逐特征详情，不删除断裂基因或非典型候选。
严格子集不是已确认 IS 集合；普通宿主基因的有无差异也可通过这些注释条件。

## 来源和检查

| 样本 | Assembly | GFF annotation source |
|---|---|---|
| H16 | GCF_051903345.1 | GCF_051903345.1-RS_2025_08_09 |
| H850 | GCF_056236025.1 | GCF_056236025.1-RS_2026_03_24 |

通过 NCBI Genomes HTTPS 下载同一 assembly 目录的 FNA、GFF3、FAA 和 MD5 清单。
6 个压缩文件均匹配 NCBI MD5；7 条 contig 的长度及序列 SHA256 与冻结输入一致；
GFF assembly accession、sequence-region 与 FNA 完整对应。URL、文件 SHA256、注释日期在
[provenance.json](results/provenance.json)。注释可以在 assembly 版本不变时更新，复跑须对比这些哈希。

[输入候选表](inputs/assembly_candidates.tsv) 来自本地 H16/H850 pilot：minimap2 2.31-r1302，
asm10 PAF CIGAR I/D gap 400–20000 bp，primary MAPQ ≥30，两侧各 500 bp 全局 identity ≥90%。
它不是 PanGraph 图拓扑调用，也未判定插入/删除的历史方向。原始 assembly 冻结依据为
仓库 commit `7ced967de214b417e50ad46a5d6ee6a94c0a8df5`。
本文件夹复现范围从已保存的 70 行坐标开始，不声称重现整个 assembly caller 或原来的 reads 搜索。
400 bp 是旧发现步骤的下限；这次没有补搜 <400 bp 的短候选。

GFF 的 1-based 闭区间转换为 0-based 半开区间。同一 ID 的多行 CDS 算一个 feature，
不按 protein_id 合并不同位点，父特征的 partial/pseudo 标记会传播。边界相交不自动向 CDS 端点修剪。
当前脚本只用于这两份 NCBI 细菌注释，非通用 GFF importer；跨 circular origin 的虚拟坐标会报错，
需要独立处理，不能静默截断。GFF ID 只在本次文件中使用，长期身份使用 assembly/contig/坐标/哈希。

## 复跑

Python ≥3.10、Biopython；先下载一次，之后可离线重跑。不要以 Python `-O` 关闭校验断言。

```bash
python experiments/gff3_short_insertion_pilot/annotate_candidates.py --download
python -m unittest discover -s experiments/gff3_short_insertion_pilot -p test_annotation.py -v
python experiments/gff3_short_insertion_pilot/verify_eligible_translations.py
```

6 个测试覆盖区间端点、相邻基因不相交、1 bp 截断、负链多段 CDS、父特征 pseudo/partial、
缺失蛋白/RNA 和转义字符。额外从基因组翻译 14 个合格 CDS，全部与 FAA 一致。

输出：

- [candidate_annotations.tsv](results/candidate_annotations.tsv)：44 个候选的基因数、protein_id、产品标签和复核原因。
- [feature_details.json](results/feature_details.json)：原始特征坐标、属性、方向和重叠关系。
- [summary.json](results/summary.json)、[translation_checks.json](results/translation_checks.json)：计数和翻译检查。

本轮没有重新验证 reads、推断转座活性或运行蛋白聚类。`eligible_1_2_CDS` 不能替代这些证据。
原始 FNA/GFF3/FAA 下载保留在被忽略的 `data/`，GitHub 只存脚本、输入元数据和小型结果。
