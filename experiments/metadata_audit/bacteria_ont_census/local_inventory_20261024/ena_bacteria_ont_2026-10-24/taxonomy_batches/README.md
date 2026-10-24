# 物种归并与逻辑 batch

2026-10-24 ENA 细菌 ONT 候选池；用 NCBI taxonomy 的 species / genus 层级归并。

原 8,605 个 taxon ID → 7,959 个 species 层级分组，另有 92 个无法归到 species 层级的保留组。

生成 1,018 个 genus 逻辑 batch，另有 252 个无 genus 的单独 batch。全部 93,682 个 runs 均保留且只归属一个组。

## 实际分批的组织方式

以 genus batch 为组织单元，多个 species / sp. 编号共享一个 batch。只有一个 run 或样本的 species group 不需要单独成批。物种归并表中的每行是分类单位，不是一个处理任务。

例如 Streptomyces sp. NBC_01197 归入 genus_1883（Streptomyces），与其他 Streptomyces 类群一起规划处理。组内每个样本的原始 run、reads 和未来 assembly 对应关系仍独立保留；同 batch 不代表拼接 reads 或合成一个生物样本。

数据密度仅作筛选标签，不是分批边界；同一 genus batch 可以同时包含 singleton、稀疏和高密度物种。较大的 genus batch 在决定下载范围后再按体积或样本量拆为执行子批。

## 数据密度标签（不是处理 batch）

以下标签按有 FASTQ 的去重 sample accession 数划分，阈值只是清单统计规则，不代表研究热度或新颖性。旧文件中的 review_batch_id / review_batches.json 是这些审阅标签的历史名称，处理归属应读取 batch_id / batch_catalog.jsonl。

| 审阅标签（旧 ID） | 每个 species 的 FASTQ 样本数 | Species 分组数 |
|---|---|---:|
| B01_very_dense | very_dense_1000_plus | 8 |
| B02_dense | dense_100_999 | 58 |
| B03_medium | medium_10_99 | 271 |
| B04_sparse | sparse_2_9 | 1,171 |
| B05_singleton | singleton_1 | 6,278 |
| B06_no_fastq | no_fastq_sample | 173 |
| B07_unresolved_species | 物种层级待核查 | 92 |

## 文件入口

- [可阅读的物种归并表](SPECIES_GROUPS.md)：名称、属、样本数、研究项目数、数据密度及原 taxon 数。
- [可阅读的 batch 概览](BATCHES.md)：每个属的物种组数、样本数和登记体积；点击 Batch 可查看成员。
- [全部 batch 的成员明细](BATCH_MEMBERS.md)：包含只有一个样本的 species group，统一按属组织。
- [七个数据密度审阅标签（旧文件名）](review_batches.json)：这些统计分层不决定处理 batch，也不是新颖性排名。
- [物种分组数据](species_groups.jsonl)；[batch 清单](batch_catalog.jsonl)；[逐 run 归属](run_membership.jsonl)；[分类映射证据](taxon_mapping.jsonl)。
- [汇总与快照来源](summary.json)。

## 怎样理解“没被挖透”

本轮只完成分类归并、数据密度和来源覆盖统计。数据少可能是研究少，也可能是上传少、缺乏组装、分类名称变化或测序平台偏好，不能直接当作新颖性。

后续应结合具体研究主题的文献、既有系统筛查、样本独立性、环境覆盖及同一 BioSample 的 assembly 可用性，再判断是否值得优先研究。没有为全表计算“研究饱和度”，没有据此选定样本。

用户于 2026-10-24 确认：降低 E. coli 等研究充分类群的优先级。E. coli（species taxon 562）已标为 low_reference（低优先级参考组），保留全部记录。其他类群的 review_priority 暂为 not_assessed，需结合文献逐项判断，不按数据量自动降级。不把整个 Escherichia 属或 Shigella 自动并入。

同一物种下菌株级 taxon 会合并；无 species 祖先的记录保留各自 taxon，不把所有“sp.”或“uncultured”误合成同一物种。

本次有 18 个原 taxon ID 在当前 NCBI 快照中已删除；这些记录保留原名称并标为待核查，没有凭名称强行匹配。

species-rank 中仍可能有未正式命名的 sp. 类群；species 层级并不保证正式双名或鉴定可靠。

Studies 数为不同 study accession 的数量，不意味着独立实验或独立样本。相同 sample 跨物种出现时，各物种内去重，但跨组统计不可直接当成全局唯一菌株数。

batch 是用于浏览和后续规划的属级集合，不是可以立即提交的下载/计算任务。超大属后续可以按所选样本量或体积拆分。原始文件和来源字段未修改。

[NCBI taxonomy 层级与名称说明](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/data-processing/taxonomy-processing/taxonomy/)
