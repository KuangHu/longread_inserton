# 按科合并的规划 batch

原 1,270 个属级/未定属 batch → 368 个科级 batch，另保留 150 个没有科级祖先的原分组，共 518 个规划 batch。

合并依据是 NCBI taxonomy 中相同的 family 祖先 taxon ID。Taxon ID 是标识符，数值大小或差值不能作为亲缘距离。

层级：科级规划 batch → 原属级 batch → 物种/编号类群 → BioSample → run。同一科里不同属可以共享一个规划 batch；原属级信息仍保留。

这是批次组织方式，不是把样本、reads 或 assembly 混成一个样本，也不表示科内所有样本都适合直接比较。下载和执行任务可在选择范围后按样本量/体积拆成子批。

E. coli 的低优先级只作用于该物种的成员，不传播到整个科。无明确 family 祖先的记录保持原 batch，不根据名称猜测。

## 文件

- [科级 batch 概览](FAMILY_BATCHES.md)
- [科级清单](family_catalog.jsonl)：科、所含属、原 batch、分类组、去重样本数及容量。
- [分类组到科级 batch 的映射](group_to_family.jsonl)
- [每个 run 的归属](run_membership.jsonl)：保留原属级 batch_id，新增 family_batch_id。
- [核对与示例](summary.json)

计数基于 2026-10-24 冻结的 ENA 清单与已保存的 NCBI taxonomy 快照。没有下载新的 reads，全部 93,682 runs 均恰好保留一次。

[NCBI Taxonomy 数据模型](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/data-processing/taxonomy-processing/taxonomy/)
