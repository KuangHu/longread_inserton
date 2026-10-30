# Nanopore organism 与环境来源目录

基于 2026-10-24 冻结的 ENA 细菌 ONT 候选池，逐 run 补充 ENA 索引中的样本元数据。

共 93,682 runs，8,605 个 taxon ID / scientific name 组合。这里不是物种数：名称可能停留在属级，或细化到菌株级。

## 环境信息覆盖

| 信息类别 | Runs |
|---|---:|
| 至少一个分离来源或结构化环境字段有有效值 | 63,859 |
| 仅有宿主/宿主部位上下文 | 11,844 |
| 上述结构化字段均缺失或为占位值 | 17,979 |

前三类互斥且覆盖全部 runs。字段缺失表示本次 ENA 索引未提供有效值，不代表原始 BioSample XML、论文或补充材料中也没有。

## 逐字段覆盖

| 字段 | 有效值 runs | 占位值 runs | 空白 runs |
|---|---:|---:|---:|
| scientific_name | 93,681 | 0 | 1 |
| strain | 57,043 | 2,066 | 34,573 |
| isolation_source | 57,603 | 4,293 | 31,786 |
| environment_biome | 9,269 | 664 | 83,749 |
| environment_feature | 6,469 | 1,550 | 85,663 |
| environment_material | 7,920 | 1,242 | 84,520 |
| broad_scale_environmental_context | 9,269 | 664 | 83,749 |
| local_environmental_context | 6,469 | 1,550 | 85,663 |
| environmental_medium | 7,920 | 1,242 | 84,520 |
| host | 52,133 | 2,768 | 38,781 |
| host_scientific_name | 51,624 | 21 | 42,037 |
| host_tax_id | 50,007 | 0 | 43,675 |
| host_body_site | 5 | 0 | 93,677 |
| country | 80,793 | 2,855 | 10,034 |
| location | 26,828 | 0 | 66,854 |
| collection_date | 75,864 | 12,857 | 4,961 |

## 如何查看每条记录

- [逐 run 完整注释](run_organism_environment.jsonl)：每行一个 run，保留 organism、taxon、strain、BioSample、环境原文、宿主、地点、采集日期和下载可用状态。
- [本地检索数据库](organism_environment.sqlite)：`run_context` 表，可按 organism、环境字段、是否有 FASTQ 等组合检索。
- [organism 汇总](organisms.jsonl)：每个 taxon/name 的 run 数、去重样本数及来源信息覆盖。
- [可直接阅读的 organism 全量表](ORGANISMS.md)；[环境与来源原文概览](SOURCES.md)。
- [来源原文汇总](reported_source_values.jsonl)：按字段和值汇总 run 数与去重样本数，未把相似词自动合并。
- [完整统计和对齐核对](summary.json)；[来源请求](provenance.json)。

## 解读原则

organism 是提交者标注的测序对象；host 是样本来源宿主，两者不混用。country/location 是地理来源，不直接等同于生态环境。

保留 environmental biome、feature、material 以及新版 broad/local context、medium 字段；不从物种名称猜测其环境。原始字段保留 not provided、missing 等占位字符串，但统计时将这些值视为缺失。

`context_summary` 只拼接有值字段及其原文，不做自动生态分类。sample title/description 保留供核查，未从自由文本猜测环境。

同一 sample 可以有多个 runs；来源汇总同时给出去重样本数。不同字段的统计可以重叠，不可相加为总样本数。

本轮没有查询原始样本 XML、做物种鉴定、判定宿主健康状态或筛选适合的样本；没有下载 reads。

## 索引对齐核对

- 冻结列表缺失的 context runs：0。
- 新增而未纳入冻结底表的 runs：0。
- sample / taxon / scientific name 变更记录：0。

[ENA 样本元数据返回字段定义](https://www.ebi.ac.uk/ena/portal/api/returnFields?result=read_run&format=json)
