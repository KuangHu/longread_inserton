# 细菌 Nanopore 全量目录与 batch 优先级

冻结日期：2026-10-24。此目录保存本地整理结果，是全量元数据候选池，不是 batch50 扩展。

## 阅读入口

- [候选池范围、下载状态与统计](ena_bacteria_ont_2026-10-24/README.md)
- [科级 batch 分组](ena_bacteria_ont_2026-10-24/family_batches/FAMILY_BATCHES.md)
- [batch 优先级规则](ena_bacteria_ont_2026-10-24/batch_priority/README.md)
- [全部优先级](ena_bacteria_ont_2026-10-24/batch_priority/PRIORITIES.md) 与 [首轮核查清单](ena_bacteria_ont_2026-10-24/batch_priority/SHORTLIST.md)
- [organism 清单](ena_bacteria_ont_2026-10-24/organism_environment/ORGANISMS.md) 与 [环境/来源字段概览](ena_bacteria_ont_2026-10-24/organism_environment/SOURCES.md)
- [其他数据来源调查](OTHER_DATA_SOURCES.md)

## 当前结果与口径

ENA 查询为 `tax_tree(2) AND instrument_platform="OXFORD_NANOPORE"`，含全部文库类型；查询前后及导出均为 93,682 个唯一 run，涉及 83,074 个 sample accession。89,134 个 run 列有 FASTQ 链接。此处只证明元数据中存在链接，未逐文件验证远端可下载性或读取质量。

按 NCBI taxonomy 祖先关系，得到 368 个科级 batch 和 150 个未能归科的保留分组，共 518 个。每个原始 run 的属、物种、样本与科归属均保留；batch 只是组织容器，不混合 reads 文件。

门槛为每个 batch 至少 **100 个去重 WGS FASTQ sample accession**，另保留至少 2 个来源项目的条件。E. coli 已降低优先级，不用于凑足门槛。40 个科级 batch 达标并完成 386 条文献查询；另外 1 个数据达标分组因分类未解决单列待核查。

P1/P2/P3 是合格科中通用文献关注度的相对分位，不是已证明的插入元件挖掘饱和程度。样本编号也不保证独立菌株；GENOMIC WGS 数量和最大项目占比另列。100 个 reads 样本尚不代表 100 对 reads + assembly。

与[上级目录原有 census](../README.md)的口径不同：原有 84,889 个 run 限制 GENOMIC；本目录先保留所有文库类型。原有 census 已做过粗粒度 BioSample–assembly 匹配，本次尚未把这些匹配逐一核实并接入 batch 优先级。

## 清单存储与本地使用

普通说明和主要汇总可直接阅读。较大的 JSONL/TXT 无损压缩；较长的逐 run 文件按 5,000 行分片。`compressed_manifest.json` 记录每片和恢复后文件的 SHA256、字节数及顺序。

在此目录运行：

```sh
python restore_snapshot.py
```

恢复过程不联网，会校验内容并拒绝覆盖不同的本地文件。恢复后可使用原始相对路径读取全部 JSONL，包括 `all_runs.jsonl` 中的 FASTQ/submitted 下载地址、MD5 与字节数，以及逐 run 的分类、环境和 batch 归属。保留 386 条本轮实际采用的文献查询原始响应。

没有上传测序 reads、SQLite 数据库、重复派生文件、大型原始 ENA 响应或 NCBI taxonomy archive。对应本地文件的 SHA256 见 `omitted_local_files.json`。这些校验值不等于文件已在远端备份。

## 方法与复现边界

所附五个 Python 脚本依次完成 ENA 导出、环境补充、物种/属归并、科级合并和优先级计算。仅使用 Python 标准库；请求还需要 curl。归并和历史名称查询需要 `taxonomy_cache/taxdump.tar.gz`，校验值见省略文件清单。

在保留原始 API 响应、taxonomy archive 和先前 run registry 的本地工作目录，可以重建相同统计；本仓库保存冻结的完整清单、分组和优先级证据，可直接恢复审阅。仅凭仓库无法保证从持续更新的公共 API 重新获得字节一致的历史响应。

`build_ena_census.py --fetch` 和 `enrich_ena_context.py --fetch` 会重新请求当时的在线数据；不可在冻结快照上当作无损恢复命令使用。刷新应在新工作目录进行，记录新日期并检查脚本中的固定路径和快照断言。`prioritize_batches.py --fetch-literature` 同样会读取或新增文献缓存。

本次验证包括所有压缩内容的往返 SHA256、93,682 个 run 的唯一性与分类归属覆盖、518 个 batch 完整性、100 样本门槛和 386 条查询证据的对应关系。没有重新运行下载 reads 或 assembly 分析。
