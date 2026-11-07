# ENA 细菌 Nanopore 全量候选目录

快照日期：2026-10-24。先保留整个可发现候选池，再从中选择样本；不是 batch50 的扩充表。

最新批次入口：[按科进一步合并的规划 batch](family_batches/README.md)。沿 NCBI family 祖先合并后，原 1,270 个属级/未定属分组整理为 368 个科级 batch 加 150 个待核查分组；原属、物种、样本和 run 归属均保留。

已按用户确认的至少 100 个 WGS FASTQ 样本门槛整理 [batch 初步优先级](batch_priority/README.md) 与 [P1 核查清单](batch_priority/SHORTLIST.md)。E. coli 不用于凑足门槛；40 个科级 batch 达标并完成文献关注度查询。文献关注度是审阅代理指标，不代表已证明的挖掘饱和度；assembly 配对仍待核查。

已补充 [organism 与环境来源目录](organism_environment/README.md)：逐 run 保留分类名称、分离来源、环境、宿主、地点与采集日期，并提供 organism 和来源原文汇总。字段覆盖与缺失值核对见该目录。

已补充 [物种归并与 batch 目录](taxonomy_batches/README.md)：依据 NCBI taxonomy 归并菌株级 taxon，按属组织 batch；多个单样本类群可在同一个 batch 中。七个数据密度标签仅供审阅，不是处理批次。数据密度不等于研究饱和度；原始候选均保留。

## 范围与完整性

范围为 ENA 公开索引中的 Bacteria 分类树（taxon 2）和 OXFORD_NANOPORE 平台：

`tax_tree(2) AND instrument_platform="OXFORD_NANOPORE"`

没有预先按物种、WGS、样本数据量、是否已有 assembly、原有 registry 或 batch50 进行筛选。
这不是所有生物的 ONT 目录：古菌、真核、病毒及归在细菌分类树之外的 microbiome/metagenome 记录不属于本次查询范围。
索引中的分类、平台和文库类型均来自提交者元数据，尚未通过序列验证。

使用 ENA `read_run` 搜索并设置 `limit=0`。查询前、查询后计数及实际导出记录数均为 **93,682**；run accession 无重复。
ENA 不提供本次查询的事务式快照；计数相符是完整性核对，不能证明查询期间每条记录都未变更。

## 当前候选池

| 类别 | Run 数 | 含义 |
|---|---:|---|
| 全部公开索引记录 | 93,682 | 保留所有文库类型 |
| 列有 FASTQ 链接 | 89,134 | 文件链接、字节数、MD5 元数据齐全 |
| 只有 submitted 文件链接 | 2,495 | 提交格式另有记录，不能默认等同于标准 FASTQ |
| 暂无上述文件链接 | 2,053 | 保留待查；不能据此断言 SRA 等渠道也不可获取 |
| 任一文件链接存在 | 91,629 | 前两种有链接类别之和；未逐文件联网验证 |

- 83,074 个非空 sample accession，6,347 个 study accession；不将 sample 数等同于独立菌株数。
- FASTQ 有 89,800 个不同文件 URL，登记的压缩体积合计 **56,626,324,129,579 bytes（约 56.63 TB，十进制）**。
- 文库类型为 WGS 的有 79,491 runs，其中 75,630 列有 FASTQ 链接；这里只提供统计，没有按此筛掉其他记录。
- 原有 registry 中的 65 个 ONT runs 均出现在全量候选池中。
- 下载地址、大小、MD5 的分号列表长度一致；本次未发现文件元数据完整性问题。MD5 仅检查格式，未下载 reads 验证内容。

## 使用哪些文件

| 文件 | 用途 |
|---|---|
| [all_run_accessions.txt](all_run_accessions.txt) | 全部 93,682 个 run ID，每行一个 |
| [runs_with_file_links.txt](runs_with_file_links.txt) | 91,629 个列有任一文件链接的 run ID |
| [fastq_run_accessions.txt](fastq_run_accessions.txt) | 89,134 个列有 FASTQ 链接的 run ID |
| [submitted_only_run_accessions.txt](submitted_only_run_accessions.txt) | 2,495 个只有提交文件链接的 run ID |
| [no_file_link_run_accessions.txt](no_file_link_run_accessions.txt) | 2,053 个暂无文件链接的 run ID |
| [run_inventory.sqlite](run_inventory.sqlite) | 本地可检索目录；`runs` 一行一个 run，`files` 一行一个文件 |
| [all_runs.jsonl](all_runs.jsonl) | 全量元数据，每行一个 JSON 记录 |
| [fastq_available_runs.jsonl](fastq_available_runs.jsonl) | 有 FASTQ 链接的元数据视图 |
| [file_inventory.jsonl](file_inventory.jsonl) | 文件 URL、大小、MD5 与所属 run；fastq/submitted 可能是同一数据的不同表示，不应全部重复下载 |
| [summary.json](summary.json) | 文库类型、来源、平台型号、容量及核对结果 |
| [provenance.json](provenance.json) | 确切 API 请求、字段及起止时间 |
| [ena_raw.json](ena_raw.json) | 未加工的 ENA 返回结果 |

主目录保存 run、experiment、sample/BioSample 及其别名、study、taxon、科学名、平台型号、文库类型/来源、read/base counts、公开/更新时间，以及 FASTQ 和提交文件的链接、大小、MD5。文件链接按原顺序与 MD5/大小逐一对应，多文件 run 不折叠成一个文件。

## 后续选择前还需要的标记

本轮完成的是 **reads 全量候选池**，并未完成全量 assembly 关联。当前 `assembly_link_status=not_assessed`，不表示没有 assembly。
后续应在保留此底表的基础上增加同一 BioSample 的 assembly accession、版本、组装级别及关联证据；同一物种或同一 BioProject 不能直接视为同一样本。
一个样本可能有多个 runs 或多个 assemblies，应保留对应关系，再决定采用哪些记录。

文件可获取状态目前指 ENA 元数据列有链接；本次没有对数万个远程文件逐一探测，也没有下载 FASTQ 或 assembly。
本地是否已有完整文件均未检查。全量元数据目录与之前的局部传输测试是不同的核对范围。

## 来源

- [ENA 高级检索及全量导出说明](https://ena-docs.readthedocs.io/en/latest/retrieval/programmatic-access/advanced-search.html)
- [ENA 文件报告说明](https://ena-docs.readthedocs.io/en/latest/retrieval/programmatic-access/file-reports.html)

`build_ena_census.py --fetch` 会更新同一目录中的快照；需要保留历史版本时，应先另建新的日期目录再运行。没有 `--fetch` 时，仅使用已保存的原始元数据重建本地索引。
