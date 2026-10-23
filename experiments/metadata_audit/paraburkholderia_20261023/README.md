# Paraburkholderia genomic Nanopore 可用性与 ANI 核查

核查日期：2026-10-23。状态：元数据与 assembly 相似度核查已执行；面板未冻结。
关联[属级跨株计划 v2](../../../docs/PLAN_PARABURKHOLDERIA_MULTI_ISOLATE.md)。
本目录没有新增菌株的 insertion discovery、宿主位置比对或 read 支持结果。

## 结论

ENA 检出 **24 个 genomic Nanopore runs，24 个 BioSamples**。
逐一按 BioSample 精确关联 NCBI current GCA，得到 **19 份 assembly**：
3 个 Complete Genome（包括已有 AR20-38）、1 个 Chromosome、15 个 Contig。
另 5 个 BioSamples 在本次 NCBI/ENA 检索中未找到关联 assembly。
这不是对未关联、误标或其它别名记录不存在的证明。

P. aromaticivorans 只有既有 SRR11947722；原定新增同种 Nanopore 菌株的目标没有数据支持。
P. fungorum 虽有 9 runs，但 7 个匹配 contig-level assembly，2 个未找到关联 assembly。
这两个未匹配的 BioSamples 都标 NBRC 102489，不可当作两个独立菌株。

已实际下载 **19 份 genomic FNA（包含各自所有记录）**，验证 ZIP CRC 和总碱基数，
记录文件 SHA256，并运行 skani 0.4.0 对 AR20-38 的一对多 ANI。
没有下载新增 raw reads，也没有运行 PILER、RepeatScout 或跨株背景归组。

| Complete Genome 候选 | Nanopore run | ANI 到 AR20-38 | AF：候选 / AR20-38 | 来源限制 |
|---|---|---:|---|---|
| AR20-38，GCA_012689525.1 | SRR11947722 | 100.00% | 100.00 / 100.00% | 既有森林土壤基准 |
| DD10，GCA_049177125.1 | SRR33434047 | 86.54% | 54.08 / 38.34% | 草甸土壤菌丝际；BioSample 为 MISAG single-amplified 模型，培养分离与扩增历史待核实 |
| CT11-68，GCA_047325075.1 | SRR32164161 | 84.88% | 52.87 / 33.15% | NCSU Phytotron 小麦根部；植物生长设施来源，不能直接称野外分离株 |

两株新增候选的 ANI 都远低于原侧翼 90% 固定阈值，但 ANI 不等于每个位点的局部 identity。
AF 表明可比同源部分有限，更不能把未比对区段自动计为新增背景。
DD10 的 NCBI 分类仍为 Paraburkholderia sp. DD10；元数据中另有 P. terricola 描述，
此处保留档案分类，不自行改名。
CT11-68 的设施来源也不能直接证明其为工程化实验室菌株，需与野外来源分层。
**当前不能宣布已具备三个合格的“野外环境培养分离株＋Complete Genome”样本。**

## 全部 24 runs 与选择顺序

[availability.tsv](availability.tsv) 每个 run 一行，包含 BioSample、菌株、assembly level、
来源原文、ANI/双方 AF、FNA SHA256、候选状态及缺口。
[ena_runs.tsv](ena_runs.tsv) 保留 ENA 查询返回的全部原始字段与 FASTQ URL。
所有记录仍是候选核查状态，不能作为冻结执行清单。

在相同来源和 assembly 合格层内，按到 AR20-38 的 ANI 降序、accession 打破平局。
Complete Genome 暂定顺序为 AR20-38 → DD10 → CT11-68；后两株需先解决来源问题。
最接近 AR20-38 的有匹配 assembly 株为 P. xenovorans R-92571 / R-92578，
ANI 分别 88.64 / 88.49%，但它们是 contig-level，不能静默替换主面板要求。
同样，SMT4a 的 Chromosome 等级不等于 Complete Genome。
环境富集培养、同一取样点和同名保藏菌株都要保留来源依赖关系，不能把每个 run 当独立环境重复。

## 核查方法与范围

1. ENA Portal 的 read_run 查询使用 Paraburkholderia taxid **1822464**：
   `tax_tree(1822464) AND instrument_platform="OXFORD_NANOPORE" AND library_source="GENOMIC"`。
   返回 24 行，均为 WGS；没有仅按物种名搜索，也没有扩展到整个 Burkholderiaceae。
2. NCBI Datasets v2 属级报告两页共 1,631 条（1000 + 631，包含 GCA/GCF）；
   跟随 next_page_token 至末页，按 current GCA 去除 RefSeq 配对重复。
   精确连接 ENA sample_accession 与 assembly_info.biosample.accession，不能用共同 BioProject 替代。
3. 对 5 个未匹配 BioSamples 另查 NCBI BioSample genome endpoint 及 ENA assembly：
   NCBI 返回 NO_ASSEMBLIES_FOR_BIOSAMPLE_IDS，ENA 无记录。
   详见 [missing_assembly_checks.json](missing_assembly_checks.json)。
   这些样本分别是 CT11-26、KKMWBUW1、NK-14 和两个 NBRC 102489 BioSamples。
4. 19 份匹配 GCA 的元数据原样保存在 [matched_assemblies.json](matched_assemblies.json)。
   用 NCBI genomic FASTA ZIP 下载全部 assembly 记录，核查 CRC、文件成员及总碱基数。
   [assembly_download_audit.json](assembly_download_audit.json) 保存 URL、长度、记录数和 SHA256。
   这些校验不等于重新证明组装正确或原始测序来源无误。
5. skani 使用 AR20-38 为 query，其余及自身为 reference；完整原始结果见
   [ani_to_ar2038.tsv](ani_to_ar2038.tsv)，无手工四舍五入后重算。

```bash
./skani dist -q genomes/GCA_012689525.1.fna \
  --rl genome_list.txt --slow -s 80 --min-af 15 -t 2 --detailed \
  -o ani_to_ar2038.tsv
```

19/19 输入在输出中有对应行，退出码为 0；AR20-38 自比为 100/100/100。
ANI 是同源部分相似度估计，两个 AF 是不同的分母，不能只报告 ANI。
skani 官方定位主要为 ANI 高于约 80% 的基因组；WTPI1 的 79.88% 位于边缘，
这类低值不应用作精细排序或物种判定。`-s 80` 是筛选参数，不能保证最终估计值必然 >=80。
本次未运行所有候选之间的 all-pairs ANI，也未运行位置/共线性分析。

## 复核与下一步

- [provenance.json](provenance.json) 保存完整查询 URL、分页记录、软件版本和二进制 SHA256；
  [ani_provenance.json](ani_provenance.json) 保存实际 ANI 命令。
- [download_assemblies.py](download_assemblies.py) 是本次执行的 assembly 下载脚本，
  读取同目录 matched_assemblies.json；仅在新的空目录中复现，避免覆盖当前校验记录。
  下载完成后先检查 19 条状态全部为 downloaded_verified_length，再运行 ANI。
  软件 release 的 latest URL 会变化，复现必须核对记录的版本和哈希，不能盲用届时的 latest。
- FNA 和二进制不上传仓库；完整属级两页 JSON 的哈希已记录，相关 19 条元数据快照已上传。
  数据库内容会变化；这些结果是上述日期和查询范围的快照。
- 冻结前核实 DD10 的 single-amplified 模型与实验历史、CT11-68 的来源层级，
  再检查 run 文库、原始命名和处理历史。精确 BioSample 匹配本身不能替代文库审计。
- 后续计划以属级比较为主；共享位置依赖宿主共线性/基因顺序，未定位归 unresolved。
  每株每家族报告新增背景、可比较分母和未解决数，不使用“+1 即成功”。

## 原始来源

- [ENA Portal API](https://www.ebi.ac.uk/ena/portal/api/)
- [NCBI Datasets API](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/api/rest-api/)
- [skani 官方代码与使用范围](https://github.com/bluenote-1577/skani)
- [DD10 BioSample](https://www.ncbi.nlm.nih.gov/biosample/SAMN47339132)
- [CT11-68 BioSample](https://www.ncbi.nlm.nih.gov/biosample/SAMN43145758)

技术核对已通过：24 run 行唯一；19 个精确 BioSample/GCA 连接；19 个下载长度校验通过；
19 个 ANI reference 与下载 accession 一致；5 个缺失 assembly 检查；Complete Genome 数为 3。
这是元数据/ANI 审计的验证，不能算 package 单元测试或生物学路线验证。
