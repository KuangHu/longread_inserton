# 全 Bacteria ONT genomic runs 普查（metadata only，雏形）

日期：2026-10-24。脚本：[ont_census.py](ont_census.py)。**只做数量普查：不下载 reads，不下载 FASTA，不冻结任何下一批。**

## 数据来源

- ENA read_run：`tax_tree(2) AND instrument_platform="OXFORD_NANOPORE" AND library_source="GENOMIC"`，共 84,889 runs。
- NCBI `genomes/genbank/bacteria/assembly_summary.txt`：流式扫描 3,420,265 行，按 BioSample 精确匹配，命中 38,467 行，只用 `version_status=latest`。
- 原始文件留在集群 `/global/scratch/users/kh36969/longread_callout/bacteria_ont_census/`；URL 与 ENA 表的 SHA256 记录在 [census_summary.json](census_summary.json)。

## 漏斗（按 run 计数）

| 步骤 | runs |
|---|---:|
| genomic ONT runs | 84,889 |
| WGS | 75,659 |
| 非 PCR | 72,077 |
| ENA 有 FASTQ MD5 和 bytes | 68,289 |
| 精确 BioSample 对上 latest GenBank assembly | 38,897（35,627 个 BioSample；872 个属） |
| Complete Genome 或 scaffold ≤10 | 30,246（其中 Complete 21,510） |
| 粗略来源筛：ENA host/isolation_source 无临床/动物关键词 | 16,225（742 个属） |

按属的计数见 [genus_counts.tsv](genus_counts.tsv)。

## 局限（这是上限估计，不是资格表）

- **最长记录未检查**：没有 FASTA，draft 的“≥3 Mb”条件没有判。
- **来源筛很粗**：只看 ENA run 的 host/isolation_source。很多 run 的这些字段为空，会被算作“通过”；batch50 用的是 BioSample attributes，这里没有逐个查询。例如 Escherichia 有 3,774 runs“通过”，大多应是缺元数据，而不是环境来源。
- **病原物种未排除**：batch50 只排除了 B. pseudomallei/mallei，这里没有物种级排除。
- **没有克隆合并**：没有 ANI，同一 BioSample 的多个 run 也未合并（run 数大于 BioSample 数）。
- 属名取 ENA scientific_name 的第一个词，所以 “Candidatus” 会被单独算作一个“属”。

## 下一步若要用它选批次

需对候选属逐个查询 BioSample attributes 和 FASTA（与 [batch50](../../callout_scale_v1/batch50/README.md) 同样流程），
再做 skani 克隆合并。本目录不作这一步。
