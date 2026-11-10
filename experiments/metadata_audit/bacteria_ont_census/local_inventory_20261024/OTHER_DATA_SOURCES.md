# ENA Nanopore 之外的数据来源

调查日期：2026-10-24。以下区分新的存储入口、新的测序平台与新的样本来源；本轮是官方资料调查，未建立新的全量 census，也未实测这些站点的 bulk 下载。

| 来源 | 对本项目的用途 | 增量与限制 |
|---|---|---|
| NCBI SRA / DDBJ DRA | 交叉检索 ONT/PacBio 元数据，补查 ENA 缺失链接，提供备用访问方式 | 与 ENA 同属 INSDC，许多记录共享同一 accession；不能把三站数量相加 |
| PacBio HiFi/CCS 数据 | 如果扩展到其他长读长平台，可在 ENA/SRA 另建目录；官方示例数据适合方法验证 | 是平台扩展，不一定是新数据库；PACBIO_SMRT 不自动等于 HiFi，需区分 CCS/HiFi、旧 CLR 和未知类型 |
| CNCB-NGDC GSA 原生提交 | 搜索公开非人类 ONT/PacBio reads，寻找现有池之外的项目 | GSA 同时整合 INSDC 记录；需区分原生提交和镜像，增量需去重后计算 |
| NCBI Datasets / GenBank / RefSeq | 为现有 reads 补齐 assembly 及 BioSample、项目、质量元数据 | 是 assembly 来源，不是新增原始 reads；GenBank/RefSeq 对应组装也需避免重复 |
| JGI Genome Portal / IMG / GOLD | 环境微生物项目、组装、注释与样本环境信息的发现入口 | JGI raw reads 也发布到 SRA；数据入口可能要求登录或项目条件，不能把全部记录视为新的长读长样本 |
| ONT/EPI2ME 与 PacBio 官方公开数据集 | 可追溯的示例与验证数据 | 精选集合适合验证，不能代替全量样本目录 |

## 依据

1. [DDBJ DRA 官方说明](https://www.ddbj.nig.ac.jp/dra/index-e.html)明确其与 NCBI SRA、ENA 的 INSDC 协作关系。[NCBI SRA 云元数据字段](https://www.ncbi.nlm.nih.gov/sra/docs/sra-cloud-based-metadata-table/)包含 SRR/ERR/DRR、BioSample、平台和项目等。备用入口可能有同步、格式和访问差异，仍需逐记录核对。
2. [PacBio 官方数据集](https://www.pacb.com/connect/datasets/)及[公开下载目录](https://downloads.pacbcloud.com/public/dataset/)提供不同应用的数据；[DDBJ 文件格式说明](https://www.ddbj.nig.ac.jp/dra/datafile-e.html)列出 PacBio 的 subreads/CCS 等文件类型。应保存平台和 read 类型证据；没有足够证据的记录保留为未知，不凭平台字段判定 HiFi。
3. [GSA 官方浏览入口](https://ngdc.cncb.ac.cn/gsa/browse/run/list)说明整合 INSDC 元数据；[NGDC 2026 数据资源论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC12807729/)区分其自有资源与整合的 INSDC 数据。因此 GSA 原生提交是值得核查的补充，不能直接把网站总量作为新增量。此次浏览入口一度返回 502；未验证稳定导出 API 或批量下载成功率。
4. [NCBI Datasets 官方下载说明](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/how-tos/genomes/download-genome/)支持按 taxon、assembly accession、BioProject 获取组装数据包，也支持仅取 metadata。这可用于优先 batch 的配对审计。
5. [JGI 官方下载说明](https://portal-web-2.jgi.doe.gov/portal/help/download.jsf)列出 RAW、assembly、QC 数据与 API/Globus 方法，并说明 raw reads 也发布至 SRA。其环境项目发现价值不能直接等同于新的独立 reads 数量。
6. [EPI2ME data releases](https://epi2me.nanoporetech.com/category/data-releases/)提供 ONT 公开发布数据；作为方法验证入口与全量档案检索分开记录。

## 下一轮目录应如何衔接

保持原先 Nanopore + assembly 主线：先利用已有 census 的 assembly 匹配结果，核查每个 batch 有多少可追溯的同 BioSample 配对。归档层面优先调查 GSA 原生 ONT 提交，SRA/DRA 用于交叉核对与备用入口。若决定扩展测序平台，再单独建立 PacBio HiFi/CCS 候选目录，保留与 ONT 的平台差异。

统一保留 `archive_origin`、原始 run/sample/project accession、BioSample 交叉引用、tax_id、平台、read 类型、文库来源、环境原文、文件 URL/格式/大小/校验值、assembly accession 与配对证据。按共享 accession 和明确交叉引用去重；不能因同名菌株或同一物种就认定同一样本。跨平台的同一 BioSample 在样本数中只算一次，多个 run 仍分别保留。

新增源暂不报“新增多少样本”或“多少 batch 达到 100 对”：这些数字必须等目录获取、去重和 assembly 配对后再给出。
