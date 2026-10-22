# 长读段可用性预检（冻结 pool 前）

检查日期：2026-10-22  
范围：候选八株的公开 reads 元数据；未下载 FASTQ、未运行发现或比对。  
关联：[候选表](CANDIDATE_ASSEMBLIES.md)、[实验计划](PLAN_CROSS_SPECIES_LONG_READS.md)。

## 结论与数据选择

- 原八株继续作为阶段 A 的 assembly 候选组成；本次未替换成员。
- 阶段 C 的配套 PacBio 优先顺序：DdD110 → RoMoAcy1 → AmHOT1A3。
  三者均找到 WGS/GENOMIC 的公开 run 和配套来源证据；下载及校验仍待执行。
- 原八株尚未确认匹配的 Nanopore run。
- 另找到环境土壤菌株 P. aromaticivorans AR20-38 的匹配 Nanopore 数据，
  登记为独立的阶段 B/C Nanopore 对照候选 PaAR2038，**不纳入 scale_2/4/8**，
  不将其 reads 关联到 BN5。它不是混合微生物样本，也不验证混合样本能力。
- pool 仅在输入版本、来源映射、注释可用性和 reads 元数据审查完成后冻结。
  目前 CP 组合输入和下载校验尚有待完成项，不宣称已形成可执行冻结 manifest。

## 八株逐项结果

| sample_id | 状态 | 已核查的证据 | 阶段 C 使用决定 |
|---|---|---|---|
| SaSB2B | no_matched_long_reads_identified | ENA 对 strain taxon 326297 的 PacBio/ONT 查询为 0 行；检索出的 SRX873015/16 是 Illumina RNA-seq | 仅 assembly 候选；不能称确定不存在长读段 |
| DdD110 | matched_metadata_download_unverified | SRR16235513；PACBIO_SMRT / Sequel；SAMN19030608；与 GCA_020510685.1 BioSample 一致 | 首个 PacBio 对照 |
| ShANA3 | no_matched_long_reads_identified | ENA strain taxon 94122 查询为 0 行；PRJNA13905 页面列 assembly 与 SAMN02598318，未列匹配长读段 | 仅 assembly 候选 |
| ScL1 | no_matched_long_reads_identified | ENA strain taxon 690566 查询为 0 行；名称搜索出现转座子库/RNA 数据，未确认匹配 WGS 长读段 | 仅 assembly 候选 |
| RoMoAcy1 | matched_metadata_download_unverified | SRR15616890；PACBIO_SMRT；SAMN12617337；原始论文将 run、该 BioSample 和 CP082158–CP082160 关联 | PacBio 对照；CP 输入版本与所属 assembly 解析仍待完成 |
| PaBN5 | no_matched_long_reads_identified | species taxon 2026199 的长读段查询只有 AR20-38 的 SRR11947722；BioSample 与 BN5 不同 | 仅 assembly 候选，禁止借用 AR20-38 reads |
| AmHOT1A3 | matched_metadata_download_unverified | SRR2084418；PACBIO_SMRT / PacBio RS II；SAMN03835502；与 GCA_001578515.1 BioSample 一致 | PacBio 对照 |
| AaH17 | no_matched_long_reads_identified | species taxon 589873 的查询为 0 行；原始组装论文使用 Illumina HiSeq 2000 | 仅 assembly 候选 |

“未找到”仅适用于以下检索范围，不是公共数据库中绝对不存在的证明。
尤其 strain taxon 查询可能漏掉登记在 species taxon、其他别名或未公开的记录。
这些行不能被解释为零连接支持；缺 reads 与比对后无支持是不同状态。

## 可定位的 runs 与下载元数据

以下均为 ENA filereport 返回的 archive 元数据。
read_count 是档案计数，PacBio subreads 不自动等于独立分子数。
MD5/bytes 是档案报告值，不是本地下载验证结果。

| sample_id | run | BioProject | BioSample | 平台/档案 instrument | reads | bases | FASTQ bytes | FASTQ MD5 |
|---|---|---|---|---|---:|---:|---:|---|
| DdD110 | SRR16235513 | PRJNA727717 | SAMN19030608 | PACBIO_SMRT / Sequel | 525543 | 2525540520 | 775843604 | 4b806b690b544dd0e89396a6bb51bc19 |
| RoMoAcy1 | SRR15616890 | PRJNA561397 | SAMN12617337 | PACBIO_SMRT / PacBio RS | 122352 | 1576526964 | 1357700967 | 15c63f8a6f101605744e5172d4bbdfff |
| AmHOT1A3 | SRR2084418 | PRJNA288688 | SAMN03835502 | PACBIO_SMRT / PacBio RS II | 163482 | 3217055373 | 2822177977 | e6984e8ad3613912959d7aea9b05a9a4 |
| PaAR2038（独立对照） | SRR11947722 | PRJNA624061 | SAMN14567048 | OXFORD_NANOPORE / PromethION | 132813 | 2112799366 | 2146409497 | 1cfd36ae0a44b7b05db73d8e0b28dee5 |

下载入口：

- [DdD110 FASTQ](https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR162/013/SRR16235513/SRR16235513.fastq.gz)
- [RoMoAcy1 subreads FASTQ](https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR156/090/SRR15616890/SRR15616890_subreads.fastq.gz)
- [AmHOT1A3 subreads FASTQ](https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR208/008/SRR2084418/SRR2084418_subreads.fastq.gz)
- [PaAR2038 FASTQ](https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR119/022/SRR11947722/SRR11947722_1.fastq.gz)

RoMoAcy1 论文写 PacBio RS II，而 ENA instrument 字段写 PacBio RS；
两者并列保存，不擅自改元数据。该 run 不是已核实的 HiFi/CCS 数据。
D110 的 Sequel 型号本身也不证明 CCS；按原始论文和文件类型核实 reads 类型后选择 preset。

## 别名和独立 Nanopore 对照

DdD110 的论文名称为 D110；NCBI assembly 的 strain 为 110，
BioSample sample name 为 D.den110。resolved manifest 必须保存
aliases=[D110, 110, D.den110]，并以 SAMN19030608 作为来源关联证据。
纯字符串 110 不能脱离 organism/BioSample 跨项目全局匹配。

PaAR2038：Paraburkholderia aromaticivorans AR20-38，
GCA_012689525.1，SAMN14567048，PRJNA624061。
NCBI Datasets 当前为 Complete Genome，sequencing_tech 包含
Illumina NovaSeq PE150 与 Nanopore PromethION；论文标题使用 Draft，
两种表述保留，正式下载时再核对复制子和组装状态。
来源为意大利 South Tyrol Alpine forest soil。
这份对照不能充当 BN5 的 read support，也不参加八株背景标签定义。

## 查询方法与范围

ENA Portal API：

- search endpoint：https://www.ebi.ac.uk/ena/portal/api/search
- result=read_run
- query=tax_tree(TAXID) AND (instrument_platform="PACBIO_SMRT" OR instrument_platform="OXFORD_NANOPORE")
- format=tsv，limit=1000；本次各查询均少于 1000 行，未触及返回上限。
- fields=run_accession,study_accession,sample_accession,secondary_sample_accession,scientific_name,sample_title,instrument_platform,instrument_model,library_strategy,library_source,fastq_ftp
- TAXID 与结果行数：326297→0；94122→0；690566→0；589873→0；
  2026199→1（AR20-38）；28108→8（其中只有 HOT1A3 行被采用）。
- 对 SRR16235513、SRR15616890、SRR2084418、SRR11947722
  分别调用 filereport，补取 read_count,base_count,fastq_md5,fastq_bytes。
- 上述四个 run 的 library_strategy=WGS、library_source=GENOMIC。

查询未穷尽所有 taxon/别名；species 查询返回其他菌株时没有自动采用。
NCBI SRA 若干网页读取失败，改用 ENA 的同一 SRR 档案记录；
不把网页失败计为没有数据。未进行 FASTQ 文件下载。

NCBI Datasets API 交叉核对：

https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/ACCESSION/dataset_report

已检查 GCA_020510685.1、GCA_001578515.1、GCA_012689525.1 的
organism、strain、assembly_info.biosample、BioProject 和 sequencing_tech。
RoMoAcy1 的来源关联依据下面原始论文，未声称解决了其 assembly accession。

补充原始来源：

- [D110 组装及 reads 论文](https://doi.org/10.1128/mra.01020-21)
- [D110 BioProject/assembly/BioSample](https://www.ncbi.nlm.nih.gov/bioproject/727717)
- [MoAcy1 原始论文](https://doi.org/10.1128/mra.00814-21)
- [HOT1A3 完整组装论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC4781885/)
- [AR20-38 原始论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC7330236/)
- [H17 原始论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC4119200/)
- [SB2B SRX873015](https://gold.jgi.doe.gov/sraexperiment?id=SRX873015)
  与 [SRX873016](https://gold.jgi.doe.gov/sraexperiment?id=SRX873016)：RNA_SEQ / TRANSCRIPTOMIC / Illumina HiSeq 2500。
