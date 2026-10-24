# 50 株 ONT 面板：来源/数据资格审计（metadata only）

日期：2026-10-24，LBL Lawrencium 登录节点。状态：**面板已按资格规则冻结，可进入下载/校验；尚未下载 reads，未运行任何 call-out。**
前置：五株验收（`nearby_competition_v2`）通过；[10 株批次](../batch10/README.md) 10/10 通过。
路线：培养分离株 assembly＋同一 BioSample 的 genomic Nanopore reads（[protocol](../../../docs/CALLOUT_SCALE_PROTOCOL.md) §2）。
PILER/RepeatScout 默认、边界推断、审核分类规则与报告定义均未改动。
**选择只用来源与数据资格；没有使用任何 candidate、tnp、read-spanned 数或 PILER/RepeatScout 输出。**

## 结论

| 问题 | 答案 |
|---|---|
| 扫描的 ENA ONT genomic runs | Burkholderiaceae（tax_tree 119060）全部 ONT runs 1,064 条，其中 library_source=GENOMIC **867 条**（817 BioSamples） |
| 精确 BioSample assembly 匹配 | **490 run×assembly 行**（490 runs，468 份 current GenBank GCA）；377 runs 无匹配 assembly |
| 通过连续性 | 去掉非 WGS（103）、B. pseudomallei/mallei（125）、ENA FASTQ 缺 bytes/MD5（1）后 261 行，连续性淘汰 20，**241 行通过** |
| 通过来源/文库 | 去掉已运行 15 株、临床/动物来源 58、ChIP 文库 1、mixed culture 1、PCR 文库 3（进 pcr_reserve）、同 BioSample 备用 run 4，**159 个新样本合格** |
| 克隆合并后的独立组 | 159 个合格样本 = **116 个 clone group**（ANI ≥99.99% 且双向 AF ≥95%）、65 个 ANI ≥98% 簇、49 个 ANI ≥95% 组 |
| 是否真有 50 个新样本 | **有**。每个 clone group 最多 2 个（已运行的成员计入上限）时，规则下最多可选 134 个 |
| 冻结状态 | **已冻结，可下载**：50 行都有 ENA bytes＋MD5；GCA/GCF 序列集核对已完成（49 份逐条一致，SFA1 没有 GCF，用 GCA）。剩余只是作业准备（见下），不是选择阻碍 |

选中的 50 株：
- 来自 **50 个不同 clone group、50 个不同 ANI ≥98% 簇**（39 个 ANI ≥95% 组，28 个 NCBI 物种名），与已运行的 15 株不共享任何 clone group 或 ANI ≥95% 组；
- 31 个 Complete Genome，19 个高连续性 draft；全部 library_selection = RANDOM；
- 来源：环境 25、富集培养 11、植物 14，无 unclear；
- 属：Burkholderia 18、Cupriavidus 16、Ralstonia 8、Caballeronia 3，Pandoraea、Imbroritus、Polynucleobacter、Limnobacter、Trinickia 各 1；
- ENA FASTQ 合计 48.5 GB（单株 0.09–6.0 GB），read bases 合计 53.0 Gb；read bases / assembly bp 为 17.2–1124.5，中位 117。

**Paraburkholderia 没有新的合格样本。** 属内 24 个 genomic ONT runs：14 个已在五株验收或 batch10 中运行，
5 个无匹配 assembly，5 个 assembly 连续性不达标（P. hospita ×2、P. xenovorans、P. tropica、P. fungorum 各一，11–81 条记录）。
因此按规则全部 50 株由 Burkholderiaceae 补足（`genus_tier` 列一律为 `2_burkholderiaceae`）。

## 规则（[build_panel.py](build_panel.py) 固定实现）

1. **扫描**：ENA read_run `tax_tree(119060) AND instrument_platform="OXFORD_NANOPORE"`，再取 GENOMIC。
   NCBI Datasets v2 对该科逐页扫描 16,844 条 genome report（17 页，全部翻完），只保留 GenBank current GCA，
   用 `assembly_info.biosample.accession == ENA sample_accession` 精确连接。
   NCBI 的 BioSample 端点对 SAMEA/SAMD 返回空，所以不用它；同一 BioProject 或相近株名都不算匹配。
2. **硬条件**（依次判定，`exclusion_reasons.tsv` 记录第一条不符合的原因）：
   WGS → 不是 B. pseudomallei/mallei → ENA 有单个 FASTQ 及其 bytes 和 MD5 →
   连续性（Complete Genome，或**下载的 GenBank FASTA**记录 ≤10 且最长 ≥3 Mb）→ 不是已运行的 run/BioSample/assembly →
   来源不是临床/人/动物宿主 → library_selection 属 RANDOM/other/unspecified/size fractionation → 纯培养 → 非 PCR。
   含 PCR 的文库进 `pcr_reserve`，不进主面板。
3. **来源分类**：对 BioSample 的 host、isolation_source、env_medium 用固定关键词规则（ENA 字段作补充）。
   人、昆虫、鱼、羊宿主，痰、血、脓、伤口、CF、医院等记为 clinical/animal 并排除；
   `Pathogen.cl` package 且无任何来源信息时也保守记为 clinical/animal。
   Rhizopus 内共生菌、工业/药液/ISS 废水、无元数据记为 unclear（可用，排序靠后）。
4. **每个 BioSample 一个 run**：非 PCR、RANDOM 优先，然后取 bases 较多者；其余进 reserve。
5. **相似度**：skani 0.4.0 `triangle --slow`，输入为合格池、PCR reserve 与已运行的 15 份 assembly（共 177 份）。
   clone group 用单链接，要求 ANI ≥99.99% 且双向 AF ≥95%；另标 ANI ≥98% 簇与 ANI ≥95% 组（仅作参考）。
6. **排序（发现之前）**：Paraburkholderia 优先 → RANDOM → 环境/植物/富集优先于 unclear → Complete 优先于 draft
   → read bases / assembly bp 降序 → run accession。
7. **多样性**：第 1 轮每个 ANI ≥98% 簇取排名最高者（已运行株所在的簇视为已有代表），第 2 轮每个 clone group 取一个，第 3 轮补每组第二个；
   每组上限 2，已运行成员计入。本次第 1 轮就选满 50，没有用到第二个 clone 成员。
   与已选或已运行株同 clone group 的其余合格样本进 `clone_reserve.tsv`，不删除。

## 注意事项（不改变资格，执行和解读时要带着）

- **来源依赖**：10 株来自 PRJEB88388（R- 系列，污染场地富集培养分离株；其中 5 株来自同一 HCH 富集样品 DE_BIT_08，2 株来自 IT_COL_04），
  6 株来自 PRJEB108492（LMG/R- 保藏的 Burkholderia cepacia complex 株）。它们彼此 ANI 不同，但不是相互独立的野外重复。
- **实验室参考株**：H16、H850、DSM 11853 及 LMG/R- 保藏株的元数据只记原始分离来源，实验室传代历史未核查。
- **基因组小**：Polynucleobacter TUM22923 是 1 条 1.86 Mb 的 Complete Genome，因 Complete 而合格（≥3 Mb 只适用于 draft）。
- **深度不均**：Q53 和 CT11-110 为 19.5× 和 17.2×（与 LjRoot223 的压力样本相近）；GP95 为 1124×（6.0 GB FASTQ）。不设深度硬门槛。
- **无 GCF**：SFA1 没有 RefSeq 配对，按 `prepare_assembly.py` 规则使用 GCA；GCA 无 CDS 注释时 tnp 状态为 unavailable。
- **被排除但可能有用**：96 个 R. pseudosolanacearum ONT runs 的 library_strategy 是 `Synthetic-Long-Read`，按 WGS 规则排除；
  Pericarpium Citri Reticulatae（陈皮）来源按关键词归为植物，但实际是药材产品。
- 来源规则是关键词，不是人工逐条判读；具体依据写在每行的 `source_note`。

## 选中的 50 株（[selection_manifest.tsv](selection_manifest.tsv)）

按排名排列；完整字段（BioProject、GCA/GCF、N50、仪器、文库、ENA counts/bytes/MD5/FTP、host、来源、地点、日期、clone/簇、核对状态）见 TSV。

| # | sample_id | organism | 等级 | 记录数 / 最长 Mb | bases/asm bp | 来源 | run |
|---:|---|---|---|---|---:|---|---|
| 1 | H16 | Cupriavidus necator H16 | Complete | 3 / 4.05 | 409.9 | environmental | SRR32057873 |
| 2 | DZH-419 | Ralstonia sp. DZH-419 | Complete | 3 / 3.70 | 317.3 | environmental | SRR39705685 |
| 3 | Tfa17 | Cupriavidus oxalaticus | Complete | 2 / 4.02 | 230.8 | environmental | SRR32057868 |
| 4 | DSM_11853 | Cupriavidus basilensis | Complete | 8 / 4.57 | 203.0 | enrichment | SRR12739613 |
| 5 | LMG_24269 | Burkholderia arboris | Complete | 5 / 3.68 | 186.0 | plant | ERR17405070 |
| 6 | TMH.W2 | Cupriavidus sp. TMH.W2 | Complete | 5 / 3.83 | 176.6 | environmental | SRR35788632 |
| 7 | PCR-C_Y9 | Burkholderia contaminans | Complete | 5 / 3.61 | 170.4 | plant | SRR31277720 |
| 8 | R-92558 | Pandoraea norimbergensis | Complete | 1 / 6.27 | 169.3 | enrichment | ERR15686594 |
| 9 | H850 | Cupriavidus necator H850 | Complete | 4 / 3.84 | 166.6 | environmental | SRR32057871 |
| 10 | SBC2 | Caballeronia sp. SBC2 | Complete | 8 / 3.99 | 165.7 | environmental | SRR11193104 |
| 11 | DG-3 | Ralstonia pseudosolanacearum | Complete | 2 / 3.94 | 158.9 | plant | SRR25270028 |
| 12 | KBS0801 | Burkholderia sp. KBS0801 | Complete | 3 / 3.53 | 117.2 | environmental | SRR10064718 |
| 13 | LMG_32992 | Imbroritus primus | Complete | 2 / 3.25 | 116.7 | environmental | ERR13929838 |
| 14 | LMG_34232 | Burkholderia sp. LMG 34232 | Complete | 3 / 3.60 | 111.9 | environmental | ERR17405047 |
| 15 | R-18634 | Burkholderia crassaminum | Complete | 3 / 3.69 | 106.8 | environmental | ERR13696569 |
| 16 | Mo-D5 | Ralstonia solanacearum | Complete | 3 / 3.11 | 104.7 | plant | SRR39418895 |
| 17 | TUM22923 | Polynucleobacter sp. TUM22923 | Complete | 1 / 1.86 | 99.5 | environmental | DRR439344 |
| 18 | R-1925 | Burkholderia arboris | Complete | 3 / 3.70 | 89.0 | environmental | ERR17405080 |
| 19 | R-50709 | Burkholderia sp. LMG 34232 | Complete | 3 / 3.67 | 76.0 | plant | ERR17405870 |
| 20 | LMG_24066 | Burkholderia arboris | Complete | 5 / 3.63 | 75.3 | environmental | ERR17405067 |
| 21 | HG_4 | Burkholderia gladioli | Complete | 4 / 4.16 | 74.1 | plant | SRR32316772 |
| 22 | R-15933 | Burkholderia mycopellens | Complete | 3 / 3.92 | 65.1 | environmental | ERR13696568 |
| 23 | R-20536 | Burkholderia sp. LMG 33096 | Complete | 3 / 3.77 | 64.5 | environmental | ERR17366079 |
| 24 | SFA1 | Burkholderia sp. SFA1 | Complete | 8 / 3.13 | 59.1 | environmental | DRR198995 |
| 25 | CS-K2 | Limnobacter thiooxidans | Complete | 1 / 3.54 | 56.7 | environmental | SRR25907068 |
| 26 | JAB02 | Burkholderia lata | Complete | 4 / 3.73 | 56.6 | plant | SRR39158507 |
| 27 | GBBC_1122 | Ralstonia pseudosolanacearum | Complete | 3 / 3.56 | 51.7 | plant | SRR32059352 |
| 28 | KL41 | Burkholderia sp. KL41 | Complete | 3 / 3.77 | 51.4 | environmental | SRR23378375 |
| 29 | GBBC_1956 | Ralstonia pseudosolanacearum | Complete | 2 / 3.62 | 30.3 | environmental | SRR32059353 |
| 30 | Q53 | Burkholderia ambifaria | Complete | 3 / 3.50 | 19.5 | plant | SRR17038419 |
| 31 | CT11-110 | Ralstonia sp. Ralssp110 | Complete | 3 / 3.67 | 17.2 | plant | SRR32164188 |
| 32 | GP95 | Ralstonia sp. GP95 | draft | 6 / 3.53 | 1124.5 | environmental | SRR31333963 |
| 33 | R-91474 | Cupriavidus basilensis | draft | 4 / 4.32 | 258.2 | enrichment | ERR15693033 |
| 34 | R-93091 | Cupriavidus campinensis | draft | 3 / 3.54 | 235.7 | enrichment | ERR15695518 |
| 35 | R-91747 | Cupriavidus basilensis | draft | 5 / 4.11 | 228.4 | enrichment | ERR15681449 |
| 36 | R-93106 | Cupriavidus pauculus | draft | 4 / 3.86 | 218.6 | enrichment | ERR15695470 |
| 37 | R-92745 | Cupriavidus necator | draft | 4 / 4.54 | 199.8 | enrichment | ERR15695471 |
| 38 | LAS2 | Burkholderia sp. LAS2 | draft | 5 / 4.46 | 189.7 | environmental | SRR15135747 |
| 39 | R-92680 | Cupriavidus nantongensis | draft | 4 / 3.61 | 184.3 | enrichment | ERR15686826 |
| 40 | R-92709 | Cupriavidus campinensis | draft | 3 / 3.50 | 181.9 | enrichment | ERR15686846 |
| 41 | R-92684 | Cupriavidus campinensis | draft | 2 / 3.62 | 162.4 | enrichment | ERR15741483 |
| 42 | FW507-4G11 | Cupriavidus sp. FW507-4G11 | draft | 2 / 4.52 | 156.4 | environmental | SRR37523914 |
| 43 | EB106-05-01-XG147 | Cupriavidus sp. EB106-05-01-XG147 | draft | 3 / 3.88 | 155.9 | environmental | SRR37523861 |
| 44 | TBRC_18363 | Burkholderia stagnalis | draft | 4 / 3.64 | 140.0 | environmental | SRR27732414 |
| 45 | R-92694 | Cupriavidus necator | draft | 6 / 4.31 | 105.8 | enrichment | ERR15686829 |
| 46 | FW03 | Caballeronia sordidicola | draft | 10 / 3.40 | 77.1 | environmental | SRR40321973 |
| 47 | HG_7 | Burkholderia glumae | draft | 5 / 3.59 | 50.8 | plant | SRR32316769 |
| 48 | LMG_27706 | Ralstonia syzygii subsp. celebesensis | draft | 3 / 3.62 | 45.4 | plant | SRR32059343 |
| 49 | LjRoot230 | Trinickia sp. LjRoot230 | draft | 2 / 4.22 | 27.6 | plant | ERR13194303 |
| 50 | LjRoot31 | Caballeronia sp. LjRoot31 | draft | 4 / 4.13 | 22.5 | plant | ERR13194358 |

## 文件

| 文件 | 内容 |
|---|---|
| [selection_manifest.tsv](selection_manifest.tsv) | 主面板 50 行（`order` 为执行顺序） |
| [ranked_candidates.tsv](ranked_candidates.tsv) | 全部 159 个合格新样本的排名及状态（selected / clone_reserve 28 / ranked_reserve 81） |
| [clone_reserve.tsv](clone_reserve.tsv) | 与已选或已运行株同一 clone group 的合格样本 |
| [exclusion_reasons.tsv](exclusion_reasons.tsv) | 其余每个 genomic ONT run（及其 assembly）的排除或 reserve 原因（700 excluded，5 reserve，3 pcr_reserve） |
| [clone_groups.tsv](clone_groups.tsv) | skani 输入的 177 份 assembly 的 clone group、ANI ≥98% 簇、ANI ≥95% 组 |
| [selection_summary.json](selection_summary.json) | 逐步漏斗计数、合格数与独立组数、选中分布 |
| [query_provenance.json](query_provenance.json) | ENA/NCBI 查询 URL、翻页计数、原始文件 SHA256、skani 命令/版本/二进制 hash、全部规则、脚本 hash |
| [panel_check.json](panel_check.json) | [check_panel.py](check_panel.py) 的结果：**pass** |
| `raw/*.gz` | ENA 返回表、匹配的 NCBI report、skani 边表（解压后的 SHA256 与 provenance 一致） |

下载的 342 份 GCA 和 68 份 GCF FASTA（2.6 GB）留在集群 `/global/scratch/users/kh36969/longread_callout/batch50_meta/genomes/`，
每份旁边有 `<accession>.json`，记录各条序列的长度与 SHA256。

## 校验

`check_panel.py` 不依赖 build 的中间状态，重新读取 ENA/NCBI 原始缓存和 FASTA 统计，检查以下各项：
- 选中行的 run、BioSample、assembly、sample_id 都唯一；
- 不含五株验收或 batch10 的任何 run、BioSample、assembly；
- 每行都是 ONT WGS/GENOMIC、非 PCR，run 与 assembly 的 BioSample 精确相同，有 ENA bytes＋MD5，不是 pseudomallei/mallei；
- 每行满足连续性，且 FASTA 统计与 manifest 一致；没有临床/动物来源；
- 每个被选用的 clone group（连同已运行成员）不超过 2 个，并且第二个成员标为 `second_allowed`；
- 所有 TSV 表头与固定列一致，行键不重复。

人为制造重复行、混入已运行的 run、改动 FASTA 统计后，检查都会报 fail。

## 复现

```bash
python fetch_metadata.py query                         # ENA + NCBI 全科扫描 -> batch50_meta/raw
python build_panel.py                                  # 缺 FASTA 时写 need_gca.txt 并退出
python fetch_metadata.py download --accessions .../need_gca.txt
python build_panel.py                                  # skani + 选择；写 need_gcf.txt
python fetch_metadata.py download --accessions .../need_gcf.txt
python build_panel.py && python check_panel.py
```

ENA 与 NCBI 都是活数据库，重查可能多出新 run 或新版本；冻结版本以 `raw/*.gz` 和 provenance 中的 hash 为准。

## 下一步（尚未做）

- 把 50 个 run 加入 `scripts/long_read_runs.json`，生成 batch50 的 `samples.tsv`，沿用 batch10 作业（`prepare_assembly.py` 再核对一次 GCA/GCF）。
- 存储：reads 48.5 GB，按 batch10 比例估计输出约 13 GB；ENA 下载按 0.3–1 MB/s 估计，合计 13–45 小时。
- **本审计没有提交 50 株 SLURM 作业。**

## 可执行化与 2 株 smoke test（2026-10-24）

选择未改动。[make_samples.py](make_samples.py) 由 selection_manifest 生成 [samples.tsv](samples.tsv)，并把 50 个 run（role `callout_scale_v1_batch50`）加入 `scripts/long_read_runs.json`。
samples.tsv 保留 batch10 的列（order…preset），另加 panel rank、ENA reads/bases/bytes/MD5/URL、选用的 FASTA 来源与 accession、期望 SHA256 和记录数、read policy。
所选 FASTA：49 株用 GCF（序列集与 GCA 逐条一致），SFA1 用 GCA（`prepare_assembly.py --gcf none`）。

[batch_sample.sbatch](batch_sample.sbatch) 的 array 序号就是 samples.tsv 的 `order`，每次提交 10 个，并发 2，例如 `sbatch --array=1-10%2 --exclude=<故障节点> batch_sample.sbatch`。
与 batch10 相比多两步：`assembly_expected` 要求集群上选出的 FASTA 与冻结的 SHA256 和来源一致；
ENA bases / assembly bp >200× 的 9 株（H16、DZH-419、Tfa17、DSM_11853、GP95、R-91474、R-93091、R-91747、R-93106）
在完整下载并通过 MD5 校验后，用 [downsample_reads.py](downsample_reads.py) 按分子降采样到 80×，只用于 minimap2 和 read 检查。
降采样规则：按 SHA256(“seed:read 名”) 排序取到目标碱基数，seed 固定为 20261024；输出的 gzip 头固定，SHA256 可复现。
原始 FASTQ 和降采样 FASTQ 留在集群，只记录 SHA256 与 seed（`results/<株>/qc/downsample.json`）。发现步骤不使用 reads。

Smoke test（只提交两株，**尚未提交其余 48 株**）：

| 株 | 角色 | 结果 |
|---|---|---|
| CS-K2（#25） | Complete、SAMN、小 FASTQ（0.19 GB，56.7×，全部 reads） | 全部阶段 exit 0；validate、checker 均 pass；11 个原生家族（0 个 tnp-supported），100 records（≥400 bp 30 条），100 条全部 audited_full_span |
| R-93091（#34） | draft 3 条记录、SAMEA、236× → 降采样到 80.0×（378,511 → 127,772 reads） | 第 2 次全部阶段 exit 0；validate、checker 均 pass；48 个家族（4 个 tnp-supported），275 records（≥400 bp 167 条）；audited_full_span 264，edge_limited 5，ambiguous_competing 6 |

- 两株在集群上选出的 FASTA 都是 GCF，SHA256 与 samples.tsv 中的期望值一致；reads 的 bytes、MD5、spots、bases 校验全部通过；均无 read 达到 secondary 上限。
- R-93091 第 1 次在 n0230 上失败：curl 加载共享库时报 Lustre `Error 108`，8 次重试都失败，未写任何结果。这与 DD10 在 n0258 上的节点故障相同，失败目录保存在集群 `batch50_failed/`。
  第 2 次加 `--exclude n0230` 重交后完成。故障节点名单现为 n0348、n0039、n0292、n0082、n0031、n0258、n0083、n0230。
- 结果见 [results/](results/)，结构与 batch10 相同（`batch_summary.tsv`、`resources.tsv`、`omitted_large_files.tsv`、`failed_attempts/`）；records 不是独立 insertion 数。
