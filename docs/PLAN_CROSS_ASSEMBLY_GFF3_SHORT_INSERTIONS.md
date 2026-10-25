# 跨 assembly 短 insertion：GFF3、长读段连接与核心蛋白归并计划

日期：2026-10-24（America/Los_Angeles）。状态：GFF3 两株注释试验已完成；跨 FNA 生产流程和蛋白家族归并待实现。

## 1. 目标与本轮决策

从多个 assembly 联合发现候选，允许同一家族在每个样本中只有一个位置；
用配套 Nanopore 或 PacBio 长读段检查连接，并寻找更多上下文。
主要分析 **core 严格 <2500 bp、主要包含 1 或 2 个完整蛋白编码基因** 的候选。
core 长度不包含 upstream/downstream。没有已知转座酶注释也保留。

本文把“core protein”解释为 **insertion 内编码的同源核心蛋白**；宿主核心基因只用于亲缘/正交位置判断。
这两个概念分别存储。先产生 protein family，再区分 insertion 架构，最后汇总各样本位点。
蛋白同源、同一 insertion 亚型、同一正交位点和独立插入事件不是同一计数单位。

这是独立扩展实验。现有单株 PILER/RepeatScout 基线及 CLI 行为不变，
不能把本计划描述为已经实现 pooled FNA 或 GFF3 输入模式。
与 [既有跨物种计划](PLAN_CROSS_SPECIES_LONG_READS.md) 并存；本计划具体化短片段和蛋白归并路线。

## 2. 本地 GFF3 试验提供了什么证据

详见 [试验、脚本和结果](../experiments/gff3_short_insertion_pilot/README.md)。
H16/H850 的同版本 FNA/GFF3/FAA 下载和校验已完成，7 条 contig 序列与冻结 SHA256 一致。

| 筛选/检查 | 结果 |
|---|---:|
| 原 assembly 差异候选 | 70 |
| core <2500 bp | 44 |
| 严格注释子集合格 | 14，均为 1 个完整 CDS |
| 需复核 | 30 |
| 边界与注释特征相交 | 25；与其他原因可重叠 |
| 完整包含 2 个 CDS 的记录 | 4，全部有其它复核原因 |
| 合格 CDS 从 FNA 翻译后与 FAA 一致 | 14/14 |

这支持采用 GFF3 做结构注释，也显示不能仅按“长度合适、内部有一个蛋白”认定完整 insertion。
主要问题是候选边界和基因边界的关系。25 个相交记录不等于 25 个假阳性：真实插入也可能打断宿主基因。
下一阶段须区分宿主基因被打断、候选 core 被截短、重叠基因和注释异常，不能一律删除或自动扩边。
14 个合格记录也可能包含普通 accessory gene，仍需多个不同上下文和 reads 证据。

本次仅重注释已有坐标；未重新进行 read 审计或蛋白聚类。
旧发现下限为 400 bp，因此 44 不覆盖所有 <2500 bp 候选，后续需单列更短片段的灵敏度测试。
未将此前 400–20000 bp 的候选产量或两个 run 的前缀 reads 结果解释为本目标的准确率。

## 3. FNA 数据源：NCBI Assembly，FNA/GFF3/FAA 同批冻结

主来源采用 **NCBI Assembly，经 Datasets 批量获取**；优先有注释的 RefSeq，GenBank 补充缺失样本。
每个 accession.version 同时取得 genomic FNA、GFF3、protein FAA、assembly/sequence reports。
这几种文件在 [NCBI genome package 文档](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/how-tos/genomes/download-genome/) 中有明确下载方式。
本次小试验使用同一 NCBI assembly HTTPS 目录；规模化采用冻结 accession 清单的批量包及可恢复下载。

每条 assembly 记录保存：accession.version、BioSample、BioProject、taxon_id、strain、
assembly level、注释来源/日期、全部 contig 的 accession.version/长度/序列 SHA256、文件 URL/哈希和抓取时间。
GCA/GCF 配对只有核对序列对应关系后才去重；不要仅删除 `NZ_` 前缀推断等价。
相同 BioSample 的重组装和重复提交不能增加独立样本数；异版本保留映射及选择理由。
无 GFF3 的 assembly 仍可用于序列发现，标注 annotation_missing；可选统一重注释须独立版本化。

候选大 batch 延续 **至少 100 个去重样本** 的要求，E. coli 等充分研究类群降低发现优先级。
这是数据池要求，不要求每个距离子组都达到 100；也不能把 100 个 assembly 自动当成 100 个配套长读段样本。
分别报告 n_unique_samples、n_assemblies、n_annotated_samples、n_linked_ONT_samples、n_linked_HiFi_samples 和项目数。
研究充分程度使用可追溯的人工标签/文献依据，样本少不等于研究少。环境字段保留缺失和来源，不从物种名推断环境。

## 4. 跨 FNA call-out：距离决定比较范围，局部同源性决定证据

1. 将大 batch 按实测 assembly ANI/AF 建立邻居关系；用近邻加代表样本降低全对全开销。
   保留 query/ref 两侧 AF、计算版本和未返回距离的状态。用少量全对全子集检查近邻抽样的漏检。
2. 首轮复用现有 pairwise assembly gap caller，与单株发现结果并列保存。
   每条 contig 保持独立身份；允许多 FASTA 作为文件容器，不拼接不同 contig 的碱基。
3. 对每个 `U-X-D` / `U-D` 候选记录双侧比对、方向、core 坐标及竞争匹配。
   反向互补坐标、串联重复、复合 gap、重排、contig 末端和 circular origin 单列。
4. 不同 pair 对同一 sample/contig 的重复调用合并为 locus，并保留所有原始 pair 证据。
   必须比较 core 重叠、边界及侧翼定位；不把邻近两个独立插入合并，也不把不同 BioSample 的同源位点计作独立背景。
5. 导出逐 locus 的 U/core/D，再进入 GFF3 和长读段检查。只有 core 匹配而侧翼不可定位的记录保留为 unresolved。

`skani_min_af >=40` 或 asm20 query 覆盖超过一半目前只是小面板的**比较优先级线索**，不是验证过的有效性阈值。
距离扫描中的 45 对共享 10 个基因组，不能当作 45 个独立重复。no-hit 和低 AF 不直接永久排除，
留作局部可靠上下文的对照组；相反，高 AF 也不能保证某个位点的两侧正交。
PanGraph 可作为多样本表示/对照，当前 gap 调用结果不应称作 graph-bubble caller。
两份 assembly 的有无差异只支持 presence/absence，历史插入/删除方向未定。

## 5. GFF3：基因数、完整性和边界分别记录

根据 [NCBI GFF3 规范](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/file-formats/annotation-files/about-ncbi-gff3/)
处理 1-based 坐标、同 ID 的多段特征和 Parent 层级；内部统一使用 0-based 半开坐标。
不以 GFF 行数当基因数，不以相同 protein_id 合并不同位点。记录 CDS phase、strand、partial、pseudo、exception、
RNA 和缺失蛋白；GFF ID 不是跨注释版本的稳定标识，必须绑定注释文件哈希。

主分析子集先采用本次严格条件；同时保留三类复核输出：边界相交、基因完整性异常、其它结构/基因数。
基因端点用于提示复核，**不自动把 DNA insertion 边界改成 CDS 边界**。
发现阶段保留全部候选；目标子集只是后处理视图，不改变上游原始输出。
小规模阶段将所有合格 CDS 从 FNA 翻译后核对 FAA；扩展 importer 前增加 phase、多 Parent、
跨原点、缺注释、重叠基因和注释版本变化的挑战样例。

## 6. 如何 link Nanopore/PacBio 到 assembly

关系应是 `assembly ↔ BioSample ↔ experiment/run ↔ file`，而非一个 FNA 固定对应一个 FASTQ。
优先记录 assembly 明确引用的 source reads；否则以同 BioSample 建立候选关联，进一步核对培养物/菌株和实验类型。
同一个 BioProject 或相似 strain 名只能进入待查列表；同 BioSample 也不自动证明是组装所用数据或独立重复。

下载清单来自 [ENA file report](https://ena-docs.readthedocs.io/en/latest/retrieval/programmatic-access/file-reports.html)，
保留 run、experiment、BioSample、platform、library strategy/source、文件 URL、bytes、MD5、
原始提交类型和 read-unit 说明。没有 ENA FASTQ 时再检查 SRA 原始提交文件；镜像同一 run 不算第二份数据。

ONT、PacBio HiFi/CCS、CLR/subreads 分开。PacBio 类型不能仅凭文件名或 Sequel 型号判断，
须记录元数据/提交文件证据；多个 subreads 不自动等于多个独立分子。
现有 [Sam8H1 HiFi pilot](../experiments/pacbio_hifi_pilot/README.md) 已报告单样本兼容性通过，
但其 assembly 来自同批 reads，且读名被改写后 ZMW 不可恢复；不能推广为全平台独立验证。

每个候选先用对应 DNA 序列/家族内已观察的 DNA 成员检索 reads，再检查两端连接、
唯一可定位侧翼、CIGAR 连续性和竞争比对。核心内部命中只算 homolog hit。
回到该 run 所属 assembly 定位；新增位置只有符合相同证据标准后才加入 locus 表，
无法定位时保存 read-context 记录，不伪造 assembly 坐标。每轮搜索固定种子库版本，新增成员复核后才进入下一轮。
完整存档下载后核对校验和；reads 前缀只用于可行性检查，不能估计缺失/频率/总拷贝数。

## 7. 基于 core protein 的再次归并

采用三层模型：`protein_family → insertion_architecture → sample_locus`。

- **Protein family**：使用完整的核心 CDS 蛋白，记录双向覆盖、identity、蛋白长度和 domain 架构。
  首先核查保守聚类，再单独评估更远同源；宽松相似性得到的只是候选同源关系。
  [MMseqs2](https://github.com/soedinglab/MMseqs2/wiki) 可实现蛋白检索/聚类，覆盖模式和聚类模式必须冻结；
  不能通过 A–B、B–C 连通就默认为 A–C 属于同一种完整结构。
- **Insertion architecture**：一个或两个 CDS 的组合、顺序、相对方向、完整性、DNA core 边界和长度。
  共享一个蛋白但第二个基因不同的记录，先保留为同 protein family 下不同结构亚型。
  普通宿主同源基因不能仅因蛋白相似而成为 insertion family。
- **Sample locus**：记录该结构出现的样本、复制子、坐标、U/D 上下文和 reads 证据。
  同一正交位点在不同菌株的保留不算新的插入上下文；两个侧翼都须评估局部同源/共线关系。

保留一基因/两基因、完整/截短、核酸相似/蛋白远同源的分层结果。
无功能标签的蛋白参与结构比较；功能名不能替代同源分析。
相似性放宽的验收依据是边界和上下文仍可信、错误归并受控，不是“候选数变多”。
不会因某一蛋白功能标签提高其优先级；本阶段评估的是结构发现与归并质量。

## 8. 最小输出表及计数规则

| 表 | 关键字段/用途 |
|---|---|
| assemblies | assembly_id、BioSample、taxon_id、FNA/GFF/FAA 哈希、版本、QC |
| assembly_run_links | assembly_id、run、link_basis、source_evidence、read_type、assembly_consistency |
| pair_calls | pair、原始坐标、CIGAR/侧翼证据、caller 版本、失败状态 |
| loci | locus_id、sample/assembly/contig/start0/end0、core 哈希、边界状态、原始 call 列表 |
| locus_features | locus_id、注释快照、gene/CDS、protein_id、完整性、架构 |
| protein_membership | protein_instance_id、family_id、覆盖/identity、归并依据 |
| insertion_architectures | architecture_id、protein family 组合、方向、边界、成员 locus |
| read_support | locus/read_context_id、run、支持单元、双侧连接、竞争/冲突、独立性 |
| family_summary | protein family/architecture 分层统计样本数、正交上下文数、位点数、平台支持 |

家族、架构、BioSample、locus、read/molecule 分开计数。汇总表保留 excluded/review/no-hit/failed，
重复运行顺序不影响确定性 ID；新增数据不静默重写旧 family membership。

## 9. 执行顺序与验收

| 阶段 | 工作 | 可验收输出 |
|---|---|---|
| P0，已完成 | H16/H850 GFF3 审计 | 44 行、30 个 review、14 个严格合格；MD5/contig/GFF 对应、6 个挑战测试、14 个翻译检查通过 |
| P1，下一步 | 复核全部 44 个短候选，特别是 25 个边界相交；画 U/core/D+CDS 图 | 每个记录有边界解释和保留/待查原因；不自动修剪；选定不同背景及同正交位点对照 |
| P2 | 先完成两株跨 pair 去重和蛋白/架构归并，再扩到少量近缘已配套样本 | 单株一拷贝、跨株多背景可找回；同一位点重复 pair 不膨胀；domain-only/普通保守基因对照不错误归并 |
| P3 | 完整 ONT 数据审计；用已有 HiFi pilot 接口做独立平台兼容性评估 | 下载完整性、read type/link evidence、双侧支持/冲突、读段或分子单位完整报告；同源内部命中不过关 |
| P4 | 扩展到去重 ≥100 样本的数据池，按实测距离分子组 | 冻结 manifest、各阶段分母、缺失与失败原因、资源用量、每个家族可追溯到位点与 reads |

只有完成 P1–P3 后再冻结生产用相似度/距离参数。更远距离和更宽松匹配作为分开的敏感性分析，
保留基线结果；不要只用候选数增加作为改善证据。此计划未提交新的集群作业，也未改动现有 package。
