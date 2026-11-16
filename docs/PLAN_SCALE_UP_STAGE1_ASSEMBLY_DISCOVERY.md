# 跨 assembly 短插入发现的两阶段扩量计划

> **2026-10-27 主线纠偏：** 发现目标是“相似 core、不同上下文”，不要求空位点。Gap 搜索仅作辅助。
> 当前主线见 [纠偏计划](PLAN_CORE_CONTEXT_DISCOVERY.md) 和 [可运行的 core-context 原型](../experiments/core_context/README.md)。
> 下文已完成的 fna_pool Stage 1/2 是 **gap-derived 辅助目录与已知位点验证**，不是新主线完成。
> 新入口已通过 H16/H850 真实 smoke（旧 R33 六位置找回），但 GFF3 完整性、独立蛋白入口和精细上下文分组仍待接入。


执行对象：负责 `KuangHu/longread_inserton` 的 agent。

本轮任务是完成第一步：冻结目标 assembly 池，下载并核验全部入选 FNA 及配套注释，完成跨 FNA callout，交付可供后续长读段检索的候选目录。请实施并验收第一步，不要只返回另一份计划。第二步仅定义接口，本轮不新增 long-read 下载、比对或新位点搜索任务。

本计划调整现有 S0–S6 的执行顺序，不另建一套 caller。已有 Comamonadaceae 结果作为回归基线；已运行的长读段任务可继续，其状态不阻塞第一步，也不因本计划自动取消。

## 两步的目标与边界

| 步骤 | 工作 | 完成产物 | 本轮执行 |
|---|---|---|---|
| 第一步 assembly 发现 | 下载目标池 FNA、GFF3、FAA；按 assembly 距离选择比较；callout；注释复核；初步蛋白家族归并 | 带版本、来源、位点、上下游和证据状态的候选目录 | 是 |
| 第二步 long-read 扩展 | 逐 batch 获取匹配的 ONT/PacBio；验证已有位点；用 core 实例检索新增 read contexts；更新位点组 | 增量 reads 证据和新增候选位点 | 仅预留接口 |

主要目标保持为 **core 长度严格小于 2,500 bp，主要包含 1 或 2 个完整编码单元**。长度不含上下游。保留没有已知功能注释的候选，不要求同一个样本中必须有多个拷贝。

这里的 core protein 指插入片段编码的蛋白；宿主侧翼蛋白用于定位和上下文比较，两者分别记录。蛋白家族、编码架构、assembly 位点和 site group 分别计数，不把任何一个直接称为独立插入事件。

## 目标范围与冻结规则

“下载所有 FNA”指 **冻结目标池中所有入选的代表 assembly**，不指下载整个 NCBI、所有历史版本或所有已整理的 ENA reads。

默认第一波沿用现有优先级目录中 `P1_lower_attention_review` 的类群，按已有 `review_rank` 推进；E. coli 等已降优先级的类群维持原标签。先导出明确的 `target_batches.tsv`，记录来源文件、哈希、taxon ID 和入选理由。不要重新按候选蛋白功能挑选类群，也不要把文献关注度代理指标称为已证实的挖掘饱和度。

已有目录可从仓库中归档的 `experiments/metadata_audit/bacteria_ont_census/` 查找；原始文件名为 `batch_priorities.jsonl`。优先使用当前归档版本，记录其路径与哈希。

**旧优先级表的样本数来自 reads 元数据，不能直接充当 assembly 数。** 对目标 taxon 重新列出 NCBI assemblies，计算去重后的 assembly 样本数、项目数和注释覆盖。主研究 batch 维持至少 100 个去重样本的要求；低于门槛的组单列，并只在有分类学依据时合并到上级逻辑 batch，保留原 taxon。不要为凑数量拼接无关物种。

逻辑 batch、实际距离比较组和下载任务批次是三个层级：逻辑 batch 可以超过 100 个样本；距离子组不要求达到 100；下载任务可拆小。冻结后分波执行，第一波未完成时不要无记录地扩大范围。若没有可用的优先级归档，先完成并报告目标目录，不擅自将范围扩展为全部细菌。

## A 核对基线并解除 reads 入选限制

1. 检查当前工作分支、未提交修改、既有任务和数据目录，保存代码版本、工具版本及参数快照。复用已修复的 S3、S4，实现文件以实际 checkout 为准。
2. 复用现有 `build_pool.py`、`download_fna.py`、`build_s2.py`、`s3_pair_calls.py`、`s4_annotate.py` 和调度脚本，新增可配置的扩量入口。
3. **修改 assembly pool 的选择逻辑。** 当前 `build_pool.py` 的主要入选方式是有 usable run 的样本加其他 RefSeq Complete Genome；这仍会漏掉无配套 reads 的其他合格 assembly。第一步必须由 assembly 范围和质量条件决定入选，reads 关联只是元数据列。
4. 对同一批旧输入和同一组参数做回归，核对 caller、S4 分类和关键位点。新池加入更多伙伴后数量可以变化，不能强行要求扩量后的 main 清单仍等于旧数量。

交付：`baseline.json`、`target_batches.tsv`、入选规则与代码差异、基线回归报告。

## B 建立 assembly 清单与下载 manifest

主来源采用 NCBI Assembly，通过现有 Datasets 工作流获取版本固定的清单和文件。优先复用已核验的本地数据。RefSeq 与 GenBank 配对时保留映射，并在序列等价核验后去重。

- 每个目标 taxon 列出当前有效 assembly，冻结 `accession.version`，保留 BioSample、BioProject、taxon ID、strain、assembly level、来源和注释状态。
- 主池沿用 isolate genome 的范围；MAG、宏基因组、明显混合样本及不满足现有质量规则的数据单列，写清原因。
- 同一 BioSample 多个 assembly 采用可重复的代表选择规则，保留其他版本及选择理由。缺少 BioSample 的条目不自动计为新的独立样本。
- GCA/GCF 等价副本、完全相同序列及近克隆分别记录。近克隆保留样本对应关系，避免其重复提交虚增样本和候选数量。
- 有效染色体和质粒均保留。assembly level 较低不自动等于无用；按冻结的质量规则处理，contig edge 在 caller 中另行报告。
- 配套文件为 genomic FNA、GFF3、protein FAA、assembly/sequence reports。GFF3 或 FAA 缺失时，合格 FNA 仍可参与序列发现，但其候选标为 `annotation_pending`，不能列入已通过编码验收的 main。
- 重用已有 ONT/PacBio 关联信息，但不为凑齐 link 阻塞下载。关联区分明确 source run、exact BioSample、未解决和本次未找到；仅同 BioProject 不自动建立配对。

manifest 至少记录：batch、BioSample、assembly accession.version、文件角色、URL、可用的预期大小与存档校验值、本地路径、SHA256、状态、尝试次数、失败原因、来源快照。

下载队列采用 `pending / downloading / verified / failed`，支持从已有完整文件或可续传文件恢复。先计算文件数量、预期容量和临时空间需求；分任务下载，原子完成后再更新可用状态。已存在报告不等于已通过验收，必须读取其状态。

**下载完成门槛：**冻结清单中的每个 assembly 都有明确最终状态。损坏或缺失 FNA 不进入 caller；缺失注释与缺失 FNA 分开。若排除条目后形成实际执行子集，保存 manifest 新版本和差集，不把“部分下载完成”写成“全部完成”。

## C 核验 FNA 与注释的一致性

核对存档提供的校验值、本地 SHA256、文件可解析性和 contig accession/长度。没有提供的校验字段记为缺失，不能伪造为已校验。

FNA、GFF3、FAA 必须来自同一 assembly 和相容的注释版本。检查 GFF3 seqid 与 FNA 的映射、坐标范围以及 protein ID 与 FAA 的关联；保留注释来源和日期。不得只删掉 `NZ_` 前缀就认定 contig 或 GCA/GCF 相同。

输出每条 contig 的稳定标识和来源索引。联合 FASTA 中各 contig 独立成 record，名称可映射回 `sample / assembly / contig`；不串接不同 contig，不把多个样本变成一个样本。逐 record 的序列哈希应与原 FNA 一致。

交付：`assembly_catalog.tsv`、`download_manifest.tsv`、`assembly_aliases.tsv`、`contig_index.tsv`、`annotation_files.tsv`、文件校验报告和失败清单。

## D 按实测距离建立比较任务

1. 在逻辑 batch 内测量 ANI 与双向 alignment fraction，沿用现有 skani 工作流，保存工具版本和参数。
2. 复用当前的 clone、species-neighbour、nearest-other-species 与审计对照策略。以实测 pair 距离记录可比性；cluster 标签不代替具体 pair 的 ANI/AF。
3. 初始使用已验证配置，不将 pilot 阈值宣称为普适生物学边界。当前 species 分组的 ANI 95、双向 AF 50 等规则属于冻结配置，变化需新版本。
4. 保留近距离、跨项目和距离分层对照；很远的 pairs 单独汇总，不强行把全部目标 FNA 做全局 all-vs-all。
5. skani 无结果或没有入选 pair 的样本标为待审查或未比较，不能算作零候选。不同逻辑 batch 的可重复样本和 pair 全局去重。

交付：距离结果、clone/species 分组、`pairs.tsv`、选择理由、未比较样本清单及每 batch 的资源估计。pair ID 和输入哈希应稳定，已完成的相同任务可复用。

## E 运行跨 FNA callout 并保留原始证据

复用当前经过回归的 pairwise caller，以 U–core–D 与另一 assembly 对应位置的差异产生候选，保留两个比对方向及所有原始 calls。

- 广义发现仍沿用既有长度范围；目标目录另筛 `core_len < 2500`。同时报告 caller 的下限，目前为 50 bp，不能声称覆盖所有小于 2,500 bp 的差异。
- 双侧锚定、flank 竞争、重复拷贝、contig edge、one-flank、replacement 和 unresolved 均保留。
- `presence_absence`、`complex_small_other` 与 `replacement` 分开。保存另一侧片段的实际长度，不能让 `two_flank_anchored` 掩盖替换。
- 在同一 carrier 内合并对同一 locus 的重复 calls，保留每个来源 pair 和边界。不同 carrier 的位点此时不压缩成一个位点。
- 保存 `all_calls_anchored / partner_mixed / direction_mixed`、逐伙伴状态与边界一致性。一个方向成功或一个伙伴成功不等于所有 calls 一致。
- 边界争议只标记，不为了包含完整基因而移动 core 边界。

每个 pair 必须有成功、失败或明确排除状态。失败 pair 不产生“零候选”结果。重试失败任务，合并时检查任务和原始 call 标识，避免重复计数。

交付：原始 calls、合并 loci、one-flank 表、逐 pair 汇总、异常清单、资源记录以及实际 U/core/D 序列。

## F 复用 S4 做编码与边界复核

对所有目标长度的 anchored loci 运行 S4，不仅检查 main 候选。

- GFF3 CDS 按 ID 归并为编码单元，处理 strand、segments、phase、Parent、partial、pseudo、翻译例外和跨圆形原点特征。
- 每个单元绑定其 FNA/GFF3/FAA 来源；将基因组翻译结果与 FAA 比较，核对起止密码子、内部终止和边界关系。
- **扩展到其他 taxon 时，读取并校验注释的 `transl_table`，不能把现有 Comamonadaceae 的 table 11 无条件用于所有类群。** 未支持或相互矛盾的翻译规则标为待审查，不能默认为通过。
- `transl_except` 只作用于其指定位置和氨基酸；程序性移码与无例外支持的破坏分开记录。保留现有反例测试。
- 蛋白翻译一致性属于序列和注释检查，不写成表达或功能已验证。

维持互斥的主分类，并另保留所有可重叠的诊断标记：

| 分类 | 规则 |
|---|---|
| main | 纯 presence/absence；所有 calls 双侧锚定；1–2 个通过编码检查的单元；无边界切断；边界 exact；非 tandem |
| partner_mixed | 其他 main 条件满足，但伙伴比较不一致；附逐伙伴证据 |
| direction_mixed | 其他 main 条件满足，但比对方向不一致 |
| short_other_segment | 另一侧短片段单列，保留所有锚定一致性字段 |
| replacement | 替换候选单列 |
| review | 边界、注释完整性、重复结构、缺失注释及其他结构问题，分别记录具体原因 |

交付：逐 locus 的 S4 表、各候选清单、编码单元表、注释来源、异常案例及摘要。

## G 建立跨 batch 候选与初步蛋白目录

本轮完成到“可交给第二步使用”的 assembly 候选目录：

1. 对候选蛋白按完整序列哈希去重，保留 protein accession、编码单元和所有 locus 的对应关系；protein ID 数不等于不同序列数。
2. 复用当前保守家族参数：90% identity、双侧 80% coverage；50% 层仅作为候选同源分组。启用最终代表重新分配，并保存每个成员对最终代表的实测 identity、双向 coverage 和通过状态。
3. 各 batch 的局部结果汇入统一蛋白目录，在冻结的总目录上复核归并。保留局部到全局映射，不把每个 batch 的家族数直接相加。
4. 家族 ID 绑定目录版本和代表序列哈希。新增 batch 或改变代表时提供 ID 映射，不能让已有 locus 的证据失去关联。
5. 保存编码单元的顺序和方向，区分相同蛋白家族与相同双基因架构；保留实际 core DNA，不用单个蛋白代表或共识序列替代所有 DNA 实例。
6. 保存每个 locus 的 U/core/D、长度、坐标、方向、哈希和 flanking context。沿用已有各侧 300 bp 的上下文导出，并保留可从完整 assembly 扩展窗口的坐标；重复和近边缘不足长度要明确标记。

候选目录至少能回答：该 core 来自哪个 assembly、与哪些伙伴比较得到、属于哪个初步蛋白家族/架构、上下游是什么、有哪些争议、是否有可关联的长读段元数据。

第一步的证据来源标为 `assembly_only`。已有 pilot 的 read 证据可以通过稳定 ID 单独关联，但不能传播为整个家族或新 batch 已验证。

## 第二步预留的增量接口

本轮只交付接口，不执行新 reads 搜索：

- `assembly_run_links`：保留 run、平台、ONT/HiFi/CLR 类型、关联依据和文件元数据。
- `read_batch_manifest`：候选 locus、carrier/partner、关联 runs、已验证本地文件、未知或缺失状态、估计下载量。元数据整理不触发文件下载。
- `candidate_core_catalog`：core 实例、蛋白家族、架构、U/D 及来源哈希。
- `read_evidence`：后续按 locus/test/run 追加支持、冲突、定位歧义和计数单位，不覆盖原始 assembly 发现证据。
- `new_read_contexts`：预留新上下文 ID、core 匹配、两侧上下文、分子支持、assembly 定位及复核状态。

第二步有两个不同任务：当前 v0.3 验证已知位点；新模块用 core 实例检索 reads 并产生额外上下文。现有 v0.3 不能描述为已经实现新增 site 搜索。最终 S6 根据背景归并为 site groups，保留重复观察，不直接输出独立进化事件数。

## 目录和交付物

复用 `experiments/fna_pool/`，新扩量结果放独立版本目录，不覆盖 Comamonadaceae 基线。具体文件名可兼容已有 schema，但需给出对应关系。

| 交付物 | 内容 |
|---|---|
| run config 与目标清单 | batch 范围、纳入排除规则、版本、参数、资源估计 |
| assembly 与下载目录 | accession.version、样本去重、文件与 contig 哈希、注释来源、失败清单 |
| 距离与 pair 目录 | 实测 ANI/AF、比较理由、审计对照、未比较样本 |
| calls 与 loci | 全部原始来源、差异类型、边界与锚定一致性 |
| S4 清单 | main 及其他类别、编码检查、逐伙伴争议 |
| 蛋白家族与架构目录 | 序列去重、成员对最终代表验收、跨 batch 映射 |
| core 与侧翼目录 | 实际 DNA 实例、U/D、完整来源和可扩展窗口坐标 |
| 第二步 manifest | 元数据关联和增量接口，不下载 reads |
| 第一阶段报告 | 各 batch 及全局去重统计、失败与排除、验收结果、下一步入口 |

大型 FNA、FAA、原始比对和序列表放 scratch；GitHub 按既有项目方式保存代码、计划、manifest、哈希和适量摘要。大型文件索引同时记录可恢复路径、生成命令、输入哈希与软件版本，只有 SHA256 不足以恢复数据。

## 验收条件与停止位置

满足以下条件，才报告第一步完成：

1. 本轮 batch 和 accession 范围明确冻结；没有因 reads 缺失而排除本应入选的 assembly。
2. 冻结清单中的条目全部有状态；用于 callout 的 FNA 全部核验，缺失/不匹配的注释有明确去向。
3. 所有计划执行的 pair 均有成功结果或记录明确的排除；未完成的下载或比较不伪装为零命中。
4. 调用、合并、S4 和编码例外的相关回归通过；新增数据问题有复现案例。更换翻译表、坐标处理等逻辑时添加针对性测试。
5. 每个 locus 可追溯到版本固定的 assembly、contig、原始 calls、注释和 core 序列。
6. main 与伙伴/方向不一致、短另一侧片段、replacement、待审查类别分开；边界诊断可重叠但分类合计闭合。
7. 每个初步家族成员对最终代表达到冻结阈值；蛋白 ID、序列、家族、架构、locus、BioSample 分别统计。
8. 跨 batch 共享 assembly、pair、蛋白和位点映射已经去重；各局部计数与全局计数关系可解释。
9. 已交付第二步所需的 core 实例和元数据接口，但没有新增 long-read 下载或搜索任务。

最终汇报逐 batch 和全局的：assembly 数、去重样本数、项目数、注释覆盖、成功/失败下载、比较 pair 数、raw calls、loci、各候选类别、不同蛋白序列、家族和架构数量，以及既有 reads 元数据覆盖。物种计数使用明确的 taxon/分组定义，organism 字符串数只能称 organism label 数。

若冻结范围仍有未完成条目，汇报已完成部分与阻塞清单，不把第一步标为完成。达到验收后停在 `assembly_catalog_ready_for_long_reads`，交付结果与第二步入口；等用户启动第二步。

## 参考实现

- [现有端到端计划](https://github.com/KuangHu/longread_inserton/blob/docs/read-unit-audit/docs/PLAN_CROSS_ASSEMBLY_GFF3_SHORT_INSERTIONS.md)
- [当前 assembly pool 入选逻辑](https://github.com/KuangHu/longread_inserton/blob/docs/read-unit-audit/experiments/fna_pool/build_pool.py)
- [现有执行结果与修复记录](https://github.com/KuangHu/longread_inserton/blob/docs/read-unit-audit/experiments/fna_pool/README.md)

执行前固定实际 checkout commit；以上分支链接用于定位实现，不替代运行时版本记录。

