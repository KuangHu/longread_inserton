# 近缘 Paraburkholderia：逐株发现与 Nanopore 核查，跨株汇总家族

日期：2026-10-23  
状态：待执行计划 v1；本次只提交计划，未选择新增 run、未下载新增数据、未运行跨株分析。  
基准结果：commit [629d37d](https://github.com/KuangHu/longread_inserton/tree/629d37da0d738b6db4bb9e981ddff95be2ebe202)，
见 [READ_CALLOUT_RESULTS.md](READ_CALLOUT_RESULTS.md)。

## 1. 本轮目标和优先顺序

沿用已跑通的路线：

1. 每株完整 assembly 分别运行 PILER-DF 与维护版 RepeatScout。
2. 每株 Nanopore reads 回贴自己的完整 assembly，核查候选连接。
3. 跨工具、跨菌株比较候选 X，建立家族与位点对应表。
4. 输出一个 X 家族对应的全部实际拷贝及 upstream/downstream，
   分开统计样本数、位点数、不同侧翼背景数和 read 支持。

主要问题是：增加近缘独立菌株后，是否能增加同一 X 家族的可解析背景？
本轮不要求 perfect empty site、tnp、TIR 或 TSD 才保留候选。
没有 tnp 支持的家族照常保留；同源空位点只提供额外证据。

本轮主路线不是 pooled DNA discovery，也不是将不同菌株 reads 混在一起。
后者会混淆来源；A 株 reads 不能补作 B 株位点的支持。
旧 reads → SV → empty/filled 主路线不恢复。
低频非参考插入（包括此前仅两条支持的 IS5 线索）不作为本轮主要验收目标。

原[联合 FNA 计划](PLAN_CROSS_SPECIES_LONG_READS.md)保留为后续对照；
本计划优先执行，不依赖 RepLong 安装或其结果。
原八株 scale_2/4/8 的组成和冻结结果不改写。

## 2. 起点与尚未证明的部分

已完成基准：PaAR2038，Paraburkholderia aromaticivorans AR20-38，
GCA_012689525.1，BioSample SAMN14567048，Nanopore SRR11947722。
GCA/GCF 序列集合在此前实验中已逐条校验；复用时仍保存实际文件哈希。
这是环境土壤分离株，不代表未经培养或传代。

PaAR2038_RS_R1 在基准中有 11 个不同候选坐标，长度 1,171–1,173 bp，
11/11 candidate 完整包含已注释 tnp CDS；55 对侧翼比较未检出对应 flank HSP。
每个成员有 145–181 个 v1 full-span 计数单位。
这些是首轮操作性证据，不是已认证的精确元件末端或独立转座事件。
没有显著 flank HSP 也不等于侧翼绝无同源性。

基准的 762 条候选记录包含跨工具/跨家族重复，不能作为独立 insertion 数。
此处 11 拷贝家族用于回归与人工审计，不用于选择“必须带该家族”的新样本。
跨株比较后的家族数、不同背景数和增量目前均未知。

## 3. P0：先核实新增样本，再冻结面板

目标为 AR20-38 加 2–4 株独立菌株，总计 3–5 株。
这只是目标规模，不表示已找到足够匹配数据。
首轮仅检索 P. aromaticivorans；若数量不足，先报告缺口。
同属近缘物种可另建 near_species 扩展层，冻结前列明名单和原因，
不把它们冒充同种重复。现有 PaBN5（GCA_002278075.1）可进入元数据核查，
但目前没有已确认的配套 Nanopore，不能借用 AR20-38 reads。

| 条件 | 纳入主面板的要求 |
|---|---|
| 样本来源 | 环境来源有 BioSample/原始论文依据，非已知工程化实验室衍生株；来源不明单列 |
| assembly | 完整、可下载，染色体与所有关联质粒一并保留；实际记录清单核对通过 |
| reads | 同一 BioSample 的公开 Nanopore genomic WGS；仍检查文库/样本描述，不能仅凭 BioProject 相同 |
| 独立来源 | 同一菌株重测、同一 BioSample 的多个 run、GCA/GCF 配对不增加菌株数 |
| 版本 | assembly/复制子 accession.version、注释、reads 文件均锁定并保存校验值 |
| 注释 | 优先使用可核对序列的统一注释来源；缺注释不等于无 tnp，不因缺注释删除原生发现 |
| 读段可解释性 | 记录 chemistry/basecaller、处理历史、原始 read ID、split/duplex 关系的可得性 |

先取得 assembly 元数据及 FNA，再依据全基因组相似度和比对覆盖报告亲缘关系。
锁定 ANI 工具版本、命令、双方比对覆盖与 taxonomy 来源；ANI 仅用于分层，
不凭高 ANI 删除可能带有不同插入的独立菌株。
在符合来源/数据要求的候选中按到 AR20-38 的 ANI 降序选取，平局按
assembly accession 排序；不依据 tnp 数、IS5 有无或候选产量选样。
保留全部候选及未入选原因，便于评估选择偏差。

至少一个新增独立菌株匹配成功，才进入双株技术检查；
有三个合格独立菌株才称为完成目标中的最小多株 pilot。
不足三株时报告有限样本结果，不静默放宽样本类型或拿重复 run 凑数。

冻结清单至少包含：
sample_id、strain_aliases、taxon/来源、BioSample、assembly/replicon accession.version、
GCA/GCF 对应、sequence/file SHA256、注释版本、run、下载 URL、档案 MD5/bytes、
测序/处理元数据、匹配证据 URL、纳入状态和原因。
面板冻结发生在新增菌株 discovery 结果查看之前。
原三株回归 isolate 不加入本面板。

## 4. P1：逐株保留原生发现和候选结构

- 使用当前 package 的独立 sample 语义，运行 --subset all。
  两软件均保留所有原生家族、成员、日志和失败；不新增发现过滤器。
- 冻结 package 与二进制版本/提交、完整命令、输入和注释哈希。
  RepeatScout 两阶段实际 l 和 FASTA 文件字节数逐株记录；
  这里没有 pool 引起的 l 变化，若不同单株默认 l 不同也必须标出。
- 固定共同的蛋白证据来源，或把参考库差异明确列为不可直接比较的因素。
  tnp-supported 原定义不变，不凭跨株同源标签把无证据的成员升级为 tnp 阳性。
- 每条记录保存 native、independent inferred、candidate 坐标，
  boundary_source、实际 insertion 序列、左右 500 bp flank 与截断状态。
  负链按家族方向导出，同时保留基因组正向坐标，内部统一 0-based 半开区间。
- 完整 tnp CDS 计数使用 candidate_overlapping_tnp_CDS 中的 fully_contained，
  不用 native 区间的首条注释替代。
- contig/origin 受限、边界未解决、无效坐标分别保留。
  不跨 circular origin 人工绕回，不把未解析边界当作完整元件。

报告全部长度范围；另列 <400 bp 与 >=400 bp 层。
400 bp 在这里是汇总分层，不是删除短候选的规则，也不是证明完整 IS 的标准。

## 5. P2：跨株家族对应与株内重复报告归并

此部分是新增的比较/报告层，尚未实现；不能描述成 PILER/RepeatScout 原生输出。
原生 family/member ID 永久保留，新增 global_family_id、sample_locus_id 和 alias 表。
tnp 产品名、单个催化域或同属一个已知 IS 类别，不足以合并 DNA 家族。

### 5.1 家族对应

输入为所有可导出的实际 candidate X 序列，不带 flank；
同时保留 native/consensus 的对应关系用于诊断边界不一致。
使用 BLAST+ blastn（-task blastn -dust no，双链）保存逐对 HSP；
版本及其余参数写入冻结配置，不只保留最佳一条 hit。

本 pilot 预先选择的严格对应标准为：
非重复计数的共线、同向 HSP 覆盖双方各 >=90%，链上 identity >=95%。
保存双方覆盖、匹配/错配/gap 的分母和方向；禁止跨重复 HSP 累加出假覆盖。
这是操作性家族定义，不是已优化的物种或 IS 分类阈值。
部分匹配、嵌套序列和边界明显不一致的成员建立 related/partial 链接，
不把一个短片段作为桥梁合并两个大元件。

满足对应标准的图连通分量，仅当所有成员两两均满足标准时才自动形成严格组。
存在 A–B、B–C 命中但 A–C 不满足的分量标 clustering_ambiguous，
保留全部链接和原生分组，不用单链接传递强行合并，也不删除记录。
如需更宽阈值，只作为预先另存的敏感性分析，不覆盖严格组。

组代表使用真实成员序列：优先边界非截断者，其次长度最接近组内中位数，
再按稳定 member ID 排序。不拼造新的“天然共识插入”。
实际长度变体和各成员序列全部保留。

### 5.2 同一菌株内的 locus 对应

先折叠相同 sample/replicon/start/end 的重复坐标，保留所有来源别名。
同一严格 X 组内，方向一致且 candidate 双向区间覆盖均 >=90% 的株内记录，
可建立同位点对应；同样要求组内两两满足，不通过重叠链把相邻位点连起来。
只重叠一部分、嵌套或边界矛盾的记录标 locus_overlap_ambiguous，单列未解决数。
相邻串联拷贝不因相邻或同属家族就合并。
保存 native 和 candidate 的两套重叠证据，不悄悄裁剪或重写原边界。

同一位点的不同软件候选可有多个候选边界版本；每个版本的 read 证据仍独立保存。
位点汇总时合并支持 read 的身份集合，不能相加两个软件的支持条数。
若没有逐 read 身份表，暂不计算合并支持数。

## 6. P3：近缘菌株的同源位置与不同侧翼背景

以下四种量必须分开：

| 量 | 定义 |
|---|---|
| raw_member_count | 原生及导出成员记录数，允许重复报告 |
| sample_locus_count | 株内已解析且去重的位点实例数；同源位点在不同菌株各有一个实例 |
| resolved_context_group_count（N_context） | 在可评估的完整侧翼中，跨株同源位置归组后得到的背景组数 |
| read_supported_context_group_count | 至少一个位点达到 P4 支持标准的背景组数；逐株支持仍分别报告 |

N_context 是本次序列/共线性规则下的背景数量，不是历史插入事件数。
菌株数、物种数、未解决位点数、侧翼不足数同时输出。
例如三株各 11/8/6 个拷贝得到 25 个实例，不代表 25 个不同背景。

每个 locus 用真实左右 500 bp flank 分别对其它菌株的完整 assembly 比较，
保存所有竞争位置。初始近缘匹配标准：每侧覆盖 >=80%、identity >=90%；
两侧映射必须方向和次序一致，并共同指向可比较的宿主区域。
遇候选边界偏移、重复侧翼或多重匹配，辅以独立提取的左右各 2 kb 上下文
和全基因组共线关系；不以单侧匹配或 tnp 名称判定同源位点。
具体比对工具版本、参数、等分命中和共线性规则在新增数据解释前冻结。

只有双侧对应与位置关系无歧义时合并为同源背景组。
整个 flank 不可用、短于 250 bp、跨记录边缘、双侧指向不同位置或无法排除
重复宿主区域时，记 context_unresolved；不能把每个未解决项都算作新背景。
图中非传递或多对多关系也保留为未解决，不用连通性强行归组。
无 flank 命中只表示在本次检索条件下未检出；区分“可评估但未匹配”和“不可评估”。
新背景的人工审计同时查看全基因组邻域，不能仅靠 contig 名或坐标变化。

存在同源位置时，允许报告 present / observed_empty / altered / unresolved：
- observed_empty 必须在另一株真实 assembly 中定位到对应连续区间，
  并由该株 reads 检查连接；不得把两段 flank 人工拼接当作观察到的空位点。
- 没有 native family call 不等于该位点无 X。
- gap、覆盖不足、复杂替换及断裂不能记录为明确缺失。
- 有/无比较不推断转移方向、发生时间或活性。

## 7. P4：Nanopore 连接证据按菌株独立核查

复用当前 minimap2 全 assembly 比对思路；冻结实际版本、preset、
索引及参考哈希、完整命令。基准 map-ont 保留；
只有元数据支持的高准确度输入才采用另行登记的 lr:hq（minimap2 >=2.27）。
保留 secondary/supplementary 与逐 read CIGAR、竞争位置、AS/NM/MAPQ。
候选 triplet 不作为唯一参考。达到 secondary 输出上限者按原计划的
敏感性分支检查；仍无法排除竞争时标 ambiguous，不凭 primary/MAPQ 单独证明唯一。

本轮 baseline 与审核版分开：
- v1：保留现有 classify_reads.py 及其结果，不能覆盖。
  现有 full_span 是至少一个支持单位跨 candidate 和左右各 100 bp，
  它不是全部 500 bp 侧翼的跨越，也不是精确末端证明。
- 审核版：至少两个身份可区分的 read 单位，各自在单条连续比对中跨过
  candidate 和两侧各 100 bp；检查跨 junction 的实际 CIGAR。
  初始操作性要求为 candidate 参考碱基的 M/= /X 覆盖 >=90%，
  两侧 100 bp anchor 各 >=90%，candidate 内 I+D 总长度 <= candidate 长度的 10%。
  保存所有实际值和不通过原因；这些是连接报告标准，不回头删除 discovery 家族。
  exact command、CIGAR 分母及端点归属规则须随实现和对照一起冻结。

Nanopore read 数不能无条件等同原始独立分子数：
检查原始 read ID、重复导出、split parent 与 duplex/parent 信息；
避免对同一原始证据的多个派生记录重复计数。
不能确认关系时标 independent_molecule_count_unverified，
同时保留 read-level 支持，不把缺少元数据解释成已独立。
同一证据单位在同一位点出现相互矛盾结构时单列 conflict，
不通过旗标集合合并成额外的完整支持。

审核版先在 AR20-38 上重算，与 v1 逐成员对照，再应用于新增菌株。
保留 full_span、左右分别支持、单侧、内部、歧义、无可用覆盖及边缘受限等状态。
来自组装所用同批 reads 的支持标 assembly_consistency，不称独立实验验证。
组装中不存在的低频插入另列，不与本轮 assembly locus 总数相加。
若报告 allele fraction，分母必须来自同一位点的可分辨支持/参考跨越单位，
不得用全基因组平均深度代替。

## 8. 可选 P5：用已发现 X 搜索其它 assembly

严格逐株 native 汇总会漏掉“某株只有一份 X，未被原生重复发现器报出”的成员。
先交付 P1–P4 基线，再可用已发现的真实 X 代表对所有面板 assembly 做同源检索，
单列 homology_recovered；它是家族引导补查，不是新的从头发现。
沿用 P2 对应标准并保留 partial/ambiguous 命中，恢复 flank 后仍逐株核查 reads。
增加的 locus/context 数与 native-only 结果分别报告。

如果所有菌株都只有一份且没有任何株发现 X，本支路仍可能完全漏掉它。
这类问题由后续真正的 pooled discovery 对照评估；
必须先通过原计划的 toy gate、基因组规模 spike-in、固定 l 与资源阶梯。
近缘侧翼共同相似可能破坏 PILER 的 isolated 条件，不能预设 pooled 一定更好。

## 9. 对照、验收和停止条件

先完成技术对照，后解释新菌株增量；对照均为计算数据，不涉及生物构建。

| 对照 | 必须满足的行为 |
|---|---|
| 同一成员重复报告、GCA/GCF 重复来源 | alias 增加，sample/locus/context 不增加 |
| 同一同源位置出现在两株 | sample_locus_count 可增加，context 组不因此翻倍 |
| 同一个 X 位于三个独立随机背景 | 在真值已知且完整的比较输入中恢复 1 个 X 组和 3 个背景组 |
| A–B–C 链式相似、共享短片段和嵌套 | 不无条件合并成一个严格 X 家族 |
| 反向互补、contig ID 重名、边缘截断 | 来源和切片回写正确；截断不算完整背景 |
| 同一 read 派生记录重复、竞争比对、CIGAR 跨界缺失 | 不增加独立支持；不能伪装成完整跨越 |
| assembly 实际有 X 但 native 无 call | 不标 observed_empty；可在 P5 单列补查 |

技术 go/no-go：所有来源映射、坐标/切片和上述已知真值计数检查必须通过；
失败时暂停跨株归组和生物学解释，保留失败记录。
AR20-38 冻结输入的原生输出若被重新生成，须解释与 629d37d 的任何差异；
新版报告分类允许变化，但必须给出逐成员差异表。

最小生物学目标：在至少三个匹配菌株的面板中，找到至少一个严格 X 组，
出现于至少两个独立菌株，且经 P3/P4 审计后的受支持背景数比 AR20-38
自身基线至少多 1。达不到不改阈值追阳性，报告“未证明新增背景”，
并区分数据不足、native 未召回、边界/归组歧义和仅共享既有同源位点。
N=50 不作为承诺；完整报告每个家族的 N、增量与未解决数。

人工审计至少覆盖：IS5 基准组、一个非 tnp 组（若有）、
一个跨株新增背景、一个共享同源背景、一个归组歧义项（若有）。
某一类不存在时明确写 none observed，不强造示例。
tnp-negative 不作为“新功能”的证据；本计划不预设功能或新颖性结论。

## 10. 交付和执行安排

计划中的新分析放在独立实验目录，不先修改 package 主发现逻辑。
先复用并核对已有 AR20-38 文件，不重复下载已验证的大文件。
每次重算使用新目录；原实验目录和 v1 小结果保留。

应交付的文件（目前尚未生成）：

| 文件 | 内容 |
|---|---|
| panel_candidates.tsv / resolved_manifest.json | 候选全表、排除原因、冻结来源及文件版本 |
| raw_member_aliases.tsv / sample_loci.tsv | 原始成员、跨工具重复、位点对应和未解决状态 |
| family_pair_matches.tsv / global_families.tsv | 比对覆盖/identity、归组规则、歧义及真实代表 |
| context_matches.tsv / context_groups.tsv | 双侧命中、共线关系、同源背景和未解决原因 |
| member_read_evidence.tsv | 逐 read 身份与来源、候选边界版本、连接/竞争/冲突证据 |
| family_context_summary.tsv | raw records、样本、位点、背景、支持背景及相对基准增量 |
| sequences/ 与 coordinates.tsv | 每成员实际 X、up/down、triplet、方向和边界状态 |
| audit.md / provenance.json | 人工案例、失败、重跑、软件/输入哈希与资源 |

导出原有 500 bp flank 并保留截断长度；如使用 2 kb 扩展上下文，另存并标注。
人工拼接 up+down 如需展示必须标 artificial，不输出为 observed empty allele。

执行顺序：P0 样本核查/冻结 → 技术对照与 AR20-38 基线审计 →
P1 逐株发现 → P2/P3 家族与背景比较 → P4 本株 reads 核查 →
汇总增量；P4 比对可与 P2/P3 独立并行，但最终统计需关联同一冻结坐标版本。
P5 和 pooled 对照在基线报告后单独编号。

初始资源预算：每株原生发现 2 CPU / 32 GiB / 2 h；
每株 reads 核查 16 CPU / 64 GiB / 6 h。跨株比较先以两株测量，
初始 4 CPU / 32 GiB / 2 h；记录临时磁盘和每任务 peak RSS、时间、退出码。
预算是本轮配置，不是运行性能保证；先检查存储可用量再下载。
OOM/超时停止该分支；任一资源到预算 80% 时复核后记录调整。
扩大到 3–5 株前先检查两株结果和实际资源，不把全数据库一次混入。

## 11. 方法依据

- [本仓库 METHOD](METHOD.md)：native/candidate 区分、边界状态和原生输出保留。
- [首轮实验代码](../experiments/read_callout/)：已跑通的单株 assembly＋reads 基线；
  不代表本计划的跨株归组或审核版 read 判定已经实现。
- [PILER 官方手册](https://www.drive5.com/piler/piler_userguide.html)：
  多序列输入与 intact/isolated 定义；近缘背景带来的限制。
- [minimap2 官方仓库](https://github.com/lh3/minimap2)：
  完整参考比对、preset 与多重比对输出。
- [NCBI Datasets](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/command-line/datasets/download/genome/)：
  assembly、注释和元数据获取。
- [Nanopore 输出字段说明](https://software-docs.nanoporetech.com/output-specifications/26.01/protocol_formats/sequencing_summary)：
  read/parent 与 duplex 来源关系。

本计划是对现有单株结果的下一步扩展，不是宣称跨株路线已验证。
