# Paraburkholderia 属级面板：逐株发现与 Nanopore 核查，跨株汇总家族

日期：2026-10-23  
状态：计划 v3；已完成属级公开数据可用性核查、19 份 assembly 对 AR20-38 及两两 ANI、
克隆归并与组装连续性分层；尚未冻结面板、下载新增 raw reads 或运行新增菌株 discovery / 位点比较。  
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

主路径改为显式的属级跨物种面板；同种比较是有数据时才报告的独立分层，
不再把 near_species 写成可有可无的后备。目标仍为 AR20-38 加 2–4 株，
但数量必须服从数据和来源约束，不能拿目标数替代已核实的可用性。

2026-10-23 的[完整可用性核查](../experiments/metadata_audit/paraburkholderia_20261023/README.md)
检出 24 个 genomic Nanopore runs；19 个 BioSample 可匹配 current GCA，
其中只有 3 个 Complete Genome（含 AR20-38）、1 个 Chromosome、15 个 Contig。
P. aromaticivorans 仅检出既有 SRR11947722。
P. fungorum 的 9 runs 中，7 个匹配 contig-level assembly，另 2 个未找到关联 assembly；
不能直接把它设为“完整 assembly＋reads”主面板中心。

| 暂定顺序 | 菌株 / assembly | ANI 到 AR20-38 | AF：候选 / AR20-38 | 来源及待办 |
|---|---|---:|---|---|
| 0 | AR20-38 / GCA_012689525.1 | 100.00% | 100 / 100% | 已有基准，森林土壤 |
| 1 | DD10 / GCA_049177125.1 | 86.54% | 54.08 / 38.34% | 草甸土壤菌丝际；BioSample 使用 MISAG single-amplified 模型，须核实是否培养分离株及扩增历史 |
| 2 | CT11-68 / GCA_047325075.1 | 84.88% | 52.87 / 33.15% | 小麦根部、NCSU Phytotron 植物生长设施；单列 plant-associated controlled-environment 来源 |

这是数据可用候选顺序，不是已冻结的三株合格面板。
DD10 的样本模型与 strain/assembly 描述存在需要解释的差异；
CT11-68 不能标为野外分离株，也不能仅凭生长设施来源断言为工程化菌株。
若坚持“野外环境培养分离株＋Complete Genome”，当前证据不足以宣布满足三株目标。
来源未解决时保留候选状态，不提前下载其 reads 或运行 discovery。
NCBI assembly level 不再单独决定纳入；按 3.1 的实测连续性分层，
fragmented 层仍只作探索，不静默并入主面板。
同种无新增合格样本的事实保留；PaBN5 不可借用 AR20-38 reads。

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
本次已用 skani 0.4.0 对全部 19 份匹配 FNA 测量 ANI，命令、二进制和输入哈希见核查目录。
ANI 是可比同源部分的估计，必须同时报告双方 aligned fraction；接近 80% 的结果谨慎解释。
锁定 ANI 工具版本、命令、双方比对覆盖与 taxonomy 来源；ANI 仅用于排序/分层，
不凭高 ANI 删除可能带有不同插入的独立菌株。
在同一来源/assembly 合格层内按到 AR20-38 的 ANI 降序选取，平局按
assembly accession 排序；不依据 tnp 数、IS5 有无或候选产量选样。
优先接近 AR20-38 是预先选择的设计，会偏向共享位置，报告中说明此偏差。
来源核实后的纳入/排除与最终顺序须在查看 discovery 产量前固定。
保留全部候选及未入选原因，便于评估选择偏差。

至少一个新增独立菌株匹配成功，才进入双株技术检查；
有三个合格独立菌株才称为完成目标中的最小多株 pilot。
不足三株时报告有限样本结果，不静默放宽样本类型或拿重复 run 凑数。

### 3.1 克隆归并、连续性分层与同种组（2026-10-23 实测）

数据：[ani_all_pairs.tsv](../experiments/metadata_audit/paraburkholderia_20261023/ani_all_pairs.tsv)
（skani 0.4.0 triangle，与对 AR20-38 核查相同的 `--slow -s 80 --min-af 15`，171/171 对有输出）、
[genome_groups.json](../experiments/metadata_audit/paraburkholderia_20261023/genome_groups.json)、
[panel_tiers.tsv](../experiments/metadata_audit/paraburkholderia_20261023/panel_tiers.tsv)。
19 份 FNA 在独立目录重新下载，SHA256 与核查记录 19/19 一致；对 AR20-38 的 ANI/AF 复现一致。

**克隆归并**（单链接，ANI ≥99.9% 且双方 AF 的较小值 ≥95%）：
PRJEB88388 的 P. fungorum R-93128、R-92973、R-91495、R-92561、R-92583 两两 ANI 100.00%、
AF 99.8–100%，视为一个克隆，不计为五个独立菌株。代表按“contig 最少→最大 contig 最长→accession”
选 R-93128；其余记为 alias，不删除原始记录。同项目 R-91497 不属于该克隆。
19 份 assembly 归并后为 **15 个独立基因组**。阈值是本面板的操作性定义，不是菌株分类标准。

**连续性分层**（以实际 FNA 计，不依赖 NCBI level 标签）：

| 层 | 规则 | 独立基因组 |
|---|---|---|
| complete | NCBI Complete Genome | AR20-38、DD10、CT11-68 |
| near_complete | ≤10 条记录且最大记录 ≥3 Mb | R-92578、LjRoot223、SMT4a、FW306-02-B、R-93128、R-91497、LjRoot267 |
| fragmented | 其余 | R-92571、PR3、LjRoot241、LjRoot22、WTPI1 |

complete＋near_complete 共 10 个独立基因组，可作为主面板候选池（仍须来源核实与 reads 文库审计）。
near_complete 的 contig 末端位点照常标 edge_limited / contig_edge_limited，
不计为完整背景；记录数阈值在查看 discovery 结果前固定，不按产量调整。

**同种组**（单链接，ANI ≥95%）：P. fungorum（7 份：克隆 5＋R-91497＋PR3）、
DD10＋P. terricola SMT4a＋LjRoot223（DD10 与二者同组，与其元数据中的 terricola 描述一致；
保留档案分类 Paraburkholderia sp.）、P. hospita（3 份）、P. xenovorans（R-92571/R-92578，97.33%）；
AR20-38、CT11-68、FW306-02-B、WTPI1 各自单独成组。AR20-38 没有同种伙伴。

由此 P3 应同时报告两种比较，二者不混合：
- **跨物种（以 AR20-38 为中心）**：complete＋near_complete 候选与 AR20-38 的双方 AF 仅约 25–57%，
  AR20-38 的大部分候选位点预计落在无共线对应的区域，按规则记 unresolved。
  可比较分母小是预期结果，不是流程失败；报告必须给出分母。
- **同种内**：terricola 组（DD10、SMT4a、LjRoot223，两两 ANI 99.22–99.54%，AF 88.9–95.4%）、
  fungorum（R-93128 克隆代表与 R-91497，ANI 98.51%，AF 84.2/88.1%；fragmented 层 PR3 另列）等，
  共线性高，最适合区分“共享同源位置”与“新背景”。这是原问题在现有数据下更可解析的形式，
  但不同于“与 AR20-38 同种”的原始目标，须单独命名。

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

## 6. P3：跨物种同源位置判定；未比对上不算新背景

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

取消“flank identity >=90% 即同源，否则视作不同背景”的规则。
P3 以双向全基因组比对中的宿主共线区段和邻近基因顺序为主要依据；
500 bp flank 及另存的 2 kb 邻域是局部证据，不能用全基因组 ANI
直接换算一个局部 identity 阈值，也不能把 ANI 当作局部同源位置的证明。

在新增 discovery 前，先用完整 assembly 做成对宿主比对及基因顺序检查。
首个技术试验使用 MUMmer4 nucmer --maxmatch，交换 reference/query 各运行一次，
保留未过滤 delta 与 show-coords 输出；1-to-1 筛选仅为辅助视图，不能抹掉竞争位置。
邻近 CDS 的同源与顺序需另存证据，不凭产品名相同判为同源。
固定实际二进制版本、其余参数、共线区段和唯一锚点的判定规则。
工具配置必须在两株技术检查和下述分化对照通过后存入冻结配置；
本次只运行了 ANI，并未实现或验证该位点比较层。
若在这两个较远物种间不能可靠定位宿主同源区域，P3 停留在 unresolved，
不以放宽一个 identity 数字代替方法验证。

预期结果须事先写明：跨物种（ANI 约 80–89%、AF 约 18–57%）比较中，
大部分位点没有可用宿主共线对应，unresolved 占多数是方法限制的正常体现；
同种组（3.1）才是检验“共享位置 vs 新背景”判定的主要可解析对照。
两者的分母、resolved 数和 unresolved 数分别报告。

每对同家族 locus 保存以下三态及原始证据：
- same_position：两侧宿主锚点在双向对应的共线区域内，方向/基因次序一致，
  且保留所有竞争位置后仍能唯一定位；允许 X 的存在/缺失及局部长度变化。
- distinct_position：有积极的位置证据，例如双方分别落在已解析且不同的
  宿主同源位置，或株内为已排除重复报告的不同物理位置。
  跨物种时必须能区分“位置不同”与“宿主区域无法比较”。
- unresolved：没有命中、只命中一侧、多重匹配、双侧矛盾、未比对的附属区域、
  flank 短于 250 bp、记录/origin 边缘受限、局部重排无法唯一定位等。
  原来的“可评估但未匹配”也不进入新增背景计数。

只有无矛盾的 same_position 关系可合并；非传递和多对多关系保持未解决。
不同组之间也不能仅因缺少 same_position 边就断言 distinct_position。
N_context 仅针对关系充分解析的可比较子集计算，同时交付成对三态表、
参与该计数的成员列表及所有被排除的 unresolved 项。
仍可能合并的孤立组只能作为候选组列出，不作为已证实不同背景。
株内不同物理位点与跨株位置对应分别报告，避免跨物种未知关系膨胀总 N。

按 same-species / cross-species 分层，逐菌株对报告可定位比例和三态数量；
没有同种比较时明确写 no additional matched same-species isolate。
人工审计查看完整邻域和竞争位置，不以名称、坐标变化或未命中证明新背景。

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
公开 SRA 导出的 Nanopore FASTQ 可能已替换原始 read ID（基准 SRR11947722 即如此），
因此预期多株会缺少 split/duplex 来源关系；这是档案格式限制，不是逐株实验失败。
每个 run 仍实际检查，不假定所有提交格式相同。
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
| 同一同源位置出现在两株，包括宿主侧翼低于 90% identity | 能凭双侧共线/基因顺序恢复时归为同位点，否则 unresolved；绝不因阈值失败增加背景 |
| 跨物种宿主区域未比对或多重匹配 | unresolved 数增加，已证实新增背景数不增加 |
| 同一个 X 位于三个独立随机背景 | 在真值已知且完整的比较输入中恢复 1 个 X 组和 3 个背景组 |
| A–B–C 链式相似、共享短片段和嵌套 | 不无条件合并成一个严格 X 家族 |
| 反向互补、contig ID 重名、边缘截断 | 来源和切片回写正确；截断不算完整背景 |
| 同一 read 派生记录重复、竞争比对、CIGAR 跨界缺失 | 不增加独立支持；不能伪装成完整跨越 |
| assembly 实际有 X 但 native 无 call | 不标 observed_empty；可在 P5 单列补查 |

技术 go/no-go：所有来源映射、坐标/切片和上述已知真值计数检查必须通过；
失败时暂停跨株归组和生物学解释，保留失败记录。
AR20-38 冻结输入的原生输出若被重新生成，须解释与 629d37d 的任何差异；
新版报告分类允许变化，但必须给出逐成员差异表。

生物学结果采用连续统计，不再设置“新增 >=1 背景即成功”的通过线。
每个严格 X 家族逐株报告：原生/补查来源、位点实例数、株内不同位点数、
跨株共享位置数、有积极证据的新增位置数、unresolved 数、可比较位点分母，
并分别给出 read-level 支持和分子独立性状态。
同时报告相对 AR20-38 的逐株增量，以及按冻结加入顺序的累计增量。
这些增量应由同一版最终对应关系计算；由于新增证据解决旧歧义产生的计数变化
单列 resolution_change，不混作新株新增位置。顺序会影响边际增量，始终同时列全体总量。

汇总每株/每家族的新增可解析背景及 read 支持背景；
零增量、家族不共享和无法定位均是有效结果，不修改阈值追阳性。
严格 DNA 对应可能在属级面板中一个共享家族也找不到，此时报告该限制。
至少三株仅是面板规模目标，不是生物学成功证据；N=50 不作承诺。
区分数据不足、native 未召回、DNA 家族分化、边界/宿主位置歧义和共享既有位置。

人工审计至少覆盖：IS5 基准组、一个非 tnp 组（若有）、
一个跨株新增背景、一个共享同源背景、一个归组歧义项（若有）。
某一类不存在时明确写 none observed，不强造示例。
tnp-negative 不作为“新功能”的证据；本计划不预设功能或新颖性结论。

## 10. 交付和执行安排

计划中的新分析放在独立实验目录，不先修改 package 主发现逻辑。
先复用并核对已有 AR20-38 文件，不重复下载已验证的大文件。
每次重算使用新目录；原实验目录和 v1 小结果保留。

P0 已生成的原始元数据、24-run 可用性表、19-assembly ANI 与下载校验见
[metadata audit](../experiments/metadata_audit/paraburkholderia_20261023/README.md)。
下面的 discovery/read/context 文件尚未生成，resolved_manifest 也尚未冻结：

| 文件 | 内容 |
|---|---|
| panel_candidates.tsv / resolved_manifest.json | 候选全表、排除原因、冻结来源及文件版本 |
| raw_member_aliases.tsv / sample_loci.tsv | 原始成员、跨工具重复、位点对应和未解决状态 |
| family_pair_matches.tsv / global_families.tsv | 比对覆盖/identity、归组规则、歧义及真实代表 |
| context_matches.tsv / context_groups.tsv | 双侧命中、共线关系、同源背景和未解决原因 |
| member_read_evidence.tsv | 逐 read 身份与来源、候选边界版本、连接/竞争/冲突证据 |
| family_context_summary.tsv | raw records、样本、位点、背景、支持背景、逐株/逐家族边际增量、可比分母及未解决数 |
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
先两株运行只用于技术/资源检查，不按其产量更换已经冻结的候选顺序。

## 11. 方法依据

- [本仓库 METHOD](METHOD.md)：native/candidate 区分、边界状态和原生输出保留。
- [首轮实验代码](../experiments/read_callout/)：已跑通的单株 assembly＋reads 基线；
  不代表本计划的跨株归组或审核版 read 判定已经实现。
- [PILER 官方手册](https://www.drive5.com/piler/piler_userguide.html)：
  多序列输入与 intact/isolated 定义；近缘背景带来的限制。
- [MUMmer4 官方仓库](https://github.com/mummer4/mummer)：宿主全基因组比对和多重匹配。
- [minimap2 官方仓库](https://github.com/lh3/minimap2)：
  完整参考比对、preset 与多重比对输出。
- [NCBI Datasets](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/command-line/datasets/download/genome/)：
  assembly、注释和元数据获取。
- [Nanopore 输出字段说明](https://software-docs.nanoporetech.com/output-specifications/26.01/protocol_formats/sequencing_summary)：
  read/parent 与 duplex 来源关系。

本计划是对现有单株结果的下一步扩展，不是宣称跨株路线已验证。
