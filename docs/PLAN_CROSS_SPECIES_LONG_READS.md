# 跨物种联合发现与长读段评估计划

日期：2026-10-22（运行前审查修订 v4）  
状态（2026-10-23 更新）：四株 assembly＋reads 连接证据首轮实验已完成，见
[READ_CALLOUT_RESULTS.md](READ_CALLOUT_RESULTS.md)；联合 FNA 与长读段从头重复发现尚未完成。
当前下一步优先采用[近缘 Paraburkholderia 逐株发现与跨株汇总](PLAN_PARABURKHOLDERIA_MULTI_ISOLATE.md)。
下文 A 的 pooled 发现及其对照保留为后续实验；历史修订记录不代表当前执行状态。  
适用仓库：KuangHu/longread_inserton（multisite-insertion）。

## 1. 要回答的问题

目标结构是多个独立背景中的：

`[upstream_A][X_A][downstream_A]`  
`[upstream_B][X_B][downstream_B]`  
`[upstream_C][X_C][downstream_C]`

X 的内部相似，两侧背景不同。拷贝可以来自同一菌株，也可以来自不同物种；
不以 perfect empty site、已知 tnp、TIR 或 TSD 为从头发现的先决条件。

本轮分别检验：

1. 多个 assembly FNA 联合搜索，能否找回逐基因组搜索遗漏的家族，尤其是每个基因组仅有一份 X 的情况？
2. 长读段从头重复发现工具能否在细菌数据上提供有效候选？
3. 原始长读段能否支持每个候选的 insertion–flank 连接，核查组装折叠、断裂或错误连接？
4. 哪些结果只有重复序列证据，哪些另有完整 tnp、独立边界或宿主归属证据？

“同一类 X 出现在不同物种背景”不直接等于观察到水平转移事件；
不推断转移方向、发生时间或转座活性。

## 2. 已实现范围与未验证部分

当前 README 和 METHOD 描述的是 assembly-based 流程：

- 每个 sample 是一个 isolate 的染色体和质粒集合。
- PALS/PILER-DF、维护版 RepeatScout 按 sample 独立发现。
- 同一次运行可共享 deposited protein 注释库，但这不是 DNA 层面的联合发现。
- 当前 package 不读取 FASTQ/BAM，不执行长读段 SV calling。
- 原生发现、注释、候选边界推断与序列导出已有分层。
- 既有三株回归验证支持重构一致性，不代表跨物种或宏基因组准确率已验证。

此前 cluster 审计结论来自用户提供的报告；本计划未重新检查 cluster 原始文件。
审计提醒必须延续：native 区间、candidate 区间、consensus、CDS 长度与重叠长度分别记录。

本计划是显式扩展实验，不恢复已放弃的旧 reads → SV → empty/filled 主路线。
旧路线不因长读段重新纳入评估而自动启用。现有 CLI 不应被描述为已经有 pool 或 reads 模式。

## 3. 工具证据与选择

| 工具 | 已有输入与方法证据 | 本项目用途 | 不能直接外推的部分 |
|---|---|---|---|
| PILER-DF | 官方支持多条 DNA 序列；寻找默认至少三个完整、isolated 的重复拷贝 | 联合 FNA 的主测试；最接近“内部相似、外部不同” | 多物种输入可行不等于细菌跨物种效果已经验证；一般重复不必是 IS |
| 维护版 RepeatScout | 种子向外延伸构建重复家族；有序列边界保护，1.0.7 可记录成员范围 | 联合 FNA 的并行对照 | 不显式要求两侧独特；原生成员范围和 consensus 不等于完整元件边界 |
| RepLong | PacBio reads 的重叠网络与群组识别，构建从头重复库；论文在人和果蝇测试 | 长读段从头发现候选工具 | 未核实其 Nanopore、细菌混合样本或现代 HiFi 适用性；安装和实际输出需先测试 |
| LongRepMarker | 有仅长读段的从头模式，使用组装和重叠/多重比对；论文包含 PacBio 测试 | 第二个长读段重复发现对照 | 广义重复库不等于逐位点 insertion＋双侧 flank；未验证跨物种宏基因组能力 |
| InMut-Finder | 已知 target 序列在 Nanopore/PacBio reads 中定位，提取 flank 并映射参考 | 可选的已知元件连接验证对照 | 需要已知插入序列和参考；不是未知家族发现器 |
| TELR | 长读段对参考调用插入，结合 TE 库、局部组装和 flank 定位 | 可选的已知 TE 位点对照 | 需要参考和 TE 库；non-reference insertion 不表示 reference-free |

不把 TELR、InMut-Finder 的已知元件定位表现作为未知家族发现成功。
不将任何上游工具描述成已经能直接完成本项目全部输出。

## 4. 共用实验约束

- 保留所有原生家族、成员和日志。先看软件自身给出的结果，不增加自定义候选淘汰规则。
- 第一次比较使用原生默认发现参数。记录二进制版本、提交、完整命令和实际默认值。
- 如果默认失败，保留失败；后续参数实验单独编号，不覆盖基线。
- 工具自带限制与本项目后处理分别列出。PILER 默认家族大小约束、PALS 相似度/长度约束都可能影响跨物种召回。
- PILER consensus 提取默认成员上限不等于家族实际成员数；清点全部原生 GFF 成员。
- 方法返回错误、超时、缺输入、注释缺失与真正零命中分开。
- 每次运行使用新目录，保存输入、注释、原生输出和导出的哈希。
- 冻结现有 tnp-supported 定义；同时保留全部家族，不能以 tnp 子集代替未知元件发现评估。
- 没有 tnp 注释或蛋白命中不代表没有相关功能；/product 与 /note 证据分别报告。

## 5. 阶段 A：多个 FNA 联合发现（保留为后续对照）

### A1. 输入及来源追踪

先使用少量已完成单株评估的环境分离株，随后扩展不同物种。
八株候选、定位 accession、来源及预注册 pool 组成已登记在
[CANDIDATE_ASSEMBLIES.md](CANDIDATE_ASSEMBLIES.md)。这是候选登记表，
下载后仍须生成含 accession.version、BioSample、文件哈希的 resolved manifest；
不能默认它们都有可用配套长读段。

**先检查 reads，再冻结输入。** 首轮元数据检索已登记在
[LONG_READ_AVAILABILITY.md](LONG_READ_AVAILABILITY.md)：DdD110、RoMoAcy1、
AmHOT1A3 已找到匹配的 PacBio run；其余五株在已查范围内未确认匹配长读段，
不等于证明不存在。保留原八株作为 A 的候选，C 优先使用这三株。
独立 Nanopore 对照 PaAR2038 不属于八株，不替代同种的 PaBN5。
BioSample 匹配与文件下载/校验是不同状态；实际下载、版本和哈希锁定后才冻结 manifest。

AltDE、UT26S、RHA1 均排除在本轮 scale_2/4/8 及其 spike-in 输入之外，
仅单独作回归检查。ScL1 与 UT26S 同属 Sphingobium，RoMoAcy1 与 RHA1
同属 Rhodococcus；亲缘关系不是重复样本，但不能据此把回归株静默加入 pool。
增加近缘株须另行预注册组成及同属背景标签。

每个 discovery pool 记录成员 assembly、BioSample、版本、物种标签及标签来源。
同一样本的不同组装或同一菌株的重复提交不能无说明地算成独立宿主。

将 FNA 合并为 multi-FASTA，保留每条 contig/复制子的边界；
不要把不同记录的碱基直接拼成一条序列，不人为创建跨 contig 的连接。

准备 source map，至少包含：

`pool_record_id, source_sample_id, assembly_accession, biosample, original_record_id,
sequence_length, sequence_sha256, taxonomy_label, taxonomy_source`

新 ID 必须全局唯一。原始序列不改动；若同步更改注释 record ID，
保存原注释和映射，逐条验证序列及 feature 坐标仍正确。
跨来源匹配按 source map 解释，不能把 pool 名当作宿主名称。

先在原生工具层运行；不要通过把多个 isolate 偷换成一个现有 sample、
丢失来源语义的方式宣称 package 已支持联合发现。
后续若增加 package pool 接口，需独立记录设计和实现。

### A2. 参数与最小对照矩阵

#### A2.1 RepeatScout 的输入规模依赖

“同样使用 native defaults”不表示逐株和 pooled 的实际参数相同。
官方文档给出默认 l 约为 ceil(1 + log4 L)；本次检查维护版源码发现，
build_lmer_table 和 RepeatScout 均将 FASTA 的文件字节数（ftell）传给 default_l，
不直接使用 A/C/G/T 碱基数。因此标题、换行和合并都会影响估计，
尤其在整数边界附近。实际 cluster 二进制仍需与其对应源码核对。

每次运行记录 input_bases、input_file_bytes、record_count、
l_table、l_repeat、l_source、版本/提交及日志路径。能从日志获取时以日志为准；
用对应源码推导的值必须标为 source_derived，并与输出的 l-mer 长度核对。
两个阶段使用相同文件和相同 l，否则停止该运行。

保留两个独立对照系列：

- native_default：每个输入自行使用原生默认值，报告所有实际 l。
- fixed_l：在查看 discovery 结果前，以冻结的最大计划输入（包括预注册 spike-in pool）的默认 l
  确定 l_fixed；逐株、2/4/8 pool 均在两个可执行程序中显式传同一 -l。
  若源码/日志核对表明值无效，修订配置后重新编号，不临时按结果选值。

两系列分开解释；只有 fixed_l 系列才能在固定 l 条件下比较 pool 效应。
不覆盖原生默认结果，不增加 family 过滤。

#### A2.2 数据对照

| 对照 | 输入安排 | 要检查的内容 |
|---|---|---|
| 单株基线 | 各 FNA 分别运行 | 重现当前发现行为 |
| 联合发现 | 同一批 FNA 作为一个序列集合，包含跨记录比对 | 比较原生家族与成员；不能仅拼接独立运行结果 |
| 技术阳性及限制对照 | 下述 2/3 记录 × 100/97/93% 相似度矩阵，各记录一份 X | 测量跨记录召回、局部命中及边界误差；不是天然 IS 证据 |
| 覆盖/重复提交对照 | 同一序列记录的重复副本，明确标记人工来源 | 检查重复记录不会在报告中变成多个独立生物学位点 |
| 保守背景对照 | X 与较长侧翼共同相似的记录 | 观察软件如何处理保守宿主片段，保留结果 |
| 天然案例 | 有独立依据支持的自然重复家族 | 检查实际边界和侧翼表现，避免只验证合成例 |

技术矩阵预注册如下：

- 基准为 1,200 bp 的随机非编码 X，每个记录只有一份；
  左右各 2,000 bp 独立随机背景；固定并发布种子 11、29、47，共三个重复。
- 记录数取 2、3；X 的目标两两 identity 取 100%、97%、93%。
  输出真实两两 identity 和错配分布，不能把各序列对祖先的 identity
  冒充序列间 identity。突变矩阵与 truth BED 在运行工具前冻结。
- 每个单记录输入也运行一次，组成逐记录和 pooled 对照。
  在本矩阵内另做 fixed_l 系列，l_fixed 从最大技术输入预先确定，
  不与自然数据最大 pool 的 l_fixed 混用。
- PILER 默认至少三个拷贝，因此两记录组是限制对照；
  三记录、100% identity 组才是强制技术阳性。
- PALS 原生默认保留至少 400 bp、至少 94% identity 的局部比对。
  97% 组评估高相似跨来源敏感性；93% 组评估低于门槛时的行为。
  全长 identity 为 93% 不保证没有 >=94% 的局部岛：
  要另报局部 PALS 命中、片段家族和完整家族，不能预设为绝对零命中。
- 额外以三记录、100% identity 的 X 长度 300/399/400/401 bp
  检查长度边界，保存命中真实长度与端点；偶然侧翼匹配造成延伸也要记录。
  此系列不作为强制成功门槛。
- 所有坐标均转换为 0-based 半开区间后与 truth 对照；
  同时报 native 端点、独立推断端点、召回拷贝及多余成员。

RepeatScout 自身的种子频次等默认限制单独记录，不能认为两拷贝组
对它必然阳性或必然阴性。不得放宽参数把限制对照改成阳性再覆盖基线。

若真实数据中未找到“每个来源一份、合并后多份”的天然案例，
明确报告；不能把合成阳性当成已经证明天然跨物种发现成功。

近缘基因组加入后，PILER 的 flank 在整个输入集合内可能不再独特。
不同 pool 组成分别运行并记录，不假设添加更多基因组必然提升召回。

#### A2.3 基因组规模的计算 spike-in（天然结果解释前必做）

A2.2 的 toy 阳性验证技术流程，不能证明基因组规模的发现敏感性。
额外建立 scale_4 和 scale_8 的计算 spike-in 对照，各自与未改动的天然输入配对。
这是 FASTA 文件上的模拟，不涉及构建生物材料。

- 沿用 A2.2 的随机非编码 X（1,200 bp），每个重复中的三个 X 完全相同；
  分别放入 SaSB2B、DdD110、ShANA3 的主染色体，每个来源恰好一份。
  使用种子 11、29、47；每个种子的 scale_4/8 共用相同改动。
- 插入位置在工具运行前按固定随机规则选择，距记录两端至少 10 kb；
  保留真实基因组背景，不额外添加 toy 随机侧翼。生成脚本、随机算法、
  X、位置及 truth BED 在查看 discovery 输出前冻结。
  冻结前检查 X 在原始输入中的既有同源命中；检查方法、阈值及遇命中时
  的确定性重生成规则也须先登记，不能根据发现成功与否选择 X 或位置。
- 只修改隔离的模拟文件，天然 FNA 不改动。保存新记录 ID、前后哈希、
  修改清单及分段坐标映射；插入片段没有原始坐标，后续区段相应偏移。
  若触及 CDS，明确原注释已被打断，不能将其描述为完整 CDS。
  本对照只测试 DNA 发现与坐标恢复，不产生 tnp 生物学结论；
  天然注释背景比例仅使用未改动的天然输入。
- 两个规模均运行 native_default 和 fixed_l，并与对应天然输入、
  逐株基线比较。fixed_l 使用 A2.1 同一个预注册值，不为 spike-in
  另调 l；native_default 则逐文件记录实际 l，包括加序列导致的变化。
- 每个工具/规模/种子/参数分支报告 X 召回数（0–3）、是否恢复共同家族、
  native 与 candidate 的端点误差和双向覆盖率、多余成员、侧翼表现及资源。
  沿用 A4 的端点/覆盖阈值作结果标记，但分别列出原生与后处理表现。

完成该对照是解释相应规模天然结果的前提，不替代 A4 的 toy 硬门槛。
若某分支未恢复三份 X，保留结果并标记 genome_scale_sensitivity_unresolved；
该分支的天然零命中不能解释成不存在重复家族。超时/OOM 记资源阻断，
不能算作生物学阴性。任何修复或参数试验另行编号，不能覆盖首次结果。

### A3. 注释、边界和计数

先比较 native discovery，再比较后处理结果。联合运行可能改变蛋白参考库：
注释比较应使用固定且有哈希的共同蛋白参考，或明确标为参考库不同、
不能直接把 tnp 标签变化归因于 discovery 改进。现有实现若不能固定参考，
先报告 native 比较，不虚构已完成的受控注释比较。

对全部家族保留清单；需要导出全部家族时明确记录使用 --subset all。
tnp-supported 作为独立报告层，保持现有定义不变。

逐成员保存 native、independent inferred 和 exported candidate 坐标。
candidate 的 tnp 完整性读取 candidate_overlapping_tnp_CDS，
列出每个重叠 CDS 及 fully_contained，不能只显示第一条注释。

原生 member 数、独立来源数、坐标去重后的位点数及物种数分别统计。
同一位置跨工具或多个范围重叠时先报告匹配关系；去重规则须预先记录，
原始行不删除。相同序列在不同样本出现不自动证明新的插入事件。

边界逐项报告 METHOD 中的五个 boundary_source 值及未解决原因：
deposited_mobile_element、provisional_homology_median、
homology_extent_censored_at_window、native_interval_only_unresolved、
invalid_native_coordinates。最后一类保留在 inventory 中，不导出无效区间。
用于评价边界时，
采用 deposited annotation 的 candidate 与独立推断结果分开，
不能用注释替换后的精确一致证明从头边界推断准确。

#### A3.1 背景重复占比：预注册报告层

在所有 native 家族上报告，另对已导出的 candidate 区间单独统计；
两种区间不混用。按原始 feature location 的每一段计算重叠，
至少 1 bp overlap 记为触及，同时保存重叠 bp 和区间覆盖比例，
避免将擦边误读成完整包含。

标签为 rRNA、tRNA、panel_single_copy_core_CDS、shared_same_genus_CDS、
other_CDS、no_feature_overlap、annotation_unavailable。
前五类可多选，不能直接相加成 100%。
no_feature_overlap 仅用于有注释覆盖的记录；缺失注释记 unavailable。

背景基因证据采用独立、预先冻结的 orthogroup 成员表；
工具版本、参数、蛋白输入、属分类依据及哈希均在 discovery 结果查看前固定。
不能按“CDS”或 housekeeping 产品名称臆测，也不能用 discovery 家族本身
循环定义保守背景。

- 主要背景标签 shared_same_genus_CDS：该 CDS 的 orthogroup 在
  本 pool 的另一个同属来源中也有成员；不要求单拷贝。
  本轮预注册两对：SaSB2B–ShANA3（Shewanella）和
  AmHOT1A3–AaH17（Alteromonas）。保存 pair_id、来源及 CDS ID。
  子 pool 仅评估两株均在其中的配对；无同属伙伴记 not_applicable，
  有伙伴但分析未完成记 not_assessed，均不能冒充零重叠。
- panel_single_copy_core_CDS 保留为次要描述：
  冻结完整八株 panel 中每个来源恰好一份的 orthogroup；
  子 pool 沿用该定义。它跨越多个属，可能很少，不能代表主要宿主背景混淆。
  未获得证据时记 core_status=not_assessed。
- 同属共享 CDS 也可能是移动元件相关基因。与 tnp-supported 标签交叉报告，
  不把 shared 或 core 标签自动解释为宿主核心功能、假阳性或排除条件。

每类至少输出：

- 有任一成员 overlap 的 family 数 / 全部 family 数；
- 已有可用注释的成员全部 overlap 的 family 数，并单列其未注释成员数；
- overlap member 数 / 全部 member 数，以及注释可评估成员分母；
- family × annotation-class 成员计数表，和每成员的证据 ID；
- 各类重叠的组合表，以及全无注释/部分缺注释的 family 数。

比较逐株与 pooled 的类别占比、绝对数量和 annotation availability，
不把 rRNA/tRNA/core/shared_same_genus 标签作为候选过滤，不从 native 输出删除任何家族。
重叠宿主基因也可能由真实插入导致，故该表仅衡量注释关联背景，
不能直接把全部 overlap 家族称为假阳性或非 IS。

### A4. 阶段 A 验收

- 所有合并记录能回溯到原始来源，序列与坐标核验通过。
- 单株和联合原生输出均完整保存，工具失败不被记录成零命中。
- **强制 go/no-go（先于任何天然 pool 运行）**：三记录、100% identity、
  X=1,200 bp 的技术阳性，在三个固定种子中，至少同一个工具的
  native-default 分支均恢复包含三个真值位置的一个家族，
  对应三个成员的 native 左右端点绝对误差各 <=10 bp，且各自与真值
  的双向区间覆盖率均 >=95%。额外成员照常报告，不删除。
  同时必须通过来源映射、正反向坐标转换和序列切片回写完全一致检查。
- 此数值门槛仅用于判断软件与输入是否可继续，不是天然候选过滤规则，
  也不承诺达到生物学精确末端。独立边界后处理不能挽救 native 分支失败来充当通过。
- 若没有工具满足强制门槛，暂停天然 pools，定位安装、输入、参数、
  原生工具或坐标实现问题；保留失败结果，不临时扩大误差容忍度。
  只有一个工具通过时可继续该工具的天然测试，另一工具单独标 blocked。
- 2 拷贝、97/93% 及短长度各组必须报告，但它们不通过不触发阳性门槛失败。
- 相应规模的天然解释前必须完成 A2.3 genome-scale spike-in；
  未通过分支标记敏感性未解决，不以 toy 通过替代基因组规模验证。
- 报告逐株独有、联合独有、两者共有的可比家族/位点及匹配依据；
  不能只减去两个家族总数推算“新增”。
- 报告完整双侧 flank、截断 flank、边界未解决、tnp 完整覆盖及跨来源分布。
- 只有天然案例通过独立核查，才声称获得天然跨来源候选；
  species 归属不清的记录保留为 unknown。

### A5. 资源阶梯：通过技术门槛后、完整天然 pool 前

先按候选登记表中的嵌套集合做 2、4、8 株阶梯，各工具分别测量。
对每个规模记录输入总 bp、文件字节、contig 数、实际参数、
wall time、CPU time、peak RSS、临时/输出磁盘、PALS 命中数、退出码。
建议用 /usr/bin/time -v 加 SLURM accounting 留存观测；不同统计口径注明。

初始每个阶梯任务的预算为 2 CPU、32 GiB RAM、2 小时 wall time、
50 GiB 工作目录；这是本轮预先选择的预算，不是工具需求保证。
应用程序内部使用线程数另记，不能把 SLURM 分配数当作实际并行度。

按 2 → 4 → 8 顺序运行，不把全部规模一次性提交。
任一步 OOM、超时、磁盘越限或未正常完成时停止升级；
RSS、wall time 或磁盘任一达到预算的 80% 时先复核，修订资源配置并记录
后才继续，不自动重试放大。完成 8 株阶梯前不启动更大的自然 pool。
scale_4/8 在资源允许后分别运行 A2.3 配对的 spike-in，使用同样预算和
升级规则；该对照未完成前不解释对应规模的天然发现结果。

all-vs-all 比较的潜在序列对随规模增加，不能用单株耗时线性外推；
也不把“时间或内存严格按总长度平方增长”当作已证实事实。
PALS 的种子、重复内容、长度/identity 阈值及构建版本会影响实际成本。

若需分块，作为新配置按官方方法保留块内 self 与所有块间比较，
再合并 hits；不能只做块内比较后声称等价联合搜索。
先在能完整运行的小规模输入上检查分块结果一致性，
保留命令和坐标核验；不静默更换 aligner 或阈值。

## 6. 阶段 B：长读段从头发现工具评估

### B0. RepLong 兼容性预检（不阻塞 A）

2026-10-22 已完成的是只读源码检查，不是安装成功：

- 上游 README 仍可读取，声明使用 Canu fork、Java 8、R/igraph、
  pyfaidx 0.4.7.1，测试环境为 Ubuntu 14.04/CentOS 7；
  这些是上游历史说明，不是本项目确认的现代兼容性。
- install.sh 会直接 clone ruiguo-bio/canu 后 make -j 4，
  未固定 fork 的 commit。该仓库目前可读取，但仍须锁定其提交；
  不仅固定 RepLong 本身。
- install.sh 检查时 blob SHA 为 1c643c526b05857fffb1bacafdf1cf3d91e5ed74。
  安装前核对实际取得版本，禁止用未记录的最新 fork 冒充可复现依赖。

先在 cluster 隔离环境预检；不得修改现有 multisite-insertion 环境。
保存系统、编译器、Java、Python、R、igraph、pyfaidx、RepLong/Canu
提交、安装日志和退出码。逐步验证依赖、构建、启动、作者小示例输出。
预检初始预算为 1 个工作日；仍无法运行则标 blocked_install，
记录错误和尝试，不让阶段 A 等待旧工具移植。

若必须改上游代码或依赖版本，保存独立补丁及环境，不把 patched run
与 upstream baseline 混为一谈。作者示例不可获取也必须单列。
只有小示例完成且输出可解析后，才安排 RepLong 的目标数据实验；
LongRepMarker 独立通过同类预检，不以 RepLong 成败代替其验证。

### B1. 工具层可运行性

RepLong 与 LongRepMarker 分别建立隔离环境，先运行作者提供且可获取的示例。
冻结版本、依赖、命令、输入类型、运行时间、峰值内存和原始输出格式。

确认每个工具到底返回 repeat library、代表性 read、区间、consensus
还是成员集合；未输出的逐位点 flank 不能描述成工具原生输出。
若官方示例不可获取或依赖无法复现，记录明确阻断原因。

作者示例只验证工具可运行；不验证细菌 IS 或跨 species 能力。
只把文献明确验证的数据类型记为已验证；尤其不把旧 PacBio、
PacBio HiFi/CCS 和 Nanopore 当作同一种已测试输入。

### B2. 目标层数据

首轮元数据预检已取得三个 PacBio 候选：
DdD110（SRR16235513）、RoMoAcy1（SRR15616890）、
AmHOT1A3（SRR2084418），优先顺序按此列出。
独立 Nanopore 对照为 PaAR2038（SRR11947722，GCA_012689525.1），
不属于八株 pool，也不是 PaBN5；仅用于 B/C 的独立环境分离株测试。
匹配证据、下载地址及未完成项见
[LONG_READ_AVAILABILITY.md](LONG_READ_AVAILABILITY.md)。

先完成文件校验、参考版本解析和输入类型确认，再安排工具目标测试。
混合微生物样本的匹配数据仍未选定，不能宣称已覆盖宏基因组场景。
优先依据数据完整性，不把不同菌株的结果解释为平台优劣。

Nanopore 和 PacBio 分别记录：

- 平台、化学版本/basecaller 或 CCS 状态（可获得时）；
- SRA/ENA accession、BioSample、是否与 assembly 同源；
- read 数、总碱基、长度分布和质量信息；
- 原始数据、纠错数据、抽样数据的关系与哈希；
- 抽样方法、随机种子和样本内覆盖不均的限制；
- 记录与分子的对应关系（subread / ZMW / CCS / Nanopore read）及依据。

抽样以分子为单位：PacBio CLR 按 ZMW（或已确认的 SRA spot）整体抽取，
不能独立抽取 subread，否则同一分子的重复读取会改变有效覆盖。
分组不可得的 run 先记录限制，再决定是否用于需要抽样的工具测试。

不同物种/不同样本之间的结果不能直接解释成测序平台优劣。
只下载/使用事先核实匹配的数据；没有配套 reads 的 assembly 留在阶段 A。

### B3. 与 assembly 路线比较

在同一数据集上比较 assembly 候选与 read-based repeat 库；
先固定并记录序列匹配规则，再报告共有、各自独有、片段化和未分类结果。

针对 repeat 库中的候选，另行恢复支持 reads 与不同侧翼背景。
这是额外分析层，不宣称 RepLong/LongRepMarker 原生已经完成。

必须区分同一位点的测序覆盖与不同位点的重复拷贝。
宏基因组中丰度不同，不能用全样本平均覆盖度直接换算各宿主的拷贝数。
如果工具无法区分高丰度单拷贝序列与真实多位点重复，记录为关键限制。

阶段 B 的成功分层记录：安装成功、工具示例成功、细菌目标成功、
混合样本可解释候选。前一层不能代替后一层。

## 7. 阶段 C：原始长读段的连接证据

这一步可在阶段 A 得到候选后进行，不必等阶段 B 成功。

### C1. 比对工具与竞争参考

使用 minimap2，对 **同一 BioSample 匹配的完整 assembly** 进行 reads 比对：
分离株包含所有染色体和质粒；混合样本包含可获得的完整 metagenome
assembly（含未分箱 contigs），不能只保留候选所在的 MAG。
若只有部分 MAG 可用，明确 reference_incomplete，不能据此确认全样本唯一定位。

候选与参考保持相同记录 ID/坐标或提供无损 source map。
只把多来源 pool 中属于该 reads 样本的完整 assembly 作为主参考；
其他样本的序列若用于竞争敏感性分析，独立列出，不能伪装成同源参考。

固定 minimap2 版本/提交、参考及索引哈希、完整命令。
预注册 preset：PacBio HiFi 用 map-hifi，旧 PacBio CLR 用 map-pb；
Nanopore 高准确度且有元数据支持时用 lr:hq，其余用 map-ont。
lr:hq 从 minimap2 2.27 起提供，使用该 preset 必须固定 2.27 或更高版本；
旧 Nanopore run 不能仅凭平台名就认定为高准确度输入。
PacBio Sequel 仪器名也不等于 HiFi/CCS，须核查 read 类型。
preset 在查看候选支持结果前根据数据类型选择并记录；
实际固定版本不支持某 preset 时先修订配置。

主映射设置显式保留 secondary：--secondary=yes -N 50 -p 0.8。
其中 --secondary=yes 和 -p 0.8 已是 minimap2 默认值，写出它们便于审计；
-N 50 将默认最多 5 条 secondary 增加到 50 条。

输出 SAM/BAM，保存 supplementary 及 unmapped 记录；
不得用只保留 primary 的处理来决定候选唯一性。
记录 primary/secondary/supplementary、MAPQ、AS/NM、CIGAR、剪切、
reference/read 起止及 read ID。junction 证据按 CIGAR 检查真实跨界覆盖，
不能仅凭 alignment 外包区间跨过两端就判定完整支持。

-N 是最多输出的 secondary 数，-p 也是候选截断条件，
保留 secondary 并不意味着枚举了所有可能位置。
对达到 50 条 secondary 上限的 reads，标记竞争位点未充分解析，
按预注册敏感性分支 -N 200 -p 0.5 重新对同一完整参考比对；
仍受上限影响则保留 ambiguous，不判 unique。
未达到上限也不能仅凭 secondary=0 就证明唯一，须结合可区分 flank。

pilot 使用容纳完整参考的单部分索引；记录 -I 并检查是否发生分片，
不要在未处理多部分索引限制时直接比较 MAPQ。
若参考过大无法单部分索引，单独设计并验证全参考竞争方案后继续。

提取的 triplet 只用于展示/局部序列核对，不作为主竞争映射参考。
允许对 reads 初步检索候选，但位点归属必须回到上述完整参考评估。

### C2. 连接分类与 circular-origin 限制

当前 package 按线性记录处理，不跨 circular origin 绕回。
已知环状复制子中，candidate 跨 origin 或所需 junction/flank 被 origin
截断者记 origin_limited；保留能观察到的单侧证据，完整双侧结论记
not_assessable，不将缺失部分解释成 no_support 或否定 insertion。
未标环状的 contig 末端记 contig_edge_limited，不能自行假定环状。

本轮不自动拼接首尾或复制参考末端以“修复”；如后续加入 circular-aware
处理，须单独验证旋转坐标、去重与 read 竞争，不能混入本轮基线。

对每个候选分别报告：

- 单条 read 同时跨越左 junction、candidate 和右 junction；
- 仅有分别支持左右 junction 的 reads；
- 单侧支持；
- candidate 内部命中但无连接证据；
- 没有可用覆盖；
- 多重比对、冲突连接或潜在嵌合，无法确定。

“两个 junction 分别有支持”与“同一分子跨越完整结构”分开。

支持计数以分子为单位，不以 FASTQ 记录为单位。PacBio CLR 中同一 ZMW 的多条
subread 对同一 junction 只计一个分子；同一分子不同 subread 的结论冲突时记为
该分子内部不一致，不拆成独立支持。分组依据按 run 预先登记：
zmw_from_read_name（DdD110）、sra_spot_as_zmw_assumed（AmHOT1A3，待确认）；
无法分组（当前 RoMoAcy1）时同时报告记录数并标 molecule_grouping_unavailable，
不能把记录数当作独立分子数。Nanopore 每条 read 视为一个分子。
宿主背景需有可区分的侧翼/更长上下文；短侧翼不足以分类时不强制给 species。

映射需考虑同一集合内的竞争位点，保留次优/多重比对信息，
不能仅把 reads 分别贴到各个 triplet 后重复计数为多个独立位点。
质量、比对和支持判定参数在分析前记录，不追着阳性结果调整。

长读段连接支持不证明精确生物学末端。
有 reads 直接覆盖空位点时可作为额外证据；没有 empty site 不排除候选。
将左右 flank 拼起来的序列只能叫人工拼接，不能叫观察到的空等位基因。

## 8. 交付与决策

建议分别建立 native per-isolate、native pooled、read discovery 和 read support
实验目录，避免与冻结回归结果混用。原始大数据留在 cluster，
GitHub 保存计划、manifest、下载说明、版本、命令、汇总和可复现脚本。

应交付：

1. 输入/source map 和固定版本清单。
2. 完整原生结果 inventory 与失败原因表。
3. 单株/联合比较表及代表性逐成员审计。
4. 长读段工具的实际输入、输出和适用性结论。
5. 每个位点的边界、flank、tnp、连接及宿主证据表。
6. 资源统计和人工核查案例；不以 family 总数代表独立 insertion event 数。

执行顺序：reads 元数据预检 → 下载校验和输入冻结 → A 技术门槛 →
2/4/8 株资源阶梯及 4/8 配对 spike-in → 天然联合解释；
B0 兼容性预检可独立进行，B 通过工具示例后再上目标数据；
有匹配 reads 即可按 C 的完整参考比对方案开展连接证据分析。

是否扩展规模取决于可解释的天然案例和资源消耗，不提前设置产量保证。
若 A 有效而 B 不适用，可采用“联合 assembly 发现＋长读段连接核查”；
若 B 有增量，则保留为独立发现支路，并记录其特有误判和遗漏。

这次文档提交不代表已执行实验、不改变默认 discovery 参数、
tnp-supported 定义、边界算法或现有导出。

## 9. 依据与查阅范围

以下论文/官方文档用于能力分类；查阅日期 2026-10-22。
软件兼容性和目标数据表现仍须按本计划实测。

1. [PILER 官方手册](https://www.drive5.com/piler/piler_userguide.html)：多序列输入、isolated family 定义、跨块比对与默认家族大小。
2. [维护版 RepeatScout](https://github.com/Dfam-consortium/RepeatScout)：种子延伸、序列边界保护、1.0.7 成员范围输出。
3. [RepLong, Bioinformatics 2018](https://pubmed.ncbi.nlm.nih.gov/29126180/)；[代码入口](https://github.com/ruiguo-bio/replong)：PacBio 从头重复库，真核数据验证。
4. [LongRepMarker, NAR 2021](https://doi.org/10.1093/nar/gkab563)；[代码](https://github.com/BioinformaticsCSU/LongRepMarker)：包含仅长读段的从头模式。
5. [InMut-Finder 官方仓库](https://github.com/jsg200830/InMut-Finder)；[论文](https://doi.org/10.1186/s12864-021-08206-9)：已知插入序列和相连侧翼的长读段定位。
6. [TELR 官方仓库](https://github.com/bergmanlab/TELR)；[论文](https://academic.oup.com/nar/article/50/21/e124/6717830)：参考、TE 库、局部组装和 flank 定位。
7. [宏基因组结构变异的长读段验证示例](https://doi.org/10.1038/s41467-022-30857-9)：支持 read-level junction 核查的价值，不证明本计划所有发现步骤已成立。

## 10. 本次审查修订对照

| 审查项 | 补充位置 |
|---|---|
| RepeatScout 默认 l 随输入变化 | A2.1：实际 l、文件字节数、固定 l 对照 |
| PILER/PALS 拷贝数与 divergence 限制 | A2.2：2/3 × 100/97/93%，附长度边界与局部命中报告 |
| 保守宿主背景占比 | A3.1：rRNA/tRNA/core_CDS 等多标签报告、注释缺失分母 |
| pooled 资源成本 | A5：2/4/8 阶梯、固定预算、停止/复核条件 |
| reads 比对器、参考与环状边界 | C1/C2：minimap2、完整匹配 assembly、secondary 上限、origin_limited |
| 八株名单不可复现 | CANDIDATE_ASSEMBLIES.md：候选、accession、pool 成员与下载锁定要求 |
| 缺数值 go/no-go | A4：三个固定种子的三记录技术阳性、端点误差及坐标回写 |
| RepLong 旧依赖 | B0：只读预检发现、隔离安装门槛与 blocked_install |

补充依据：
- [PALS 官方手册](https://www.drive5.com/pals/pals_userguide.htm)：400 bp / 94% 是局部 hit 默认门槛。
- [RepeatScout build_lmer_table 源码](https://github.com/Dfam-consortium/RepeatScout/blob/main/build_lmer_table.c)
  和 [build_repeat_families 源码](https://github.com/Dfam-consortium/RepeatScout/blob/main/build_repeat_families.c)：
  本次检查的 blob SHA 分别为 b0f22590437d9dc7d45ab8925ceb999995bc0442
  和 0f6a3cd76e59bb920d9f9ae2cb0685c771e290fb；实际运行仍需锁定对应 commit。
- [minimap2 官方手册](https://github.com/lh3/minimap2/blob/master/minimap2.1)：
  preset、secondary 输出、-N/-p 限制与多部分索引的 MAPQ 注意事项。
- [RepLong README](https://github.com/ruiguo-bio/replong/blob/master/README.md)
  和 [安装脚本](https://github.com/ruiguo-bio/replong/blob/master/install.sh)。


## 11. 第二轮审查修订（v3）

| 审查项 | 修订 |
|---|---|
| toy 阳性不代表真实规模敏感性 | A2.3：三来源计算 spike-in，scale_4/8 × 两个 l 系列，天然解释前必做 |
| 全 panel 单拷贝 core 太严格 | A3.1：同属共享 CDS 为主要背景标签，global core 保留为次要描述 |
| reads 核验时机过晚 | A1/B2：输入冻结前核查；三株 PacBio 已匹配，另列独立 Nanopore 对照 |
| 回归株和近缘株边界不明确 | A1 与候选表：三回归株排除在所有本轮 pools 及 spike-in 外 |
| minimap2 默认参数与版本 | C1：-p 0.8 为默认，-N 50 改变默认上限；lr:hq 需 >=2.27 |
| 边界状态数量 | A3：明确列出 METHOD 五个 boundary_source 值 |
| D110 / 110 别名 | 候选表与 reads 核查：以 BioSample 关联 D110、110、D.den110 |

本轮新增核查只涉及公开元数据与软件文档，未下载 reads、未执行发现或比对。
[minimap2 NEWS](https://github.com/lh3/minimap2/blob/master/NEWS.md)
记录 2.27-r1193（2024-03-12）加入 lr:hq。

## 12. 第三轮审查修订（v4）

| 审查项 | 修订 |
|---|---|
| CP 定位株的 assembly 未解析 | 候选表：ScL1 → GCA_000147835.3，RoMoAcy1 → GCA_019856255.1，含复制子版本与长度 |
| PacBio 读段单位与 HiFi 假设 | LONG_READ_AVAILABILITY 读段单位预检：三 run 均为 CLR；D110 质量为占位符 |
| subread 重复计数 | C2：以分子（ZMW/spot）计连接支持；无法分组的 run 显式标记 |
| 抽样单位 | B2：按 ZMW/分子抽样，不独立抽 subread |
| ENA read_count 含义 | 与 SRA spot 数而非 FASTQ 记录数核对 |
| 下载核验可复现 | `scripts/fetch_long_reads.py` 与 `scripts/long_read_runs.json` |

本轮只读取了各 FASTQ 前 4 MiB 和公开元数据；完整 reads 尚未下载。
