# 跨物种联合发现与长读段评估计划

日期：2026-10-22  
状态：计划，尚未执行本计划中的联合发现或长读段实验。  
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

## 5. 阶段 A：多个 FNA 联合发现（优先）

### A1. 输入及来源追踪

先使用少量已完成单株评估的环境分离株，随后扩展不同物种。
既有八株候选名单仅作为 assembly 扩展候选，不能默认它们都有可用配套长读段。

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

### A2. 最小对照矩阵

| 对照 | 输入安排 | 要检查的内容 |
|---|---|---|
| 单株基线 | 各 FNA 分别运行 | 重现当前发现行为 |
| 联合发现 | 同一批 FNA 作为一个序列集合，包含跨记录比对 | 比较原生家族与成员；不能仅拼接独立运行结果 |
| 技术阳性 | 已知边界的合成重复，三个不同记录各一份、外侧独立 | 验证跨记录发现和坐标回溯；不是天然 IS 证据 |
| 覆盖/重复提交对照 | 同一序列记录的重复副本，明确标记人工来源 | 检查重复记录不会在报告中变成多个独立生物学位点 |
| 保守背景对照 | X 与较长侧翼共同相似的记录 | 观察软件如何处理保守宿主片段，保留结果 |
| 天然案例 | 有独立依据支持的自然重复家族 | 检查实际边界和侧翼表现，避免只验证合成例 |

若真实数据中未找到“每个来源一份、合并后多份”的天然案例，
明确报告；不能把合成阳性当成已经证明天然跨物种发现成功。

近缘基因组加入后，PILER 的 flank 在整个输入集合内可能不再独特。
不同 pool 组成分别运行并记录，不假设添加更多基因组必然提升召回。

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

边界优先报告已有四类状态及未解决原因。用于评价边界时，
采用 deposited annotation 的 candidate 与独立推断结果分开，
不能用注释替换后的精确一致证明从头边界推断准确。

### A4. 阶段 A 验收

- 所有合并记录能回溯到原始来源，序列与坐标核验通过。
- 单株和联合原生输出均完整保存，工具失败不被记录成零命中。
- 技术阳性给出明确的跨记录结果；若未恢复，记录是工具限制、输入问题还是实现错误。
- 报告逐株独有、联合独有、两者共有的可比家族/位点及匹配依据；
  不能只减去两个家族总数推算“新增”。
- 报告完整双侧 flank、截断 flank、边界未解决、tnp 完整覆盖及跨来源分布。
- 只有天然案例通过独立核查，才声称获得天然跨来源候选；
  species 归属不清的记录保留为 unknown。

## 6. 阶段 B：长读段从头发现工具评估

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

先找一个环境分离株的“原始长读段＋同一 BioSample 配套 assembly”，
再找一个混合微生物样本的配套数据。优先依据数据完整性，
而不是先宣称某个平台更好。

Nanopore 和 PacBio 分别记录：

- 平台、化学版本/basecaller 或 CCS 状态（可获得时）；
- SRA/ENA accession、BioSample、是否与 assembly 同源；
- read 数、总碱基、长度分布和质量信息；
- 原始数据、纠错数据、抽样数据的关系与哈希；
- 抽样方法、随机种子和样本内覆盖不均的限制。

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

对每个候选分别报告：

- 单条 read 同时跨越左 junction、candidate 和右 junction；
- 仅有分别支持左右 junction 的 reads；
- 单侧支持；
- candidate 内部命中但无连接证据；
- 没有可用覆盖；
- 多重比对、冲突连接或潜在嵌合，无法确定。

“两个 junction 分别有支持”与“同一分子跨越完整结构”分开。
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

执行顺序：A（联合 FNA）优先；B（长读段从头工具）先小规模；
有匹配 reads 即可开展 C（连接证据）。

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
