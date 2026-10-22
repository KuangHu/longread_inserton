# 联合发现的八株 assembly 候选

日期：2026-10-22  
状态：候选输入登记表；已完成首轮 reads 元数据预检，尚未完成下载、序列版本锁定或数据校验。  
关联：[跨物种与长读段评估计划](PLAN_CROSS_SPECIES_LONG_READS.md)。

此名单用于复现候选选择，不是已经可以直接运行的 package manifest。
GCA 是 assembly accession；CP 是单个复制子的 accession。
CP 条目在此保留已核实的 accession 主干，不猜测版本号；
必须在下载时解析为 accession.version 并冻结记录后才能运行。

| 顺序 | sample_id | organism / strain | 来源 | 定位 accession | 预期复制子 | 依据 |
|---|---|---|---|---|---|---|
| 1 | SaSB2B | Shewanella amazonensis SB2B | 巴西亚马孙河三角洲沉积物 | GCA_000015245.1 | 1 染色体 | [分类与 assembly](https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2019.01861/full)、[采样来源](https://pmc.ncbi.nlm.nih.gov/articles/PMC3147445/) |
| 2 | DdD110 | Dechloromonas denitrificans D110（NCBI strain: 110） | 德国巴伐利亚含水层地下水 | GCA_020510685.1 | 1 染色体 | [原始论文](https://doi.org/10.1128/mra.01020-21) |
| 3 | ShANA3 | Shewanella sp. ANA-3 | 美国 Eel Pond 半咸水河口木桩 | GCA_000203935.1 | 2 个复制子 | [BioProject](https://www.ncbi.nlm.nih.gov/bioproject/13905) |
| 4 | ScL1 | Sphingobium chlorophenolicum L-1 | 土壤 | CP002798; CP002799; CP002800 | 2 染色体＋1 质粒 | [原始论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC3318906/) |
| 5 | RoMoAcy1 | Rhodococcus opacus MoAcy1 / DSM 44186 | 德国 Tübingen 土壤 | CP082160; CP082159; CP082158 | 1 染色体＋2 质粒 | [原始论文](https://doi.org/10.1128/mra.00814-21) |
| 6 | PaBN5 | Paraburkholderia aromaticivorans BN5 | 韩国汽油污染土壤 | GCA_002278075.1 | 2 染色体＋6 质粒 | [BioProject](https://www.ncbi.nlm.nih.gov/bioproject/398462)、[原始论文](https://doi.org/10.1038/s41598-018-36165-x) |
| 7 | AmHOT1A3 | Alteromonas macleodii HOT1A3 | 夏威夷附近太平洋表层海水 | GCA_001578515.1 | 1 染色体＋1 质粒 | [完整组装论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC4781885/) |
| 8 | AaH17 | Alteromonas australica H17 | 澳大利亚 Port Phillip Bay 表层海水 | GCA_000730385.1 | 1 染色体 | [BioProject](https://www.ncbi.nlm.nih.gov/bioproject/246143) |

“预期复制子”来自公开记录，用于下载核对；不能以它代替实际文件盘点。
首轮长读段检索结果如下，详细查询范围、BioSample、平台、下载地址和来源见
[LONG_READ_AVAILABILITY.md](LONG_READ_AVAILABILITY.md)。
“未确认”仅表示在已查范围内未找到匹配记录，不等于证明不存在。

| sample_id | 长读段元数据状态 | run / 未完成项 |
|---|---|---|
| SaSB2B | 未确认匹配长读段 | 已见 RNA-seq 不作为 WGS 长读段 |
| DdD110 | PacBio 已匹配 | SRR16235513；SAMN19030608 |
| ShANA3 | 未确认匹配长读段 | 仍可用于阶段 A |
| ScL1 | 未确认匹配长读段 | assembly/BioSample 仍需解析 |
| RoMoAcy1 | PacBio 已匹配 | SRR15616890；SAMN12617337；CP→assembly 仍需解析 |
| PaBN5 | 未确认匹配长读段 | 同种 AR20-38 的 run 不能算作 BN5 |
| AmHOT1A3 | PacBio 已匹配 | SRR2084418；SAMN03835502 |
| AaH17 | 未确认匹配长读段 | 原组装论文使用 Illumina，不是长读段证据 |

据此保留原八株用于 A，C 优先 DdD110 → RoMoAcy1 → AmHOT1A3。
另列 PaAR2038（GCA_012689525.1 / SRR11947722）为独立 Nanopore B/C 对照；
它不进入本轮八株 pools，也不替换 PaBN5。以上是元数据匹配，不是已下载验证。
环境来源不表示从未实验室培养或传代。

## 下载后必须冻结

建立机器可读的 resolved manifest 与 source map，保存：

- 此表 sample_id、实际 assembly accession.version、BioSample、复制子 accession.version；
- assembly level、染色体/质粒标记、原始标题、长度、序列 SHA256；
- 精确下载 URL、下载时间、文件 SHA256、注释来源和版本；
- 原始文件及 pooled FASTA 的 ID 映射，annotation 与序列匹配状态；
- reads accession（若有）和是否匹配同一 BioSample；未知保持未知。

将所有复制子作为同一来源的记录保留；不只下载主染色体。
GCA assembly 中各 contig 的实际 accession 列表必须写入 resolved manifest。
CP 定位的两株应反查所属 assembly/BioSample；不能因 accession 相邻就假定同一来源。
解析失败、版本变化或复制子数不符时停止该条准备并报告，不静默换 assembly。

先与 cluster 的完整 panel.tsv 按菌株别名、BioSample 及 assembly 去重。
AltDE、UT26S、RHA1 是单独的既有回归对照，
明确排除在本轮 scale_2/4/8 及对应 spike-in pools 外。
ScL1–UT26S 同属 Sphingobium、RoMoAcy1–RHA1 同属 Rhodococcus；
亲缘关系不等于重复样本，但不能静默加入回归株而改变同属背景。
其他已跑菌株以 cluster panel 实际记录为准。

DdD110 的 resolved manifest 必须保留别名 D110（论文）、
110（NCBI strain）和 D.den110（BioSample sample name），关联
SAMN19030608。匹配 panel.tsv 时结合物种、BioSample 和 assembly，
不能把裸字符串“110”当作全局唯一标识。

## 预注册 pool 组成

资源阶梯用相同顺序的嵌套集合：

- scale_2：SaSB2B, DdD110
- scale_4：SaSB2B, DdD110, ShANA3, ScL1
- scale_8：上述全部八株

这些是预注册候选组成，尚不是已冻结的可运行输入。
首轮 reads 元数据核查在本次修订已完成；实际下载、版本、复制子、
注释与 source map 检查完成后才冻结 resolved manifest。
任何替换必须发生在查看 discovery 结果之前并另行记录。

同属共享背景预注册配对为 SaSB2B–ShANA3 与 AmHOT1A3–AaH17；
子 pool 仅在配对双方均存在时评估。所有回归株均不参与配对。

每个规模同时记录总碱基、记录数和数据相似性；规模之间的生物学组成也变化，
因此资源曲线仅是这一批数据的测量，不是算法复杂度的证明。
scale_2 不是 PILER 三拷贝技术阳性，不能据其零家族结果判断工具失败。

完整自然数据实验至少包括上述全体 pooled 运行和相应逐株基线。
这些候选不保证具有跨物种同一家族，也不保证每株仅一份 X。
天然“每株一份”的案例若不存在，按计划报告未观察到。

若有候选已在 cluster panel 中、下载失败或注释不足，保留该事实。
在任何 discovery 结果查看前提交新的明确 pool 清单及变更原因；
不自动替换菌株、不把不同组成的运行放进同一对照对。
