# 联合发现的八株 assembly 候选

日期：2026-10-22  
状态：候选输入登记表；尚未完成下载、序列版本锁定或配套 reads 核验。  
关联：[跨物种与长读段评估计划](PLAN_CROSS_SPECIES_LONG_READS.md)。

此名单用于复现候选选择，不是已经可以直接运行的 package manifest。
GCA 是 assembly accession；CP 是单个复制子的 accession。
CP 条目在此保留已核实的 accession 主干，不猜测版本号；
必须在下载时解析为 accession.version 并冻结记录后才能运行。

| 顺序 | sample_id | organism / strain | 来源 | 定位 accession | 预期复制子 | 依据 |
|---|---|---|---|---|---|---|
| 1 | SaSB2B | Shewanella amazonensis SB2B | 巴西亚马孙河三角洲沉积物 | GCA_000015245.1 | 1 染色体 | [分类与 assembly](https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2019.01861/full)、[采样来源](https://pmc.ncbi.nlm.nih.gov/articles/PMC3147445/) |
| 2 | DdD110 | Dechloromonas denitrificans D110 | 德国巴伐利亚含水层地下水 | GCA_020510685.1 | 1 染色体 | [原始论文](https://doi.org/10.1128/mra.01020-21) |
| 3 | ShANA3 | Shewanella sp. ANA-3 | 美国 Eel Pond 半咸水河口木桩 | GCA_000203935.1 | 2 个复制子 | [BioProject](https://www.ncbi.nlm.nih.gov/bioproject/13905) |
| 4 | ScL1 | Sphingobium chlorophenolicum L-1 | 土壤 | CP002798; CP002799; CP002800 | 2 染色体＋1 质粒 | [原始论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC3318906/) |
| 5 | RoMoAcy1 | Rhodococcus opacus MoAcy1 / DSM 44186 | 德国 Tübingen 土壤 | CP082160; CP082159; CP082158 | 1 染色体＋2 质粒 | [原始论文](https://doi.org/10.1128/mra.00814-21) |
| 6 | PaBN5 | Paraburkholderia aromaticivorans BN5 | 韩国汽油污染土壤 | GCA_002278075.1 | 2 染色体＋6 质粒 | [BioProject](https://www.ncbi.nlm.nih.gov/bioproject/398462)、[原始论文](https://doi.org/10.1038/s41598-018-36165-x) |
| 7 | AmHOT1A3 | Alteromonas macleodii HOT1A3 | 夏威夷附近太平洋表层海水 | GCA_001578515.1 | 1 染色体＋1 质粒 | [完整组装论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC4781885/) |
| 8 | AaH17 | Alteromonas australica H17 | 澳大利亚 Port Phillip Bay 表层海水 | GCA_000730385.1 | 1 染色体 | [BioProject](https://www.ncbi.nlm.nih.gov/bioproject/246143) |

“预期复制子”来自公开记录，用于下载核对；不能以它代替实际文件盘点。
所有行的配套 Nanopore/PacBio reads 状态均为 **not_checked**。
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
AltDE、UT26S、RHA1 是单独的既有回归对照，不混计为八株新候选；
其他已跑菌株以 cluster panel 实际记录为准。

## 预注册 pool 组成

资源阶梯用相同顺序的嵌套集合：

- scale_2：SaSB2B, DdD110
- scale_4：SaSB2B, DdD110, ShANA3, ScL1
- scale_8：上述全部八株

每个规模同时记录总碱基、记录数和数据相似性；规模之间的生物学组成也变化，
因此资源曲线仅是这一批数据的测量，不是算法复杂度的证明。
scale_2 不是 PILER 三拷贝技术阳性，不能据其零家族结果判断工具失败。

完整自然数据实验至少包括上述全体 pooled 运行和相应逐株基线。
这些候选不保证具有跨物种同一家族，也不保证每株仅一份 X。
天然“每株一份”的案例若不存在，按计划报告未观察到。

若有候选已在 cluster panel 中、下载失败或注释不足，保留该事实。
在任何 discovery 结果查看前提交新的明确 pool 清单及变更原因；
不自动替换菌株、不把不同组成的运行放进同一对照对。
