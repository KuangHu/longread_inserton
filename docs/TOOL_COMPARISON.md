# TLDR、DeBreak、Palidis、DIVE：真实 CCS pilot 对比

日期：2026-10-21。

## 本轮目标与可比范围

目标是在近缘非 lab-strain genomes 不足时，恢复 insertion 及尽可能准确的两端边界。真实 empty allele 是加分项；没有 empty 时可以保留 insertion + occupied upstream/downstream。**本轮没有把两种状态同时存在设为 INS 的核验条件。**

复用同一个真实海洋 CCS run SRR28395891 和同研究 MAG 池。冻结上次选出的10个MAG、235条contig、12,825,193 bp，原竞争比对参考为全部418个MAG。read-only输入为60,415条完整reads、280,787,461 bp；不是按插入结果筛选的reads。原始BAM中的MAPQ保留。

整个 archived run 的中位读长3,410 bp、N50 5,309 bp，header带CCS15、质量字符统一Q30，不能当作真实逐碱基准确率。它不代表所有现代长HiFi数据。MAG到源样本的关系未逐一证明；这次是同研究参考池pilot，不是严格同样本benchmark。未在真实Nanopore数据上重跑这四种方法。

比较范围500–50,000 bp；单独列出INS和DEL。长度是分析分层，不是IS判定。DIVE focal区间和TLDR focal区间按之前已知的1.95 kb案例选择，不能用局部结果估计全样本产量。Palidis仅完成原生IR模块及明确标记的探索性后处理，不与完整SV caller混排。

## 运行与数量

<!-- RUN_SUMMARY -->

| 方法/配置 | ≥500 bp INS候选 | 其中native PASS | 有显式INS序列 | ≥3条occupied read序列支持 | DEL候选 |
|---|---:|---:|---:|---:|---:|
| Sniffles2原候选基线 | 87 | 87 | 79 | 2 | 75 |
| DeBreak默认深度过滤 | 20 | 20 | 20 | 2 | 32 |
| DeBreak关闭深度上限（诊断） | 44 | 44 | 43 | 3 | 43 |
| TLDR无家族库/每簇最多20 reads | 20 | 10 | 20 | 0 | 不适用 |
| TLDR局部诊断，不能与全范围比产量 | 1 | 0 | 1 | 0 | 不适用 |

INS的3-read序列支持不依赖empty证据，也不代表IS鉴定。TLDR数量包括非PASS候选，需同时查看native filter。

| 其他路线 | 实际完成的范围 | 本次输出 |
|---|---|---|
| DIVE | 427条focal reads完成；60415条targets reads触及3 GiB限制失败 | 1个25 bp内部anchor；未恢复完整插入 |
| Palidis | 原生IR检测/切片/映射完成；完整Palidis4未运行 | 探索性配对得1个790 bp结构候选，非已确认IS |

长度分层（INS；TLDR含其原生非PASS记录）：

| 方法 | 500–999 bp | 1,000–1,999 bp | 2,000–4,999 bp | 5,000–50,000 bp |
|---|---:|---:|---:|---:|
| Sniffles2 | 18 | 28 | 34 | 7 |
| DeBreak_default | 9 | 7 | 4 | 0 |
| DeBreak_no_depth_cap | 18 | 16 | 9 | 1 |
| TLDR | 5 | 7 | 7 | 1 |

TLDR全范围重跑状态：`completed`，退出码0，耗时441.7秒；时间不是受控性能benchmark。
原生TLDR表共95条记录（含长度0及ShortIns等失败结果）；上表仅统计500–50,000 bp部分。

<!-- END_RUN_SUMMARY -->

所有“候选”均未做IS分类。这里没有真值集合，不能计算precision/recall，不能以call数量决定胜负。

## read核验标准

独立核验使用原有竞争比对的primary、MAPQ≥20 reads。INS要求完整CIGAR插入、100 bp双侧anchor分别达到95% identity，插入段与caller序列的全长identity≥95%，按read名去重。搜索范围±100 bp、长度差≤25%。允许相近断点分组，并明确报告坐标范围；至少3条occupied reads时记录序列支持，**不要求empty reads**。

这只是可复现的诊断检查：没有处理split-only支持、完整局部组装、菌株相位或跨MAG误归属，也没有证明单碱基边界。未通过不等于假阳性；通过不等于IS。DEL中“occupied支持”仅说明reads跨过参考中的候选区间，需另看支持gap的reads，不能直接视作边界已验证。

导出的500 bp参考flanks明确标记`reference_derived_context`；TLDR输出的flanks标为`occupied_consensus`。没有把拼接参考flanks冒充真实empty allele。局部reference-like reads也只记录支持计数，没有输出已完成phasing/polishing的empty allele。

## 能改变方法选择的具体发现

### TLDR无家族库模式在已知位点也有失败

首次全范围尝试在空共识FASTA进入重比对时终止；日志为`ValueError: file does not contain alignment data`，不把未写完的table当成零检出。之后使用原生`--max_cluster_size 50`和较高内存上限重跑，但进程在C11处理中断，未留下最终表或退出状态，具体中断原因未确定。最后使用`--max_cluster_size 20`且不启用额外detail BAM输出再试。簇大小限制会影响reads组成，三次配置和状态均独立保存。

对C11的[80000,110000)进行独立focal诊断，完成后输出4条记录，长度为1、4,043、0、0 bp；均非PASS。4,043 bp记录被`TSDOver100bp`过滤，未恢复已知的1,950 bp插入。局部区域是按已知案例选择的，这证明该具体配置在该位点失败，不能推算整个工具的recall，也不能外推所有HiFi或Nanopore数据。

最后的20-read簇上限全范围运行完成，原生表95条记录中20条≥500 bp、10条PASS；0条达到本轮3-read完整CIGAR序列核验条件。这里的0不是“全部假阳性”：该检查不包含split-only证据，且运行使用了资源限制参数。在已知1.95 kb位点，全范围输出的是1,459 bp非PASS记录，也未正确还原该等位型。

### DeBreak的默认全局深度过滤不适合直接套用这个MAG混合物

本轮以全局平均约22.49×得到默认`highcov=65`。即使使用`--no_genotype`，源码仍运行深度过滤；不同MAG覆盖度差异较大，真实候选可能被一起丢弃。

默认运行输出20 INS、32 DEL。诊断性设置`--maxcov 1000000000`后输出44 INS、43 DEL，其中43个INS有显式序列、1个仍为符号等位型。此设置只用来识别损失原因，不是推荐生产上限。与Sniffles2的87个INS候选相比，44个中31个按相同类型/contig、断点差≤100 bp、长度比≥0.75匹配。

### 1.95 kb：边界小幅偏移不应直接拆成不够支持的两个组

Prochlorococcus-A-LRMED-150 C11，参考插入位置约95,195（0-based interbase）。3条通过锚定的read插入长1,950、1,946、1,953 bp，位置95,193–95,195。Sniffles序列和DeBreak的1,950 bp输出均得到这3条reads的序列支持。

DeBreak已在中间结果检出它，却在默认深度过滤时丢失。关闭上限的诊断运行保留它，但原生标记是IMPRECISE。这里更合理的产物是带**2 bp边界范围**的序列候选，而不是要求所有read的CIGAR位置精确相同。局部38条reference-like reads提供额外证据，不是本次保留occupied候选的必要条件。

### 1.27 kb：共识校正有实际收益，但需要保留混合等位型疑问

Prochlorococcus C4，位置303,760。DeBreak给出1,274 bp共识；3条read的完整插入段与其identity为95.60%、95.60%、96.23%，都通过双侧anchor检查。旧Sniffles代表序列未达到相同的3-read核验条件。

这是共识校正的一个具体收益，不证明17条caller-support reads都属于同一等位型；原审计仍提示此位点有明显序列异质性。

### 3.1 kb复杂位点：PRECISE不能替代序列一致性检查

Prochlorococcus C11，约275,364。Sniffles给出的代表插入长3,109 bp；DeBreak共识长2,638 bp且标记PRECISE。后者有5条通过双侧锚定的相容gap reads，但**0条**完整插入段与该共识达到95% identity。旧审计也发现不同插入段两两identity约52%–54%。

因此，局部共识输出“更精确”不是充分证据。这里需要按序列拆分等位型，不能把邻近断点的不同序列合成一个可靠完整插入。

### DIVE找到相关内部anchor，但没有输出整段insertion

427条完整focal reads包含1,995,644 bp。提供`collections.Iterable`和`np.Inf`旧接口兼容别名后，原生DIVE运行完成，输出1个25 bp显著anchor，位于1,950 bp插入内部的offset 1,201（0-based）。它不是插入两端，也不是完整元素。

全部60,415条reads的尝试在3 GiB每进程地址空间上限下报MemoryError。这个失败说明当前配置没有完成全量分析，不是算法零检出，也不能据此断言更大机器必然失败。DIVE在本轮更适合作为参考依赖较小的候选信号来源，尚不能独立交付本项目需要的insertion+双侧flanks。

### Palidis的结构模块提供了一个不依赖empty的候选

旧Palidis仓库已宣布退休并指向Palidis4。检查实际源码发现Palidis4先用Prodigal/InterProScan选出transposase contigs，最终要求候选IR间存在已注释transposase；不能把完整新版称为不依赖已知家族的发现流程。

完整流程未运行：其固定版本InterProScan压缩包本身为7,973,063,848 bytes，已超过当时剩余磁盘空间；默认流程部分步骤配置16–64 GB内存，本机约10 GB。没有伪造注释来绕过这一点。

实际跑了pal-MEM原生unpaired输入、原生clip_reads和Bowtie2。120,134条IR片段记录中88,643条有比对，68,539条为完整零错配命中。探索性后处理将500–3000 bp范围内精确反向互补IR配对，得到1个结构区间：

- Actinomarina-LRMED-045 `contig_27`，区间[139003,139793)，790 bp。
- 两端28 bp精确反向互补，37条primary MAPQ≥20 reads跨越整个区间及两侧100 bp，identity≥95%。
- 已导出区间序列及两侧flanks；没有要求存在empty allele。
- 未鉴定为IS；可能是其他反向重复结构。严格映射和单个primary命中也可能漏检，**不能把“1个”当成完整Palidis的检出能力或召回率**。

这个案例表明“结构边界候选+occupied证据”路线有实现意义，但本次尚未展示大量可靠新IS。

## 对本项目的判断

**本轮没有证明换一个现成caller就能直接得到大量完整IS，也没有证明某个工具普遍最好。** DeBreak能够补充序列共识和候选，尤其需要避免全局深度规则错误处理MAG覆盖差异；Sniffles保留较宽候选集。两者的交集和差集都值得保留，后接统一的局部序列及边界核验。

结构IR路线可以独立补充参考里已经存在、没有empty对照的候选。DIVE适合提供内部可变信号，但需要后续局部序列重建；不能把短anchor数量当成完整insertion数量。TLDR的具体运行状态及可核验产出见上表。

下一步实现重点应是：候选合并后按插入序列拆分等位型、occupied侧两端重新比对、报告重复导致的边界区间；真实empty只作可选独立证据。当前比较脚本没有把这些尚未实现的步骤伪装成完成。

## 文件与复现

代码和本报告提交到GitHub。候选序列、逐read结果、原生工具输出和日志目前保留本地，尚未加入GitHub分支；结果数据上传被自动审核拦截，理由是此前只明确授权代码和Markdown。以下结果路径是本地路径说明，不是已发布下载链接。

- [复现说明和参数](../validation/tool_comparison/README.md)
- 按caller统计：`validation/tool_comparison/results/comparison/summary.json`
- 逐候选核验：`validation/tool_comparison/results/comparison/candidates.tsv`
- 逐read证据：`validation/tool_comparison/results/comparison/read_evidence.jsonl`
- 候选序列及flanks：`validation/tool_comparison/results/comparison/candidate_sequences.fasta`（含未通过核验记录）
- 完整工具版本与环境：`validation/tool_comparison/results/provenance.json`；主要源码commit也列于复现说明

所有失败尝试的日志和退出状态也保存在本地`validation/tool_comparison/results/runs/`。原pilot FASTQ/BAM、原正式accepted结果和主恢复算法未改动。

上游代码：[TLDR](https://github.com/adamewing/tldr)、[DeBreak](https://github.com/Maggi-Chen/DeBreak)、[Palidis4](https://github.com/blue-moon22/palidis4)、[DIVE](https://github.com/jordiabante/biodive)。具体commit以本次provenance文件为准。
