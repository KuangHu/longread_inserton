# 长读段证据实验结果：不同数据类型给出哪类 insertion callout

日期：2026-10-23  
状态：首轮实验已完成（四株、各一个公开长读段 run）。阈值均在查看结果前固定；
下文“已知局限”中的修订尚未执行，若执行将作为 v2 另行报告，不覆盖本轮。  
关联：[计划](PLAN_CROSS_SPECIES_LONG_READS.md)、[reads 可用性](LONG_READ_AVAILABILITY.md)、
代码与原始汇总 [experiments/read_callout/](../experiments/read_callout/)。

## 1. 设计

每株三层证据，均针对同一份完整 assembly：

| 层 | 输入 | 输出 |
|---|---|---|
| assembly-only | RefSeq GCF FNA＋PGAP GenBank 注释 | multisite-insertion `--subset all`：原生家族、成员、边界状态、flank |
| reads → 已知候选 | reads 以 minimap2 比对到完整 assembly | 每个候选按**分子**计的 full span / junction / empty allele 等类别 |
| reads-only | 同一比对中的 split read 与大 CIGAR 插入 | 不在 assembly 中的插入：可归入候选家族者、无法归类者 |

另以确定性分子抽样（种子 11/29/47）在 5/10/20/40/80× 和全深度重算，得到深度曲线。
GCF 与 GCA 的序列集合逐条 SHA256 一致；用 GCF 是为了四株使用同一 PGAP 注释流程
（AR20-38 的 GCA 无注释）。

| sample | 平台 / 仪器 | run | reads N50 / 中位 | 分子深度 | 分子分组 |
|---|---|---|---|---:|---|
| DdD110 | PacBio Sequel CLR | SRR16235513 | 7.9 / 3.4 kb | 119× | ZMW（read 名） |
| RoMoAcy1 | PacBio RS CLR | SRR15616890 | 18.2 / 11.4 kb | 177× | 每 spot 一条，按分子计 |
| AmHOT1A3 | PacBio RS II CLR | SRR2084418 | 7.4 / 1.8 kb | 265× | SRA spot 视作 ZMW（假定） |
| PaAR2038 | Nanopore PromethION | SRR11947722 | 19.7 / 13.3 kb | 195× | 每条 read 一个分子 |

四个 run 均通过 ENA 字节数、MD5、碱基数核对。DdD110 约 21% 记录内含 SMRTbell adapter，
比对前在 adapter 处切开（共切 43,168 条记录）。分子深度 = 各分子最长片段之和 / 基因组长度。

## 2. 结果

| | DdD110 | RoMoAcy1 | AmHOT1A3 | PaAR2038 |
|---|---:|---:|---:|---:|
| assembly 原生家族 / 候选成员 | 19 / 164 | 94 / 646 | 30 / 188 | 149 / 762 |
| 其中 ≥400 bp 候选 | 64 | 360 | 134 | 526 |
| tnp-supported 家族成员 | 0 | 31 | 81 | 279 |
| 全深度 full span（单分子跨越 flank–X–flank） | 164 | 644 | 188 | 761 |
| edge_limited | 0 | 2 | 0 | 0 |
| empty allele（≥2 分子） | 0 | 0 | 0 | 1* |
| 非参考位点：候选家族插入（≥2 分子；双侧） | 1 (1) | 2 (2) | 3 (1) | 1 (1) |
| 非参考位点：无法归类插入（≥2 分子） | 99 | 109 | 60 | 3 |

\* 44 bp 的 RepeatScout 片段；2/197 个分子带 ≥22 bp 缺失，符合 Nanopore 缺失误差，不作多态解释。

### 深度曲线：full span 占全部候选的比例（三种子平均）

| 深度 | DdD110 | RoMoAcy1 | AmHOT1A3 | PaAR2038 |
|---|---:|---:|---:|---:|
| 5× | 97.6% | 85.3% | 66.3% | 98.6% |
| 10× | 99.6% | 90.7% | 90.4% | 99.9% |
| 20× | 100% | 96.0% | 96.8% | 100% |
| 40× | 100% | 97.3% | 99.3% | 100% |
| 80× | 100% | 99.2% | 100% | 100% |

### 深度曲线：无法归类插入位点数（≥2 分子）

| 深度 | DdD110 | RoMoAcy1 | AmHOT1A3 | PaAR2038 |
|---|---:|---:|---:|---:|
| 20× | 3.0 | 1.0 | 0.3 | 0 |
| 40× | 12.3 | 6.7 | 1.0 | 0 |
| 80× | 46.0 | 25.7 | 2.7 | 0.3 |
| 全深度 | 99 | 109 | 60 | 3 |

## 3. 解释：数据类型 → callout 类型

1. **Assembly-only** 给出候选家族、成员、边界状态和 flank，即“哪里可能有 X”。
   它不能说明该拷贝是否在群体中固定，也看不到组装折叠或遗漏的拷贝。
2. **长读段回贴同一 assembly** 给出每个候选的结构连接证据。
   四种数据在全深度都支持 99–100% 候选的完整 flank–X–flank 结构；
   这些 assembly 很可能就是由同一批 reads 组装，因此主要是一致性检查，不是独立验证。
   所需深度主要取决于 read/subread 长度和重复复杂度，而非平台名：
   长分子（Nanopore、Sequel，N50 ≥ 8 kb）5× 已近 98–99%；
   AmHOT1A3 subread 中位仅 1.8 kb，5× 只有 66%，需约 20×；
   RoMoAcy1 reads 虽长，但多拷贝短重复多，同样约需 20× 达 96%。
   三株均未见可信的 empty allele，即在这些位点未观察到群体内多态。
3. **reads-only 的非参考插入** 是长读段相对 assembly 的增量，但强烈依赖读段准确度。
   Nanopore 背景很低（195× 时 3 个无法归类位点）；
   PacBio CLR 的无法归类位点随深度近线性增加（DdD110：10× 为 0，全深度 99），
   符合 CLR 插入误差/嵌合背景，而非真实插入。以“≥2 分子”这一绝对阈值，
   CLR 数据的非参考插入不可直接采信。
4. 唯一值得人工核查的线索：PaAR2038 在 NZ_CP051514.1:2,173,086 的
   IS5 家族（PaAR2038_RS_R1，1,171 bp，tnp-supported）非参考插入，
   插入长度 1,176 bp、双侧连接，但仅 2/≈195 分子；可能为低频亚群插入或嵌合，未核查 TSD。

本实验中四个平台来自四个不同物种，结果反映“这类数据能给出什么”，不是平台基准比较。
未测试：reads-only 从头重复发现（RepLong/LongRepMarker，计划 B0）和 HiFi 数据（无匹配 run）。

## 4. 已知局限（预注册阈值在本数据上的问题）

- 非参考位点用固定 ≥2 分子；在 100–265× 下过宽。需加入 allele fraction
  （支持分子 / (支持分子＋跨越参考位置的分子)）。
- empty allele 的缺失阈值为候选长度的 50%；短候选（如 44 bp）只需 22 bp 缺失即触发。
- 家族归属允许“插入包含 ≥80% 成员”；<400 bp 的 RS 短家族会被大量无关插入包含，
  RoMoAcy1 两个家族位点即属此类。应对 <400 bp 候选单列或设最短长度。
- AmHOT1A3 的分子分组假定 SRA spot＝ZMW，未用 SRA 元数据单独确认。
- split-read 路线在模拟中未被触发（minimap2 将 ≤7.2 kb 插入均表示为 CIGAR I）；
  真实数据中该路线有输出，但未单独验证。

以上修订若执行，作为 v2 另行编号；本轮结果和阈值保持不变。

## 5. 运行记录

- 集群：LBL Lawrencium，account pc_rubinlab，partition lr6；工具版本见各样本目录。
- PaAR2038 第一次下载 MD5 不符（gzip CRC 错误），保留审计记录后重新下载通过。
- AmHOT1A3 首次核对因 SRA spot 数少 23 被判 mismatch；MD5 与碱基数完全一致，
  缺少的 spot 不含碱基。核对规则改为：仅在碱基数一致时允许 FASTQ 中缺少空 spot。
- 节点 n0348.lr6 与 n0039.lr6 出现文件系统读写异常，重跑时排除。
- 每样本资源：assembly callout 13–42 min；reads 步骤 16–31 min（16 CPU）。
- 原始 reads、PAF 与中间文件留在 `$PROJ=/global/scratch/users/kh36969/longread_callout`，
  不入库；入库结果中该前缀已替换为 `$PROJ`。
