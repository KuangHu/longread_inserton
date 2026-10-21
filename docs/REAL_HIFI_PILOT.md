# 首个真实 CCS pilot：2026-10-21

**后续审核：按完整插入元件目标统计，≥1 kb 且双状态支持的结果仅 3 个，≥2 kb 为 0；20 个结果不等于 20 个可用 IS。见 [长插入审计](HIFI_PILOT_AUDIT.md)。**

**已跑通下载、竞争比对、覆盖度筛选、INS/DEL calling 和序列恢复。恢复 20 个序列表示，其中 3 个通过当前严格标准，17 个存在重复导致的等价断点。**

但这次尚未完成“严格同一样本 assembly/MAG”的验证：参考是作者公开的同研究 MAG 池，单个 MAG 到源样本的对应关系未得到完整证明。因此它是一个真实数据流程 pilot，不能当作严格同样本 benchmark，更不是新移动元件的发现结论。

## 数据与实际规模

| 项目 | 实测结果 |
|---|---:|
| 自然样本 | 地中海海水 Med-OCT2021-15m |
| CCS run | SRR28395891 |
| 完整压缩 FASTQ | 547,419,476 bytes |
| reads | 418,396 |
| 总碱基 | 1,764,906,237 bp |
| 中位读长 / N50 | 3,410 / 5,309 bp |
| 参考 MAG 池 | 418 MAG；17,441 contigs；594,150,817 bp |
| mapped primary reads | 194,311（46.44%） |
| MAPQ ≥20 的 mapped primary reads | 178,918（42.76%） |
| 全池平均有效深度 ≥5× / ≥10× 的 MAG | 30 / 9 |
| 本轮真正进行 calling 的范围 | 10 MAG，235 contigs |

这是**完整 archived run**，不是之前的 500-read 下载测试。全文件大小与 ENA MD5 均验证通过；整个 FASTQ 已读到 EOF。所有 read header 含作者的 CCS15 标记，质量字符全部固定 Q30，不能据此推断真实逐碱基准确率。论文 Table S1 的分析数据量与归档 run 不同，本轮不宣称复现原文过滤或检出数量。

参考：Zenodo 14650439 的 MAG-LRMED.tar，完整 MD5 校验通过。原包有 419 个 FASTA；按论文 Table S2 保留 418 个。另有 14 个 MAG 的 archive 长度或 contig 数与表格不一致：它们保留为竞争参考，但不作为本轮 calling 目标。覆盖度分母以实际 FASTA 为准。

## 目标选择

先将全量 reads 比对到完整 MAG 池，再选目标；不先将 reads 强行分配给单个 MAG。目标依据：平均有效深度 ≥5×、覆盖至少 70% 的 MAG 碱基、论文 completeness ≥75%、contamination ≤5%，且 archive/table 无上述差异。在符合条件的 12 个 MAG 中按深度选择前 10 个；没有依据 SV 数量挑目标。

| MAG | 平均有效深度 |
|---|---:|
| Pelagibacter-LRMED-227 | 74.1× |
| Prochlorococcus-A-LRMED-150 | 40.3× |
| Actinomarina-LRMED-045 | 22.5× |
| Parasynechococcus-LRMED-158 | 9.1× |
| SAR86-LRMED-384 | 8.4× |
| Alphaproteobacteria-LRMED-289 | 7.9× |
| Marinisoma-LRMED-168 | 7.8× |
| Pelagibacter-LRMED-234 | 7.2× |
| Puniceispirillales-LRMED-245 | 6.5× |
| Pelagibacter-LRMED-233 | 6.4× |

深度使用 primary、非 duplicate、非 QC-failed、MAPQ ≥20 的对齐 M/=/X 碱基；不把 deletion/skip 算成覆盖。MAG 平均深度不能代替具体断点的支持数量。

## Calling 与恢复结果

Sniffles 共报告 312 个 SV：138 INS、138 DEL、35 BND、1 DUP。本流程只尝试 INS/DEL，另外 36 个保留为 unsupported。

| 恢复结果 | 数量 | 含义 |
|---|---:|---|
| read_supported_pair | 3 | 两种结构都有足够 read 支持，flanks 完整，未检出等价重复断点 |
| ambiguous | 17 | 恢复了序列，但重复允许多个等价断点；其中 13 个两种结构均有足够支持 |
| unresolved INS/DEL | 256 | 暂不输出猜测序列 |
| unsupported BND/DUP | 36 | 不属于本流程目标 |

因此，20 个恢复结果中，**16 个具有双状态 read 支持**；只有其中 3 个同时满足当前“不存在重复等价断点”的严格分类。其余 4 个 ambiguous 还缺少足够 occupied reads，不能当作双状态验证。

严格通过的 3 个事件均为参考 occupied 的 DEL 路线：插入区段从 assembly 取出，empty reads 跨越删除接头；DEL 是相对参考的描述，不代表已确定进化上发生了删除。

| MAG / 原始 contig | 参考插入区间 [start0,end0) | insertion 长度 | occupied / empty reads |
|---|---|---:|---:|
| Marinisoma-LRMED-168 / C118 | [92114,92253) | 139 bp | 3 / 4 |
| Pelagibacter-LRMED-227 / C250 | [73443,73820) | 377 bp | 17 / 6 |
| Pelagibacter-LRMED-227 / C3161 | [25556,25968) | 412 bp | 10 / 4 |

这三个结果的 empty-site 均为 2,000 bp（上、下游各 1,000 bp），内部插入坐标为 1,000。完整 contig ID、read 名及序列见结果文件。

17 个 ambiguous 结果含 3 个 read-derived INS consensus，长度分别为 108、262、1,161 bp；其他为 reference-derived DEL 区段。它们不应仅因等价断点而被视为“没有恢复出序列”，但不能宣称唯一的精确切点。

256 个未恢复 INS/DEL 的主要原因：225 个没有符合要求的完整 CIGAR 事件和双侧锚定；29 个同一精确断点支持不足；2 个插入序列不一致。本轮保持原阈值，没有为了增加产量降低条件。

## 本地运行与限制

- 环境：8 GiB 内存上限；使用 4 threads。
- minimap2：249.1 秒，峰值 RSS 4.403 GB；排序后 BAM 约 478 MiB。
- 235 contigs 的 Sniffles calling：12.83 秒。
- 最初连续下载在 413,138,944 bytes 时超时；通过 HTTP Range 补齐，再校验全文件 MD5。
- 最初默认大 batch 比对被主动中断，改用 `-K 50M`。未发生记录到的 OOM kill；比对阈值保持不变。
- 全池 broad calling 被中断，保留已完成的完整竞争比对，改为上述 coverage-based 目标集。没有完成全 418 MAG 的 SV 扫描。
- source-sample 配对仍是核心限制：这是同研究参考池，而非已验证每个 MAG 与 reads 严格来自同一个采样样本。
- empty-site 是 assembly flanks 的表示，局部锚定为 100 bp；没有对整段 2 kb flanks 做菌株分型和 read consensus 校正。
- 未做独立真值验证，不能由这些结果估计 precision/recall，也不能据此判定移动性、功能或同菌株来源。

## 文件

运行脚本和完整参数见 [pilot 目录](../validation/real_hifi_pilot/README.md)。

- [总体统计](../validation/real_hifi_pilot/results/pilot_summary.json)
- [全部候选、坐标与 MAG 归属](../validation/real_hifi_pilot/results/events_with_mag.tsv)
- [全参考池覆盖度](../validation/real_hifi_pilot/results/mag_coverage.tsv)
- [严格结果 insertion FASTA](../validation/real_hifi_pilot/results/accepted/insertion.fasta)
- [严格结果 empty-site FASTA](../validation/real_hifi_pilot/results/accepted/empty_site.fasta)
- [歧义候选 insertion FASTA](../validation/real_hifi_pilot/results/candidates/insertion.fasta)
- [20 个恢复结果的详细记录](../validation/real_hifi_pilot/results/recovered_records.jsonl)

参考来源：<https://doi.org/10.1186/s40793-025-00773-8>；<https://zenodo.org/records/14650439>；<https://mgegroup.com.es/datasets/>。
