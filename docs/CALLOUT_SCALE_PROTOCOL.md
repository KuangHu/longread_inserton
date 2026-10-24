# Insertion call-out：扩容前固定规范 v1

日期：2026-10-23。状态：**输入、发现与报告规范已固定；固定五株验收在相邻竞争比对修复（nearby_competition_v2）后已重算并通过，`ready_for_10_sample_batch = true`**，见 [callout_scale_v1 结果](../experiments/callout_scale_v1/README.md)。
优先级：先把单株 `[upstream][candidate X][downstream]` 做成可批量复现的输出。
跨株同源位置、empty-site reads、低频非参考插入、跨物种 pooled discovery 均为后续独立模块，
不作为这一轮单株 call-out 的依赖。本规范取代此前“先扩展跨株比较”的执行优先级。

## 1. 要交付的是什么

对每个样本，输出真实 assembly 上每条候选 X 及其两侧各 500 bp，带来源、方向、
native/inferred/candidate 三套坐标、边界状态及本株 reads 支持。
同家族多拷贝是候选集合；不同坐标、不同侧翼序列、不同同源背景、历史插入事件四者不混用。
**不要求空位点、不要求 tnp、不要求 TIR/TSD。**
高 read 覆盖证明所提取结构与 reads 一致，不能独自证明元件精确起止或转座活性。

保留三层，避免一种 PASS 覆盖所有问题：

| 层 | 含义 | 用途 |
|---|---|---|
| 全部原生记录与导出候选 | 两软件原样输出；短、边界未解和无 tnp 均保留 | 可追溯的发现全集 |
| 优先检查候选表 | X ≥400 bp、两侧完整 500 bp、triplet 为 A/C/G/T，边界来源为 deposited 或 provisional median | 控制人工检查工作量，仍不是精确元件认证 |
| read 支持层 | 对以上及其它候选分别附 v1 或审核版支持，永远标 classifier 版本 | 判断结构证据；不反向删除发现 |

400 bp 是固定的报告分层，不能成为软件输出的淘汰器。其它层和失败原因必须同时交付。
N 分别列 raw_member_records、exact_coordinate_groups 和尚待解析的重叠记录。
目前 package 不会自动解析跨工具非完全重合边界，不宣称 exact groups 是独立 locus/event 数。
同一物理位置两软件的支持不能相加；没有逐 read 身份集合就不计算合并支持数。

## 2. 固定 datatype 与样本资格

**首批：培养分离株的 genomic assembly FNA＋同一 BioSample 的 genomic Nanopore FASTQ。**
每株染色体与关联质粒全部保留，每株一个发现任务。不同菌株的 FASTA/reads 不合池。
数据库入口为 NCBI Assembly/Datasets＋BioSample 关联的 ENA/SRA reads；
同一 BioProject、菌株名字相近或同种，均不能替代 BioSample/文库匹配。

| 项目 | 本轮固定处理 |
|---|---|
| 连续性 | Complete Genome 或高连续性 draft：FASTA 记录 ≤10、最长记录 ≥3 Mb。另列 N50、N/gap、总长、记录清单；这些阈值不证明完整度 |
| 局部结构 | 候选必须逐位点报告侧翼长度、gap/模糊碱基及边缘限制，不能只凭全基因组等级判断 |
| 来源 | 环境/植物相关分离株优先；富集培养可用并保留历史；工程化、来源未知、PCR 文库单独标记，不冒充野外独立重复 |
| 原始 reads | 保存下载校验、全文件长度分布、bases、read/spot 单位、原始命名与处理历史 |
| 深度 | 记录总 reads bases / assembly bp 和有效局部跨越；不设未经校准的全局深度硬门槛，LjRoot223 的约 14×是压力样本 |
| 注释 | 优先同序列版本 GenBank/EMBL；GCA/GCF 必须核对所有记录序列；不强行把旧 GCF 注释套到新 GCA。无注释照常 --subset all，tnp 状态写 unavailable |
| 读段过滤 | 新标准不额外删 <500 bp reads；历史过滤保留原始日志及删除 records/bases，重算用新目录，不覆盖旧结果 |
| PCR | SMT4a 用作 PCR 文库压力样本；其 read 数不认证原始分子数，低频等位型不属于本轮主输出 |
| 去冗余 | 可按整体 ANI 选代表减少主批重复，但保留各样本身份；不能据此断言插入结构相同或借用另一株 reads |

HiFi、PacBio CLR、MAG/混样和 reads-only 从头发现暂不混入首批。
这是统一验收范围，不是宣称 Nanopore 在所有物种上优于其它平台。
注释和 BLASTX 参考蛋白来自当前 package 所提供样本；固定一株一任务，归档实际蛋白库和哈希，
不把不同注释覆盖造成的 tnp 数差异当作发现能力差异。

## 3. 固定 pipeline

1. **Fetch/verify**：FNA、注释、reads 分开校验；保存 accession.version、BioSample、run、
   URL、档案 MD5/bytes、实际 SHA256、文件记录数、过滤历史及文库来源。
2. **Native discovery**：PALS→PILER-DF 与维护版 RepeatScout 原生默认；逐株运行
   `multisite-insertion run --methods PILER RS --subset all`。
   固定 source commit、所有工具版本/二进制 hash；记录 RepeatScout 两阶段实际 l 与输入字节数。
   不加 TIR、TSD、tnp、拷贝数、复杂度和短片段发现过滤。
3. **Boundary/export**：沿用 package 的 ±2000 bp 同源范围估计和现有 annotation 关联规则；
   保存 native、inferred、candidate；导出 X、up500、down500、triplet 以及 genomic-forward 版本。
   五类 boundary_source 保持不变；窗口截断、未解和 origin/contig 边缘受限不能被 read 支持升级。
4. **Assembly QC**：`multisite-insertion validate --run ...`，逐条校验来源、切片、方向、拼接、
   原生文件 hash、全部成员保留。失败时停止该样本，不能当成零候选。
5. **Filled-site read check**：本株 reads 映射到本株全部 assembly 记录，保存全部实际命令与 reference hash。
   固定当前 ONT 基线：`minimap2 -x map-ont -c --secondary=yes -N 50 -p 0.8 -I 8G`。
   `-p 0.8` 本来就是默认；`-N 50` 仍是上限。达到上限的候选必须重新检查竞争比对，不能自动判唯一。
   旧 v1 原样保留，审核版另存目录和逐 read 证据；无需调用非参考插入/空等位型模块才能完成 filled-site 检查。
6. **Report/check**：运行下述 bundle checker，输出每株方法/家族/长度/边界/read 支持分层，
   每个统计附分母；独立汇总资源和失败，最后才决定能否扩大批次。

现有 `03_reads.sbatch` 还包含非参考插入步骤，且会删除旧 evidence 目录；
它是历史实验脚本，不是本规范的安全重跑入口。集群执行时必须使用新的 sample 输出目录，
并拆出只做 filled-site 的 wrapper 后再扩大任务数组；不要直接覆盖 round-1。

## 4. Read 规则的冻结边界

当前 v1 判据主要是 primary/MAPQ＋参考跨度＋大 indel 排除。
它可用于探索性支持和回归，但并未完成计划中的所有竞争位置、CIGAR 和来源关系审核。
不得仅把 v1 输出列名改成 audited。

审核版采用已登记的初始规则：每个 read 单独在一条连续比对中跨过 candidate 和两侧各 100 bp；
M/=/X 参考覆盖在 candidate 与每个 100 bp anchor 中分别 ≥90%；candidate 内 I+D 合计 ≤10% X 长度；
至少两个身份可区分的支持单位。明确 CIGAR N/D、边界插入归属、反链 query 坐标和冲突处理。
两条各跨一侧的 reads 不能拼成一条 full-span；不能只看 primary/MAPQ 忽略竞争定位。
缺失原始 parent/split/duplex 或 PCR 分子信息时标 independence_unverified，保留 read-level 支持。
同批用于组装的 reads 标 assembly_consistency，不称独立实验验证。

**审核版已实现并在相邻竞争定位修复后归档五株重算，classifier revision 为 `nearby_competition_v2`**。
五株分类均未因该修复改变；R-93128 只有 4 条短 tandem-repeat records 的支持 read 数减少 1。
empty-site reads 实现不阻塞上述 filled-site 验收。

## 5. Test run 固定五株

锁定列表见 [test_panel.tsv](../experiments/callout_scale_v1/test_panel.tsv)。
同种比较结果只是辅助观察，不要求其阳性才能通过单株 call-out 测试。

| 样本 | 角色 | 核对重点 |
|---|---|---|
| AR20-38 | 已归档阳性回归 | RS_R1：11 条、1171–1173 bp、v1 full 145–181；坐标/边界/导出逐成员对照 |
| SMT4a | PCR 文库、重叠记录压力 | 完整 triplet、分子独立性标记、4 条重叠记录不得算 4 个事件 |
| LjRoot223 | 较低深度、历史短 read 过滤 | 支持不足与不存在分开；过滤前后分母完整保存 |
| R-93128 | 高连续性 draft、ONT-only assembly | 多记录参考、同批 reads 一致性、重复组来源保留 |
| R-91497 | 较高深度、另一组菌株 | 低候选产量也必须完整报告；不为追阳性放宽规则 |

固定顺序：AR20-38 → SMT4a → LjRoot223 → R-93128 → R-91497。
本次 GitHub 基准 ee229861 只有 round-1 samples.tsv 和文库审计，尚无其结果/比较脚本。
用户汇报的 terricola 12 / fungorum 2 个株间差异位置记为 **reported_pending_archive**，
不能填进机器回归 expected 文件当作已经独立核验的输出。
先把集群的原始小结果、curl 修复、完整命令和退出码归档，再做审核版对照。

## 6. Check result 与扩大批次的验收

### T0：软件与计算对照

- 包测试全通过；真实环境 doctor 依赖可用，版本/二进制 hash 齐全。
- 合成 demo 用真实 PILER/RS 运行并 validate，检查已知重复及切片。它只验证工具连通和坐标，
  不代表基因组规模准确率；原生召回和边界偏差保留，不能改输出凑精确边界。
- Read 审核版需有明确真值：完整跨越、删除 X、大量小 indel 累积、anchor 缺失、CIGAR N、
  等分竞争位置、secondary 上限、反链、重复派生 read、两条单侧 read、边缘截断。
  挑战例不能被误报为审核通过；故障退出与合法零命中必须可区分。

### T1：固定五株结果复核

- 全部 5 株状态 complete、全体原生记录可追踪、坐标/方向/序列/hash 检查 PASS；不能丢失失败样本。
- AR20-38 冻结输入的 native/导出逐条保持一致；审核版 read 分类允许变化，但附逐成员差异及原因。
- 各株读段校验、文库审计和报告分母齐全；方法、长度、tnp、边界与 read 层分别统计。
- 人工核查 AR20-38 的 11 个 IS5 成员；其余每株按固定 member ID 顺序取前 5 个优先候选，
  不足 5 个则全部；另各取首个 short、edge/censored、non-tnp、conflict/ambiguous（若存在）。
  检查 assembly 切片、peer 边界依据与 read 实际 CIGAR/竞争位置，保存 member ID 和判定表。
  这是一致性审计，不估计总体 FDR；发现系统性错误必须修复并重跑受影响样本。
- 零候选、零新增背景及无 tnp 均不算技术失败。没有强行要求 N=50 或新增 ≥1 的生物学通过线。

### T2：10 株受控批次，再决定 50 株

T0/T1 后按来源和数据资格预先选 10 株，不按候选产量选；新样本先每株独立作业。
初始最多 2 个并发；沿用资源上限：发现 2 CPU / 32 GiB / 2 h，reads 16 CPU / 64 GiB / 6 h。
记录 wall time、CPU time、peak RSS、临时磁盘、退出码及每个阶段输出 hash。
先用五株实测校准资源；OOM/超时记录原因后显式调整，不能无记录跳过。
10/10 都有终态，成功样本全部通过 QC，失败已定位且无共享流程缺陷后才扩到 50；
若同类系统错误影响 ≥2 株，暂停扩容修复。生物学阴性无需换样。
预估 CPU/存储用 pilot 每株实测分布，不从候选数推算；raw reads 与 PAF 留在集群。

## 7. 现有可运行检查命令

```bash
python -m unittest discover -s tests -v
multisite-insertion doctor > doctor.json
multisite-insertion run --manifest /path/to/one_sample.json \
  --out /path/to/new_sample/callout --methods PILER RS --subset all --threads 2
multisite-insertion validate --run /path/to/new_sample/callout
```

本株 v1 evidence 已生成后，运行本次新增 checker：

```bash
python scripts/check_callout_bundle.py \
  --sample-dir /path/to/new_sample \
  --reads-verify /path/to/reads/RUN/verify.json \
  --out /path/to/new_sample/qc/bundle_v1.json
```

Live 模式要求 `callout/` 是完整 package run，`evidence/` 包含 v1 的
member_read_support.tsv 和 summary.json。它重新调用 validate 并检查表间 ID/坐标/计数、
下载校验状态和样本身份；不重新验证原始 FASTQ 或重新比对 reads。
它不会调用原生发现，也不会删除、过滤或重写发现输出。
`scale_release_approved=false` 是有意的：这个检查器不能替代 T0/T1 和 read 审核。

仓库小结果的可复现检查（不需要原生工具或大文件）：

```bash
python scripts/check_callout_bundle.py --snapshot \
  --sample-dir experiments/read_callout/results/PaAR2038 \
  --reads-verify experiments/read_callout/results/PaAR2038/reads_verify.json \
  --out /tmp/PaAR2038_snapshot_new.json
```

Snapshot 模式只查小表和归档状态，无法复验 FASTA 切片/原始 read/完整 hash 链；
不能把 snapshot PASS 写成重跑发现或测序验证 PASS。输出已存在时拒绝覆盖。

## 8. 本次交付与未完成项

已交付：本规范、固定五株 test_panel、bundle checker、破坏输入/错配回归测试，
以及[四株旧归档检查结果](../experiments/callout_scale_v1/snapshot_checks.json)。
最初本地检查包含合法零命中 live run 和输入 hash 被改写后拒绝通过；之后集群已完成并归档原 T0/T1。
最新复核发现相邻位置竞争例的误报，修复后 45 个本地测试通过；五株分类重算尚待集群执行。

下一执行任务：按 [相邻竞争修复复核](../experiments/callout_scale_v1/NEARBY_COMPETITION_REVIEW.md)
重算固定五株审核分类、检查器、报告和差异表，再更新 release_check.json；原生发现无需重跑。
现有发现与边界算法保持不变；修复必须解释具体失败例，而非换路线或为增加产量改阈值。
