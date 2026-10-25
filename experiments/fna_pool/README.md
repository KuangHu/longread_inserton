# FNA 数据池：S0 清单与关联、S1 批量下载

计划见 [PLAN_CROSS_ASSEMBLY_GFF3_SHORT_INSERTIONS.md](../../docs/PLAN_CROSS_ASSEMBLY_GFF3_SHORT_INSERTIONS.md) 的 S0–S1。
本轮决定（2026-10-24）：

- FNA 来源：NCBI Assembly，优先 RefSeq，GenBank 补充；
- 首个 ≥100 样本数据池：Comamonadaceae（taxid 80864）；
- 内容：所有与 ONT/PacBio 通过 exact BioSample 关联的 assembly，加上科内全部 RefSeq Complete Genome（assembly-only 背景）；
- 批量下载在 lrc-xfer 数据传输节点上运行。

脚本只用 Python 标准库（3.6 可运行）和 NCBI datasets CLI。原始响应和下载文件留在
`/global/scratch/users/kh36969/longread_callout/fna_pool/<pool>/`，git 只存表格、校验值和小结。

## S0：[build_pool.py](build_pool.py)

```sh
W=/global/scratch/users/kh36969/longread_callout/fna_pool/comamonadaceae_20261024
/usr/bin/python3 build_pool.py fetch --taxon 80864 --work $W
/usr/bin/python3 build_pool.py build --taxon 80864 --work $W --out comamonadaceae_20261024
```

- `fetch`：NCBI Datasets genome summary（RefSeq、GenBank 当前版本）；ENA read_run（ONT + PacBio）；
  PacBio run 的 SRA XML（原始提交文件名，作为 read type 证据）。每个原始文件的命令/URL、行数和 SHA256 在 `raw/provenance.json`，并抄入 `summary.json`。
- `build`：
  - 可用 run：GENOMIC、WGS、非 PCR、ENA FASTQ 有 MD5 和 bytes、平均读长 ≥500 bp（排除标成长读段的短读段数据，如 PacBio Onso、150–302 bp 的 R1/R2 文件）。
  - read type：ONT 按平台；PacBio 只按原始提交文件名判断 HiFi/CCS 或 CLR/subreads，ENA 的 `_subreads.fastq.gz` 是 SRA 命名，不作证据；判断不了标 `pacbio_unresolved`。
  - 关联：run 的 sample_accession 与 assembly BioSample 精确相等 → `biosample_matched`。是否是组装所用 reads（`assembly_source_run`）本轮**未核查**。
  - 每个 BioSample 一个 assembly：当前版本 → 组装级别 → 最新发布 → accession；有 RefSeq 配对且 NCBI 未标“不同”时下载 GCF，同时下载 GCA 的 FNA 做序列一致性核对。
  - 不做来源（临床/宿主）筛、物种排除或克隆合并；这些在分 batch 时处理。

### Comamonadaceae 结果（[summary.json](comamonadaceae_20261024/summary.json)）

| 项 | 数 |
|---|---:|
| NCBI 当前 assembly | RefSeq 2,400；GenBank 6,819 |
| ENA run | 422（ONT 267，PacBio 155） |
| 可用 run | 373：ONT 243，HiFi/CCS 21，CLR/subreads 17，PacBio 类型未定 92 |
| 可用且同 BioSample 有 assembly 的 run | 258 |
| 数据池 BioSample（= assembly） | 499：linked 236，RefSeq Complete 背景 263 |
| linked 样本按 read type | ONT 163，HiFi 21，CLR 12，PacBio 未定 40（一个样本可有多种） |
| linked assembly 级别 | Complete 155，Chromosome 5，Contig 76 |
| 序列来源 | RefSeq 485，GenBank 14 |
| BioProject | 295 |
| NCBI atypical 标记 | 5 |
| 下载清单 | 984 行：主包 499（FNA+GFF3+FAA+seq report），GCA 核对 485（FNA+seq report） |

表格：`pool.tsv`（每样本一行）、`download_manifest.tsv`、`assembly_run_links.tsv`、`runs.tsv`（含不可用原因）、`assemblies.tsv.gz`（全科目录）。
linked 的 Contig 级 assembly 有 76 个，跨 FNA 调用时边界可靠性需单独评估。

## S1：[download_fna.py](download_fna.py)

每块（默认 50 个 accession）：dehydrated 包 → `datasets rehydrate` → 逐 accession 核查 → 通过后才移入 `data/genomes/<accession>/`。
每次尝试追加到 `data/state.jsonl`；重跑只处理未通过的，失败最多重试 3 次；同一时间只允许一个实例（`data/.lock`）。

核查：包内 MD5 全部一致；FNA 记录与 NCBI sequence report（accession、长度）一致且全部保留；
assembly report 的 accession、BioSample 与清单一致；GFF3 `##sequence-region` 与 FNA 一致；
GFF3 每个 CDS 的 protein_id 都在 FAA，FAA 每个 id 都在 GFF3。没有注释记为 `annotation_missing`，不算失败。
`report` 再逐对比较 GCF 与 GCA 的记录序列 SHA256，写 `pair_checks.tsv`；不一致的 pair 不去重。

小试验（4 个 BioSample、6 个 accession，登录节点，系统 Python 3.6.8）：6/6 通过，2 对 GCF/GCA 序列一致
（含一对版本号不同：GCF_964252695.1 / GCA_964252695.2），49 MB；重跑不重复下载；不存在的 accession 记为 failed 并保留 NCBI 错误信息。

### 在 lrc-xfer 上运行

lrc-xfer 需要 OTP 交互登录，只能由用户启动。它与集群共享 home 和 `/global/scratch`。

```sh
ssh lrc-xfer.lbl.gov
tmux new -s fna        # 断线后用 tmux attach -t fna 回来
cd ~/middle_insertion_callout/longread_inserton/experiments/fna_pool
P=/global/scratch/users/kh36969/longread_callout/fna_pool/comamonadaceae_20261024
/usr/bin/python3 download_fna.py run --manifest comamonadaceae_20261024/download_manifest.tsv \
    --data $P/data 2>&1 | tee -a $P/download.log
# Ctrl-b d 脱离；中断后重跑同一命令即可续传
```

估计：约 7 GB（主包约 10 MB/个，GCA FNA 约 4 MB/个），按小试验速度约半小时。完成后在登录节点运行：

```sh
/usr/bin/python3 download_fna.py report --manifest comamonadaceae_20261024/download_manifest.tsv \
    --data $P/data --out comamonadaceae_20261024/s1
```

### Comamonadaceae 下载结果（2026-10-24，lrc-xfer，[s1/](comamonadaceae_20261024/s1/)）

| 项 | 结果 |
|---|---|
| 清单 | 984/984 verified，0 failed（`s1_report.json`） |
| 主包 | 499，全部有注释：RefSeq PGAP 485，提交者注释 14 |
| GCF/GCA 序列核对 | 485/485 记录序列 SHA256 完全一致（`pair_checks.tsv`） |
| 记录 | 1,320 条（每个 assembly 中位 1 条，最多 168 条）；153 个 assembly 含质粒，共 268 条质粒记录 |
| 体积 | 2,966 个文件，8.6 GB（`files.tsv` 记录每个文件的 bytes、MD5、SHA256） |
| 用时 | 两次运行共约 19 分钟（20:01–20:11，20:24–20:31） |

- 第 1 次运行在第 9 块中断：GCF_049860745.1、GCF_050310035.1 的 sequence report 有一条记录缺 `length` 字段，脚本报 KeyError。
  修正（9fb95a2）后这类记录只按 accession 对应，并记录缺失数；总长度与 GFF3 sequence-region 长度仍核对。已通过的 422 个保留，其余重跑。
- sequence report 的 `assignedMoleculeLocationType` 对 draft 的未定位 contig 也写 Chromosome，所以 “Chromosome” 记录数不等于染色体数。

reads（S5）不在本步下载。linked 样本的可用 FASTQ 共约 213 GB（ONT 151，PacBio CLR 25，HiFi 12，类型未定 25），等候选和优先级明确后再分批下载。

## S2：[build_s2.py](build_s2.py)（[s2.sbatch](s2.sbatch)）

`sbatch s2.sbatch comamonadaceae_20261024`（job 27262241，n0323，16 分钟：catalog 14:16，skani 1:32，pairs 0:01）。表格在 [s2/](comamonadaceae_20261024/s2/)。

- **catalog.tsv**：每个 BioSample 一行，含 FNA 路径/SHA256、记录与质粒数、taxon、linked read type、clone group、species cluster，
  以及 batch50 的来源规则结果（`source_class`/`source_note`，取 BioSample host、isolation_source、env_medium）。来源只作标记，不从发现中移除。
- **联合 FNA**：集群 `s2/joint.fna.gz`，表头 `sample_id|assembly|contig`，每条 contig 独立，碱基不改；
  [contig_index.tsv.gz](comamonadaceae_20261024/s2/contig_index.tsv.gz) 记录每条的长度和 SHA256，且每条都与 S1 的记录 SHA256 逐一核对（[joint_index.json](comamonadaceae_20261024/s2/joint_index.json)）。
- **距离**：skani 0.4.0 `triangle -E --slow`（与 batch50 相同），边表留在集群，SHA256 见 [skani_run.json](comamonadaceae_20261024/s2/skani_run.json)。
  clone group：ANI ≥99.99 且双向 AF ≥95；species cluster：单链接，ANI ≥95 且双向 AF ≥50；属取 organism 名第一个词。
- **pair 表**（[pairs.tsv](comamonadaceae_20261024/s2/pairs.tsv)）：每个 assembly 取 1 个 clone 伙伴、最近 3 个同 species 非 clone 邻居、1 个最近的其他 species（须有 skani 输出）；
  另取 1 个 10–20 成员的 cluster 全对全作为审计，用来衡量邻居选择的漏检。这些是本池的初始设置，不是校准过的阈值。

| 项 | 结果 |
|---|---:|
| 样本 / contig / bp | 499 / 1,320 / 2,915,207,232 |
| 属 | 28（Variovorax 171、Comamonas 63、Acidovorax 59、Delftia 53 …） |
| skani 报告的 pair | 54,461 / 124,251（只报 AF ≥15） |
| 与任何样本都无 skani 输出 | 12 |
| clone group | 426（其中 37 个多于 1 个成员） |
| species cluster | 215（154 个单例；8 个 ≥10 成员：最大 51，Delftia acidovorans/tsuruhatensis） |
| 来源标记（全部 / linked） | environmental 174/66，plant 117/72，enrichment 21/19，clinical/animal 27/11，unclear 160/68 |
| S3 pair | 1,063：同 species 邻居 551，clone 伙伴 71，最近其他 species 425，审计全对全 45；两端都 linked 401 |
| linked 样本的 pair 覆盖 | 179/236 有同 species 伙伴；52 只有其他 species 伙伴；5 没有任何 pair |

- 审计 cluster 是 10 株 Variovorax sp. OK*（全部 linked），但都来自同一 BioProject PRJNA1327437，作为漏检审计代表性有限。
- 只有其他 species 伙伴的样本，跨 FNA 比较的局部共线性可能不足；没有 pair 的 12 个样本（5 个 linked）无法进入 pairwise 路线，留给单株发现或多 FNA 图方法。
- unclear 主要是 BioSample 无 host/来源字段；clinical/animal 的 27 个仍在池中，报告时分层。

## S2 复核修正（2026-10-24）

- 无 pair 的样本记为 `no_pair_selected_pending_review`：skani 无输出（约 80% ANI k-mer 预筛、只报 AF ≥15）不等于没有可用的局部同源区域。现 10 个（linked 3 个）。
- species cluster 是单链接计算分组；每个 pair 另记 `pair_meets_species_threshold`。同组但本身不满足 ANI ≥95、双向 AF ≥50 的 pair 有 60 个（如 P00570，94.92%）。S3 只按该 pair 的实际比对判断。
- 来源规则 v2：plant 关键词加词边界（原 `stem` 会命中 system），“<设施> plant”（污水处理厂、炼油厂等）不算植物；原始 host/isolation_source/env_medium 保留。7 个标签改变（plant→environmental 6，plant→unclear 1）。batch50 冻结面板的 plant 标签都不来自这些词，不受影响。
- 审计集：原单项目 Variovorax 组（一致性检查）之外，加跨项目组（19 株 Comamonas thiooxydans/testosteroni，19 个 BioProject，171 对），以及同属、不同项目的距离分层对照（ANI 99.99–100、99–99.99、97–99、95–97、90–95、<90、无 skani 输出，各 5 对）。pair 共 1,221，77 对有多个理由。

## S3 原型：[s3_pair_calls.py](s3_pair_calls.py)

- `call`：两个基因组互相比对（minimap2 2.30 `-x asm20 -c --eqx --secondary=yes -N 20 -p 0.5`），记录 ≥50 bp 的差异：同一比对内的 I/D 段，或同一 contig 对、同链的两个共线 block 之间的 gap。
  每条原始 call 保存两侧坐标、比对方向和链、两侧 200 bp flank 的 identity/覆盖、contig 末端、tandem 样、边界滑动窗口、另一侧 gap（replacement）。只写“carrier 有 X、空位点无 X”，不推断插入/删除方向。
- probe：每条 call 的左右 flank 和 core 作为短序列分别比对两个基因组（`map-pb`；asm/HiFi 预设会隐藏短序列的等分次级比对）。任一 flank 在任一基因组有 >1 个位置 → `flank_not_unique`，不锚定；记录 core 在两侧的拷贝数。
- 一侧 flank：两个基因组在比对末端之外都还有未比对序列（如紧邻 contig 末端）→ `one_flank.tsv`，unresolved，不进 locus 表。
- `merge`：同 carrier contig 上，相互重叠 ≥0.8 或两端都在滑动窗口+20 bp 内才合并；相邻但不重叠的 core 不合并。`boundary_consistency`：exact / within_shift / partner_dependent（不同伙伴边界不同）/ direction_dependent（同一伙伴、两个比对方向边界不同）；原始 call 全部保留。`anchored_support`：all_calls_anchored / partner_mixed（另一伙伴的 call 是 replacement/unresolved）/ direction_mixed（同一伙伴的另一方向不锚定）。`max_other_side_gap(_anchored)`：replacement 判定允许另一侧 ≤ max(10 bp, core 的 10%)，这一列显示锚定 call 实际带了多少另一侧序列。
- 阈值（MIN_CORE 50、FLANK 200、identity/覆盖 0.9、EDGE 1 kb 等）是原型设置，不是校准值。

挑战测试 [test_s3.py](test_s3.py)：13/13 通过——单一插入（两个方向）、反向互补、相距 300 bp 的两个插入、tandem 复制、重复 flank、contig 末端、replacement、同一 core 在不同位置、两个 pair 的同一位点合并、60 bp 调用而 30 bp 不调用、20 kb、60 kb、两 block PAF（split-block 路径；minimap2 把 60 kb 也放在一个 CIGAR 里）。

### 12 对小规模试跑（[s3_pilot/](comamonadaceae_20261024/s3_pilot/)，job 27262435，6 分钟，每对 19–34 s、<1 GB）

| pair | 选取理由 | ANI | loci | 锚定 | 其中 <2500 bp | replacement | unresolved | one-flank |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| P00595 | 距离层 99.99–100 | 100.0 | 1 | 0 | 0 | 0 | 1 | 0 |
| P00351 | 99–99.99 | 99.47 | 65 | 21 | 15 | 15 | 29 | 12 |
| P00290 | 97–99 | 97.69 | 541 | 157 | 113 | 115 | 269 | 9 |
| P00034 | 95–97 | 95.49 | 761 | 140 | 111 | 151 | 470 | 10 |
| P00111 | 90–95 | 90.64 | 1,057 | 129 | 105 | 131 | 797 | 107 |
| P00384 | <90 | 88.68 | 1,905 | 53 | 40 | 201 | 1,651 | 127 |
| P00063 | 同属无 skani 输出 | – | 34 | 0 | 0 | 1 | 33 | 102 |
| P00136 | clone 伙伴 | 100.0 | 3 | 3 | 1 | 0 | 0 | 0 |
| P00033 | 跨项目同 species 邻居 | 98.88 | 250 | 74 | 57 | 70 | 106 | 10 |
| P00022 | 最近其他 species | 93.28 | 569 | 88 | 54 | 121 | 360 | 17 |
| P00007 | 跨项目审计组 | 95.53 | 437 | 74 | 48 | 113 | 250 | 37 |
| P01143 | 单项目审计组 | 100.0 | 0 | 0 | 0 | 0 | 0 | 0 |

合计 10,198 条 call → 5,521 个 locus：锚定 704（<2500 bp 525），replacement 881，unresolved 3,936，另有 one-flank 431。
（loci 按 carrier 分；同一 carrier 位点被多对看到会合并，所以各行之和大于合计。）

- 距离越远，差异越多但能锚定的越少：<90% ANI 的 pair 有 1,905 个 locus，只有 53 个锚定；unresolved 主要是 flank identity/覆盖不足（3,842 个，flank identity 中位 0.81）。
- 同一 pair 的两个比对方向给出相同边界；锚定 locus 中 644 个两个方向都调用，25 个只有一个方向。
- 锚定 locus 的边界：exact 453，within_shift 211，direction_dependent 33，partner_dependent 7（**更正**：之前把 40 个都记作 partner_dependent、且说 <2500 bp 有 18 个；拆分后 <2500 bp 的 15 个全是同一伙伴两个方向的差异，差 20–40 bp，来自 P00290/P00111 这类较远配对的分歧连接区；7 个真正 partner_dependent 都是 10–77 kb 的岛）。
- flank 不唯一而不锚定的有 233 个，tandem 样锚定 9 个。
- 预览（不是 S4）：525 个锚定 <2500 bp core 中，68 个恰含 1–2 个完整 CDS 且边界不切断任何注释特征，产物多为 IS3/IS481/ISL3/IS5 transposase、type II toxin-antitoxin、hypothetical protein；其余多数是 0 个完整 CDS 或边界切断基因。
- 图：`plots/`（每对按状态各取几个）和 `plots_short_cds/`（上面 8 个 1–2 CDS 示例），SVG，carrier 窗口（GFF3 CDS/rRNA/tRNA，U/core/D）＋伙伴空位点窗口。

### 复核后补充（2026-10-25）

- **锚定不等于每个 call 都锚定**：`status` 取最好的 call。704 个锚定 locus 中 all_calls_anchored 634，direction_mixed 56，partner_mixed 14；<2500 bp 的 525 个中分别为 479 / 42 / 4。
- **10% replacement 容许**（v0.2 起另有 difference_class，见下）：锚定 call 的另一侧 gap，<2500 bp 中 0–10 bp 352，11–50 bp 144，51–250 bp 29；全部锚定中 >250 bp 的 36 个都是长 core（如 91 kb 岛带 4,673 bp 另一侧序列仍算锚定）。阈值未改，只把数值写进 loci.tsv。
- **实际渲染检查**（rsvg-convert 转 PNG 逐张看）发现并修正：轨道标题压住基因名；负链伙伴未翻转、与 carrier 基因顺序相反；短 core 的 U/core/D 标签重叠；横跨空位点的伙伴基因被拉伸成整段 core 长（如 44 kb 噬菌体位点的 tRNA-Leu）；GFF3 产物名未 URL 解码。现在伙伴轨道按 carrier 方向画，U'/D' 对齐 U/D，每个伙伴/边界组合一条轨道。
- 新图：`plots_boundary/`（7 个 partner_dependent ＋ 15 个 <2500 bp direction_dependent），`plots_readcheck/`（下面 7 个 mixed/contradicted）。

### 长读段等位检查 v0.1（job 27263902；**已被下面 v0.2 取代**，结果保留在 `readcheck/` 供对照）

v0.1 只看比对起止坐标、不看 CIGAR，也不让重复位点在全基因组竞争，所以下面的 spanning 数不能当作完整结构支持数（见 v0.2）。

每个锚定 locus 的每条原始 call 一个检验：在**各样本自己的组装**上构建 present（U+core+D）和 absent（U+另一侧+D）两个局部等位（各 2 kb flank），用该 BioSample 关联的 run 比对（ONT map-ont，HiFi map-hifi，其余 PacBio map-pb），只数**跨过两个连接点各 300 bp（加滑动窗口）**的读段；两个等位都跨过时按 AS 差 ≥20 判定。判定规则预先固定（≥3 条跨越读段；另一等位 ≤10% 为 consistent），没有按试跑调。

**这些读段通常就是组装用的读段，所以这是组装一致性检查，不是独立测序验证。** 12 个 run 下载后 MD5/字节/碱基数全部核对通过（verify_*.json）。4 个 pilot 样本没有关联 run（P00290、P00111、P00384、P00007 只有一侧有读段）。构建检查：1,938 个伙伴等位中 1,925 个与同一检验的 carrier 对应等位最匹配，13 个（12 个检验）无命中或不一致。

| 检验（969 个，carrier/伙伴） | 数 |
|---|---:|
| 两侧都 consistent | 294 |
| 只有一侧有读段且 consistent | 617 |
| carrier 跨越读段不足（长 core），伙伴无读段或 consistent | 51 |
| mixed | 6 |
| contradicted | 1 |

- <2500 bp 的 525 个锚定 locus：两侧 consistent 220，一侧 consistent 301，mixed 3，contradicted 1；有读段的一侧没有一个是读段不足（跨越读段中位数 carrier 110、伙伴 111）。15 个 direction_dependent 都只有一侧有读段，读段支持 carrier 的 core，但区分不了两个方向的边界。
- contradicted：L_dc88a873e166，216 bp tandem 样 core，在 beta-strand repeat 蛋白基因内；伙伴 SAMD00407106 的 240 条跨越读段全部支持“有 core”，伙伴组装很可能塌缩了一个重复单元——这是组装假象，不是真实差异。
- mixed：L_02729f48424a（318 bp，carrier 有 6 个拷贝，carrier 读段 51 present / 29 absent）和 L_43ea79084943（1,301 bp，33/17）可能是样本内亚群或重复相关；L_128138f94820（54 bp）carrier 209 present / 193 absent / 54 ambiguous（**更正**：之前只报了 54 条 ambiguous）；另外 3 个是 4–39 kb 长 core，跨越 present 等位需要更长读段，计数天然偏向 absent，不能直接解读。
- PacBio 计数单位：SRA 规范化的 subread 同一 spot 同名，按一条计；DRR396223 名字里没有 ZMW 信息，可能把同一分子的多个 subread 重复计数。

### 长读段等位检查 v0.2（复核 0ed1ca4 后；job 27272478，commit 830efde，10 min，读段已在本地）

复核指出并确认的问题：(1) 只看比对起止，CIGAR 中整个 core 被删除仍算 present；(2) 读段只和候选局部等位比，重复位点之间不竞争，一条读段可同时支持强、弱两个位点；(3) “一侧 consistent”混了两种不同证据；(4) 长等位需要更长读段，计数偏向短等位；无分子身份的数据只能算 read record。另外本轮发现：ONT duplex run（ERR15680828、ERR15686901）里 duplex 读段（`id1;id2`）和它的两条 simplex 亲本同时存在（约 6,000 组），同一分子最多被计 3 次。

v0.2 改动（判定阈值 ≥3、≤10% 不变）：
- **CIGAR**：检验区间内合并后的 indel 事件、以及扣除该读段区间外自身 indel 偏差后的净 indel，都必须小于两个等位长度差的一半（至少 10 bp）。在 v0.1 的比对上预检：支持组装所示等位的跨越读段因此丢掉 0–4%（CLR 最高），而支持另一等位的跨越比对有 22–100% 不通过。
- **全基因组竞争**：命中等位的读段再比对到该样本整个组装（保留 secondary）；只有基因组位置覆盖该位点 ±300 bp、且别处没有覆盖同一读段区段（≥50%）、得分差 <20 的比对时才计（否则 off_locus / locus_ambiguous）。一条读段只有基因组位置同时覆盖几个位点时才同时支持它们（相邻真实位点）。
- **分子计数**：PacBio 按 ZMW 或 SRA spot，ONT 按 read id 并把 duplex 与亲本合并；DRR396223 无分子信息，标为 `read_records`。
- **长度匹配**：判定只用读长 ≥ 两个检验区间中较长者的分子（原始数也列出）；读长接近该长度的读段仍更容易跨过短等位。
- **类别名**：检验分 both_samples_consistent / carrier_only_consistent / partner_only_consistent / mixed / contradicted / not_tested；locus 证据等级 E2（所有检验两侧）/ E2s（部分检验两侧）/ E1c（只有 carrier 读段）/ E1p（只有伙伴读段；carrier 组装未被检查）/ E0 / X（mixed 或 contradicted）。
- **difference_class**（`merge --fasta`，status 不变）：边界规范化后（另一侧序列若正好等于 core 两端则吸收）分 presence_absence、complex_small_other（另一侧 <50 bp，即 MIN_CORE）、replacement（≥50 bp，本身就够一个事件）。规范化实际一个 call 也没吸收：锚定 call 1,369 条中 0 bp 453，1–2 bp 104，3–10 bp 296，11–49 bp 324，≥50 bp 192。locus 取锚定 call 中最复杂的一类：704 个锚定 locus 中 presence_absence 214，complex_small_other 386，replacement 104（<2500 bp：168 / 327 / 30）。

| locus（704 / <2500 bp 525） | 全部 | <2500 bp |
|---|---:|---:|
| E2 所有检验两侧 consistent | 267 | 218 |
| E2s 部分检验两侧 consistent | 7 | 3 |
| E1c 只有 carrier 读段 consistent | 171 | 139 |
| E1p 只有伙伴读段 consistent | 222 | 163 |
| E0 未检验（无读段或长度匹配后不足） | 33 | 0 |
| X mixed / contradicted | 4 | 2 |

<2500 bp 按 difference_class：presence_absence E2 87 / E1c 30 / E1p 46 / X 1；complex_small_other E2 118 / E2s 3 / E1c 101 / E1p 109；replacement E2 13 / E1c 8 / E1p 8 / X 1。carrier 计数单位为 read_records（DRR396223）的 <2500 bp locus 有 35 个（E2 17，E2s 3，E1c 15）。

和 v0.1 比：
- L_128138f94820（54 bp）：v0.1 carrier 209/193/54 是 mixed；v0.2 中 220 条 carrier 读段 off_locus——它们只在连接点附近约 600 bp 有比对（AS 500–700），基因组位置在下游 4–10 kb 的另一个重复拷贝，其中 199 条原先被算作 absent。剩下 204 present / 9 absent，两侧 consistent（E2）。v0.1 的 mixed 是近旁重复造成的错分配。
- L_02729f48424a（318 bp，6 拷贝）：29 条 absent 全部 off_locus（其他拷贝），剩 51 present → E1c。
- L_dc88a873e166（216 bp tandem 样）：CIGAR 和全基因组竞争后仍是 carrier 54 present、伙伴 240 present（0 absent）→ contradicted；**疑似**伙伴组装塌缩一个重复单元，仍未独立证实。
- 仍 mixed：L_43ea79084943（1,301 bp，replacement 类，carrier 匹配后 33/12），L_e4aff037c528（4.2 kb，17/21），L_647ba3ad9f35（11.6 kb，16/16）。L_d704e2a166a1（38.7 kb）长度匹配后 carrier 不足 → E1p。
- 22 个 24–97 kb core 从“一侧 consistent”变为 E0：伙伴的 absent 读段都短于 present 等位的检验区间，长度匹配后不再计数——这正是长度偏倚。
- 有等位命中的分子中 16,292 个同时支持 ≥2 个位点；按规则它们的基因组位置覆盖每个位点（相邻位点，或不同 carrier 的 locus 落在同一个伙伴空位点）。排除的分子-检验对：cigar_fail 1,832，off_locus 538，locus_ambiguous 195，分子内记录不一致 3。

仍是组装一致性检查，不是独立验证。E1p 只说明伙伴空位点有读段支持，不能排除 carrier 组装假象。

### 长读段等位检查 v0.3（复核 1e6243e 后；commit b9b2e7b，复用 job 27272478 的比对重新汇总，结果在 `readcheck_v03/`）

复核用两个构造记录复现了 v0.2 的定位漏洞：v0.2 只要求某条基因组比对**碰到**位点窗口，(a) 读段的另一段碰到窗口 1 bp、实际支持等位的那段在别处，仍算 ok；(b) 更强的重复拷贝比对只碰到窗口边缘，就被豁免竞争。v0.3：局部比对保留链和两个锚点（检验区间两端，U 侧和 D 侧，都在 core 之外）对应的读段位置；每个锚点必须由同链的基因组比对放在组装中预期位置 ±50 bp（LOCUS_TOL）内，且没有别的比对把同一读段位置放到别处、得分差 <20。锚点在 core 外，所以另一等位的读段不会因为不匹配组装里的 core 被排除。测试加入这两个反例（并确认 v0.2 规则接受它们）和一条负链读段，24/24。

| locus（704 / <2500 bp 525） | v0.2 | v0.3 | v0.3 <2500 bp |
|---|---:|---:|---:|
| E2 | 267 | 269 | 218 |
| E2s | 7 | 7 | 3 |
| E1c | 171 | 172 | 140 |
| E1p | 222 | 222 | 163 |
| E0 | 33 | 33 | 0 |
| X | 4 | 1 | 1 |

只有 3 个 locus 改变等级，都是 v0.2 的 mixed：L_647ba3ad9f35（11.6 kb）、L_e4aff037c528（4.2 kb）→ E2，L_43ea79084943（1,301 bp）→ E1c。它们的 absent 读段（47 / 35 / 17）现在都是 off_locus：一个锚点在预期位置，另一个锚点在基因组上偏开约一个 core 长度（−11,616 / +4,224 / +1,271 bp）。即这些读段是从一侧伸进 core 的 present 等位部分读段，core 末端与对侧 flank 相似（像移动元件的末端重复），局部 absent 等位因此能被冒充。剩下唯一的 X 是 L_dc88a873e166（216 bp tandem 样，伙伴 240 条读段全部有 core），仍是**疑似**伙伴组装塌缩。排除的分子-检验对：cigar_fail 1,832，locus_ambiguous 1,456（v0.2 为 195），off_locus 691。

计数比例的限制不变：长度匹配只去掉“根本跨不过去”的读段，接近门槛时短等位仍更容易被完整跨越，所以 mixed 的生物学解释和“另一等位 ≤10%”都带着这个偏倚。

**E2 不等于目标 insertion。** `s3_pilot.py crosstab` 写出 `short_loci_review.tsv`（<2500 bp 锚定 locus）。CDS 按 GFF3 的 ID 归并为编码单元（多段 CDS，如程序性移码的 IS3 transposase，算一个单元），每个单元带 partial（`partial=true`/`start_range`/`end_range`）和 pseudo 标记、protein id。cds_class：边界切断 CDS/rRNA/tRNA → boundary_cuts_feature；否则 core 内有 pseudo → contains_pseudo，有 partial → contains_partial，再按单元数 0 / 1–2 / >2。“complete”只表示注释既非 partial 也非 pseudo 且整体落在 core 内，**没有检查起止密码子**。

525 个短 locus：boundary_cuts_feature 287，0_CDS 164，1–2_complete_CDS 67，contains_pseudo 5，>2 2（**更正**：上一版按坐标数完整 CDS、不看 pseudo/partial、不归并分段，得到 72；更早的预览 68 是未入库的临时计算）。218 个短 E2：

| 短 E2 | 1–2 complete CDS | contains_pseudo | 0 CDS | 切断特征 | >2 |
|---|---:|---:|---:|---:|---:|
| presence_absence（87） | 13 | 1 | 36 | 37 | 0 |
| complex_small_other（118） | 21 | 0 | 33 | 63 | 1 |
| replacement（13） | 3 | 0 | 0 | 10 | 0 |

严格复核清单（E2、presence_absence、1–2_complete_CDS、边界 exact、非 tandem、carrier 按分子计数）：**13 个 locus，15 个编码单元，12 个不同 protein id**；其中 8 个 locus 含 transposase（IS30 1、IS481 1、ISNCY 1、ISL3 2、IS3 3——3 个 IS3 都是同一 WP_503779473.1 的两段移码单元，2 个 ISL3 是同一 protein），其余为 Fic family 2、HD domain 1、AbiEii/AbiGii 毒素-抗毒素 1、hypothetical 1。另列 1 个 contains_pseudo（L_13321b5da82e，IS3 transposase，frameshifted）。**更正**：上一版说 14 个、其中 10 个 transposase；按 locus 计 transposase 应为 9（含 pseudo 那个），现排除 pseudo 后为 8。locus 数、编码单元数和以后的蛋白家族数分开统计。这是复核清单，不是 S4 结论。17 个短 E2 的 carrier（DRR396223）只能按 read record 计数。

## S3 全量候选发现（1,221 对）

`s3_pilot.py allpairs` 把 s2/pairs.tsv 全部 1,221 对（489 个样本）写成同样格式；`s3_full_calls.sbatch`（25 个 array 任务，每对调用后删除可重建的整基因组 PAF）→ `s3_full_merge.sbatch`（合并、difference_class、逐对汇总、短 locus 复核表；大表留在 scratch，repo 只放汇总、短锚定 locus 表和 sha256）。距离分层对照照常调用，汇总按 reasons 分开报告。发现阶段不要求有读段。为全量合并，merge 只和仍可能匹配的 cluster 比较（结果与 pilot 合并逐列相同）；summary 改为按 pair 索引（pilot 输出逐字相同）。

尚未做：全量 read check；S4 的正式注释/边界复核；core 蛋白合并。
