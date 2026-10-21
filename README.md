# insertion-site-recovery

真实 CCS pilot 已运行：结果与配对限制见 [实测报告](docs/REAL_HIFI_PILOT.md)。后续 [长插入审计](docs/HIFI_PILOT_AUDIT.md) 显示，现有流程尚不能直接用于完整 IS 挖掘。

公共长读段下载方法与 8 个真实文件的试下载记录见 [下载规范](docs/DOWNLOADING_LONG_READS.md)，配套脚本位于 `tools/downloads/`。

独立的最小原型：输入 **MAG／assembly + HiFi／CCS 或 Nanopore reads**，恢复插入序列与空位点表示。
不需要一个物种的大规模 genome group；需要同一位点有可定位的长读段证据。

本项目参考 `KuangHu/nanopore_insertion_callout` 的“比对 → SV 候选 → 提取”流程，重新实现序列与坐标恢复。
原仓库未修改，参考版本为 `b5101b490e833a220626a15e677c10ccae4b8d9c`。
不包含 tail–head/circle 筛选、蛋白注释、RNA 筛选、自动下载或按属名选择参考。

## 输出的含义

对于每个事件，输出：

- `insertion.fasta`：插入序列；INS 来自同位点 reads 的 POA consensus；DEL 来自 assembly 中被跳过的区间。
- `empty_site.fasta`：`[upstream][downstream]`，两段均来自指定 assembly。
- `upstream.fasta`、`downstream.fasta`：可单独使用的两侧序列。
- `events.tsv`：参考坐标、插入位置、序列长度、两种状态的 read 数和结果状态。
- `records.jsonl`：完整序列、read 名、序列来源、替代断点、重复引起的边界范围。
- `manifest.json`：输入 SHA256、实际命令、软件版本、运行状态。失败写入 `failed`，不会伪装成功。

**empty-site 是参考序列上的重建表示，并有单独的局部 read 支持字段；不是菌株分型后的 empty-allele consensus。**
默认只验证断点两侧各 100 bp 的局部锚定；不能由此宣称整段 2 kb flanks 已被 reads 校正。
`read_supported_pair` 表示局部 occupied 和 empty 两种结构都有 read 支持，不保证属于同一菌株或同一细胞。
INS/DEL 是相对参考的描述，不能判定进化上是插入还是删除，也不自动证明该序列是可移动元件。

## 两条路线

| 输入 reads | 比对 | 候选 | 序列恢复 |
|---|---|---|---|
| PacBio HiFi／CCS | minimap2 `map-hifi` | Sniffles2 INS + DEL | 同一个位点锚定恢复模块 |
| Nanopore | minimap2 `map-ont` | Sniffles2 INS + DEL | 同一个恢复模块，使用 ONT identity 门槛 |

参考为空位点时，INS 的 read 序列给出 insertion；参考已占据时，DEL 的 reads 给出空连接，reference interval 给出 insertion。
后者不能因为 VCF 写着 DEL 就丢掉。

## 安装

建议 Linux、Python 3.10+。以下以 Python 3.12 为例，minimap2 需要在 PATH 中：

```bash
conda create -n insertion-recovery -c conda-forge -c bioconda python=3.12 minimap2=2.31 pip
conda activate insertion-recovery
cd insertion_site_recovery
python -m pip install '.[caller,test]'
```

也可使用已有 minimap2，或通过 `--minimap2 /path/to/minimap2` 指定。
`--sniffles /path/to/sniffles` 同理。`requirements-tested.txt` 记录本次验证的 Python 依赖版本。
Sniffles 固定为 2.8.1；其他版本的行为需要重新验证。

## 用法

HiFi／CCS：

```bash
insertion-site-recovery run \
  --reference assembly.fasta --reads sample.hifi.fastq.gz \
  --platform hifi --sample sample_hifi --threads 8 --out out_hifi
```

Nanopore：

```bash
insertion-site-recovery run \
  --reference assembly.fasta --reads sample.ont.fastq.gz \
  --platform ont --sample sample_ont --threads 8 --out out_ont
```

`--platform ccs` 是 `hifi` 的别名。`--reads` 接受多个 FASTQ/FASTA 文件，支持压缩 reads。
参考必须是未压缩 FASTA。HiFi／CCS 和 ONT 应分别运行，不混为同一个平台处理。

已有同一 assembly 的坐标排序、索引 BAM 及 VCF 时：

```bash
insertion-site-recovery recover \
  --reference assembly.fasta --bam sample.bam --vcf calls.vcf \
  --platform ont --sample sample_ont --out out_recovered
```

`recover` 检查 BAM 的参考名称、长度和已有 M5。缺少 M5 时无法证明序列完全一致，必须由输入提供者保证 BAM 与 FASTA 相符。
输出目录必须不存在或为空，不自动复用旧文件。输入 reference 复制到输出目录后再索引，不修改原输入。

`run` 保留 BAM、VCF 和运行日志。比对暂存 SAM，排序成功后删除，因此需要为 SAM/BAM/reference 预留空间。
建议使用样本匹配的 assembly。完整 metagenome assembly 可保留更多竞争比对位置；单个 MAG 无法排除未提供的近缘基因组。
目前没有 MAG bin/taxonomy 自动归属模块：以输出 contig ID 对接自己的 bin 表。

## 坐标约定

全部采用 **0-based、half-open 区间、interbase 插入位置**，方向为 assembly forward。
VCF POS 只用于定位候选附近，最终断点由支持 reads 的 CIGAR 确定，不直接照抄 VCF。

| 参考状态 | `reference_start0` | `reference_end0` | 意义 |
|---|---:|---:|---|
| empty，INS call | p | p | 插入在参考第 p 个与第 p+1 个碱基之间（按 1-based 描述） |
| occupied，DEL call | s | e | insertion 对应参考区间 `[s,e)` |

`empty_site_insertion_pos0 = len(upstream)`，是 **empty-site FASTA 自己的坐标**，不是原参考坐标。
例如上下游各 1,000 bp，则 empty-site 长 2,000 bp，插入点为 1,000：

```python
occupied_representation = empty_site[:position0] + insertion + empty_site[position0:]
```

重复序列造成等价边界时，JSON 的 `breakpoint_equivalence` 报告替代 start 的范围。
当前输出是其中一个参考方向表示；不能把它当成唯一的生物学切点。

## 结果分层

| status | 含义 | 序列目录 |
|---|---|---|
| `read_supported_pair` | INS/DEL 事件至少 3 条独立 read，另一状态至少 2 条；flanks 完整，边界无检测到的等价重复 | `accepted/` |
| `reference_backed_only` | 可恢复序列，但另一状态的 read 支持不足 | `candidates/` |
| `ambiguous` | 可以构建序列表示，但存在重复导致的等价边界 | `candidates/` |
| `partial` | 两种状态有支持，但 contig 边缘截短了要求的 flanks | `candidates/` |
| `unresolved` | 证据不足、多个受支持断点、非 INS/DEL、非 PASS 等 | TSV/JSON 中保留原因，不输出猜测的序列 |

状态有优先级；完整细节还应查看 `support`、`full_requested_flanks`、`breakpoint_equivalence`。
例如 `ambiguous` 不保证另一状态已经达到支持门槛。

默认参数：flank 1,000 bp、anchor 100 bp、MAPQ ≥20、事件长度 50–100,000 bp、
事件支持 ≥3、另一状态支持 ≥2；局部 identity 为 HiFi 95%、ONT 90%。
锚定和参考状态验证允许长度 ≤10 bp 的局部 I/D，并继续检查 identity；N/skipped region 和较大 gap 不能当成桥接证据。
候选搜索 ±100 bp，长度容差 25%。仅在同一个精确断点组内生成 consensus，最多使用 100 条确定性选择的 reads。
所有该组 read 的插入序列仍需与 consensus 达到 identity 门槛；这不能完全排除相近混合等位型。

Sniffles 使用 `--all-contigs --mosaic --mosaic-include-germline`，包含短 contig 和高频事件。
`--mosaic-af-min` 默认 0.05；这是候选筛选参数，**不是**保证可检测到 5% 的宏基因组等位型。
Sniffles 的人类 SV/mosaic 模型不等于已验证的宏基因组菌株模型。

## 当前适用边界

- 这是 full-CIGAR、双侧局部锚定的第一版。仅 split-read、超长插入、复杂重排和需要局部组装的候选会 unresolved；没有伪装成成功的 fallback。
- 参考与 reads 必须在该局部足够接近；read 需跨过插入及两侧锚点。DEL 的 occupied 支持也要求跨过整个被删除区间。
- 重复、低覆盖、稀有等位型、MAG 错拼、跨物种误比对、污染和参考中缺失的近缘序列仍可影响结果。
- 完全固定的 insertion，若样本没有 empty 等位型、参考也不给出有效对照，不能凭 reads 凭空恢复“实际存在过的空位点”。
- 不进行 strain phasing 或 empty-site SNP polishing；不能当作完整双等位型组装器。
- 若同位点有多个足够支持的断点，不静默挑一个。重复造成的 1 bp 等价边界也会降为 ambiguous，结果可能偏保守。
- 分子独立性目前按 read 名去重；无法替代实验级去重。
- `accepted` 是本工具的筛选结果，不是实验验证结论。

## 验证与复现

```bash
python -m pytest -q
python examples/make_synthetic.py demo
insertion-site-recovery run --reference demo/assembly.fasta --reads demo/hifi.fastq \
  --platform hifi --sample synthetic_hifi --threads 2 --out demo/result_hifi
insertion-site-recovery run --reference demo/assembly.fasta --reads demo/ont.fastq \
  --platform ont --sample synthetic_ont --threads 2 --out demo/result_ont
python examples/check_synthetic.py demo/truth.json demo/result_hifi/records.jsonl
python examples/check_synthetic.py demo/truth.json demo/result_ont/records.jsonl
```

本次验证：Python 3.12、minimap2 2.31-r1302、Sniffles 2.8.1，23 项单元/回归检查。
两条完整流程各恢复 1 个 INS 和 1 个 DEL；本组合成数据的 insertion identity 均为 100%，empty-site 和坐标与真值相同。
人工序列无生物学注释，每平台 400 条约 8 kb reads。模拟 HiFi 替换率 0.1%、indel 率 0.01%；
模拟 ONT 替换率 1.5%、indel 率 0.3%，不代表真实平台的完整误差分布。
没有进行真实宏基因组 benchmark，不能据此推断真实 recall、precision 或可发现数量。

本次受限进程环境中，Sniffles 可选 psutil 资源监控报 `NoSuchProcess`；在独立测试虚拟环境移除 psutil 后，使用其原生无监控分支完成测试。
未修改 Sniffles 算法或原仓库代码。正常计算环境不需要因此移除 psutil。

官方工具：<https://github.com/lh3/minimap2>、<https://github.com/fritzsedlazeck/Sniffles>。
