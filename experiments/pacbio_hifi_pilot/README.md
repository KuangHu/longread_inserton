# PacBio HiFi 单样本兼容性试跑（Sam8H1）

日期：2026-10-24。作业：[pilot.sbatch](pilot.sbatch)（SLURM job 27260852，n0162.lr6）。结果：[results/Sam8H1/](results/Sam8H1/)。

**范围**：这是兼容性和组装一致性检查。该 assembly 是用同一批 reads 组装的（Canu 2.0, PacBio），
所以 read 支持不是独立测序验证；一个样本跑通也不代表 HiFi 平台已全面验收。

## 样本与输入校验

| 项 | 值 |
|---|---|
| 菌株 | *Bacillus velezensis* Sam8H1（saltpan，Botswana；论文 DOI 10.1128/MRA.00175-21） |
| run | SRR13837247，BioSample SAMN17825431，Sequel II，WGS / GENOMIC / RANDOM |
| 原始提交文件 | `200909_Cell1_CCS-bc1016_BAK8B_OA--bc1016_BAK8B_OA-m54271_200909_114912.Q20.fastq.gz`（SRA XML；1 个 SMRT Cell 的 CCS Q20） |
| ENA FASTQ | `SRR13837247_subreads.fastq.gz`（ENA 文件名写 subreads，实为 CCS） |
| 字节 / MD5 | 269,402,451 / `fd288a5f52bf55203482cffca98ac1f1`，均一致 |
| SHA256 | `57d4a020eab74a9faf9cd1923c4031058c4093d311cdfccfc35e5213c049151e` |
| reads / bases | 35,763 / 277,511,453，与 ENA 一致（约 67×） |
| 读长 | 最短 124，中位 7,357，均值 7,760，N50 7,714，最长 20,786；95% ≥5,883 bp |
| 质量 | 每条 read 的平均 Q（由平均错误率换算）：最低 20.2，中位 28.3；≥Q20 100%，≥Q30 38.8%。碱基 QV 均值 72.7，中位 80，38% 为 QV93（上限）；无占位质量 |
| assembly | GCA_016904925.1 / GCF_016904925.1，Complete Genome，两者序列逐条一致，使用 GCF 序列和 RefSeq 注释（PGAP 6.10，4,026 CDS） |
| 记录 | 染色体 NZ_CP069391.1（CP069391.1）3,983,519 bp；质粒 NZ_CP069392.1（CP069392.1，unnamed1）148,194 bp；合计 4,131,713 bp，全部保留 |
| FASTA SHA256 | `725026a3d9f42a1a46403e7da8673dbeb7c36505dbcffb8a02575c68bb126869` |

[verify_assembly.py](verify_assembly.py) 按 NCBI sequence report 和 assembly data report 检查 BioSample、总长度、每条记录、染色体/质粒都保留、有注释，全部 pass（`results/Sam8H1/qc/assembly_verify.json`）。
读长和质量由 [read_stats.py](read_stats.py) 统计（`qc/read_stats.json`），不过滤 reads。

## 流程与必要修改

与 [callout_scale_v1/batch50](../callout_scale_v1/batch50/batch_sample.sbatch) 相同：assembly 发现（PILER + RepeatScout，`--subset all`）→ validate → minimap2 → cap reads → audited classifier（`nearby_competition_v2`）→ 独立 checker。
所有阈值、发现参数、分类规则不变。修改只有：

1. minimap2 预设：主比对与 recheck 都用 `-x map-hifi`（其余参数不变：`-c --secondary=yes -N 50 -p 0.8 -I 8G`；recheck `-N 200 -p 0.5`）。
2. 新增 `read_stats.py`（读长/质量）和 `verify_assembly.py`（NCBI 报告核对）两个检查阶段。
3. classify 记录 `--assembly-consistency assembly_built_from_same_reads`。
4. `scripts/long_read_runs.json` 登记 SRR13837247（role `pacbio_hifi_pilot`，platform `PACBIO_SMRT`）。

ONT 假设检查：

- 发现步骤不读 reads，与平台无关。
- `audit_filled_site.py` 和 `check_audited_evidence.py` 只读 PAF/CIGAR，没有平台常数；`scripts/` 和 `src/` 中没有 `map-ont` 或 nanopore 字样。唯一的 ONT 设定是 sbatch 里的预设。
- `fetch_long_reads.py` 已有 PacBio 读名规则（subread / ccs / zmw）。但 SRA 去掉了读名（`SRR13837247.N N/1`），所以 `molecule_grouping=molecule_grouping_unavailable`、`read_unit_counts={name_stripped: 35763}`，ZMW 无法从读名恢复。
  “每条记录一个分子”的依据只是 SRA 记录的原始文件为 CCS、每个 spot 1 条 read，而不是读名。这一点对支持单元很重要：若是 subreads，同一 ZMW 的多条 subread 会被当作多个单元。
- 降采样脚本要求 `one_record_per_molecule`，对这类 PacBio 数据会拒绝运行；本样本 67×，未降采样，不受影响。

## 结果

| 指标 | 值 |
|---|---|
| 所有阶段 | exit 0（11 个阶段） |
| validate | pass（22 个家族，136 records，0 条无效，native 输出未变） |
| checker | pass |
| 家族（tnp-supported） | 22（0）：PILER 1 个、RepeatScout 21 个 |
| records（≥400 bp） | 136（51），全部在染色体上，质粒上 0 条 |
| audited_full_span | 135（≥400 bp 的 51 条中 51 条） |
| ambiguous_competing | 1：`Sam8H1_RS_R7_m1`（NZ_CP069391.1:3,980,392–3,980,660，268 bp），39 条 read 中 38 个单元有近等分竞争比对 |
| 其他类别 | 0（无 edge_limited、one_side_only、conflict、cigar fail） |
| cap reads | 0，所以 map-hifi 的 recheck 路径**未被实际执行** |
| 比对 | 35,762/35,763 条 read 有比对；35,454 条 primary MAPQ 60；质粒覆盖深度约为染色体的 2.6 倍（167× 对 63×） |

records 是有 read 支持的候选结构，不是独立 insertion 数。

人工粗查（不属于流程）：按 RefSeq 注释重叠，候选主要是 rRNA operon（PILER 4.8 kb 家族、RS_R0/R2/R10）、tRNA 簇和 PKS/NRPS 重复结构域（RS_R1/R4/R5/R11），没有 transposase 支持。
这与 tnp-supported = 0 一致，不能当作可移动元件证据。

## 版本与资源

- 代码：`2a3261a`（`results/Sam8H1/source_commit.txt`）；multisite-insertion 0.1.0；minimap2 2.30-r1287；PALS v1.0；RepeatScout 1.0.7；blastn 2.17.0+；NCBI datasets 18.24.0；Python 3.11.17（`qc/tool_versions.txt`）。
- 最后一次作业 5 分 46 秒（16 CPU，CPU 时间 6:56，峰值 RSS 628 MB）。其中 callout 4:05（287 MB），map 0:12（474 MB），read_stats 1:00；reads 在第 1 次尝试中下载，fetch 用时 20:59。
- 集群占用：reads 269 MB，assembly 37 MB，输出 53 MB。PAF、read_evidence、member_boundaries 留在 scratch，SHA256 见 [omitted_large_files.tsv](results/Sam8H1/omitted_large_files.tsv)。

## 失败的尝试（[results/failed_attempts/](results/failed_attempts/)）

1. job 27259947（n0038）：reads 下载和统计完成后，在 `prepare_assembly` 开始时卡住 47 分钟，进程处于 D 状态，连列目录都挂起，是节点 scratch 故障。已取消；n0038 加入故障节点名单。
2. job 27260806（n0162）：`verify_assembly` 报 `KeyError: 'refseq_accession'`，datasets 的 JSON 键是 camelCase。脚本修正（2a3261a）后重跑。
3. job 27260852：通过。失败目录在集群 `pacbio_hifi_pilot_failed/`。

## 兼容性结论

现有流程在这份 HiFi CCS 数据上只改 minimap2 预设即可完整跑通，validate 和 checker 均通过，没有放宽任何阈值。
尚未验证的部分：map-hifi recheck 路径（本样本无 cap reads）；读名被去掉时的分子分组（只能靠 SRA 元数据）；IS 丰富或 draft assembly 的 HiFi 样本；与组装 reads 独立的 HiFi 数据。
