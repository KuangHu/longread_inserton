# Call-out scale v1：固定五株验收结果

主规范：[CALLOUT_SCALE_PROTOCOL.md](../../docs/CALLOUT_SCALE_PROTOCOL.md)。集群交接：[CLUSTER_TASK.md](CLUSTER_TASK.md)。
本目录是单株候选导出与本株 Nanopore filled-site read 证据的扩容验收，不是新的 discovery 算法。
运行：2026-10-23，LBL Lawrencium（pc_rubinlab / lr6）。发现与边界算法、PILER/RS 默认参数均未改动。

**当前判定：[release_check.json](release_check.json) `ready_for_10_sample_batch = false`**：
下述 b772dc7 归档验收曾通过，但补充反例暴露了同 contig 相邻竞争比对被忽略的问题。
代码已修复为 `classifier_revision=nearby_competition_v2`，本地 45 个测试通过；
需用五株原始 PAF 重算审核分类、checker、报告和差异表后恢复扩容判定。
详情及执行步骤见 [NEARBY_COMPETITION_REVIEW.md](NEARBY_COMPETITION_REVIEW.md)。
下文五株数字保留为修复前结果，不代表修复后产量。
非参考事件、空位点和候选产量不是验收门槛，本轮也未据此判定。

## T0：软件与计算对照

| 项 | 结果 | 证据 |
|---|---|---|
| 计算节点 doctor | pass（8 个原生工具） | [results/t0/doctor.json](results/t0/doctor.json) |
| 包测试 | 计算节点 42 通过（修正前代码）；最终代码 43 通过 | [unit_tests.log](results/t0/unit_tests.log) |
| 原生 demo（真实 PALS/PILER/RS/BLAST） | run＋validate pass | [demo/validate.json](results/t0/demo/validate.json) |
| 审核版挑战测试 | 17 个真值例全部通过，并经独立 checker 复算 | [tests/test_audit_filled_site.py](../../tests/test_audit_filled_site.py) |

挑战例覆盖：完整跨越、单 read、反链、X 缺失（单独及与完整共存）、secondary 跳过 X、小 indel 累积、
anchor 缺失、CIGAR N、等分竞争、candidate 处为 secondary、secondary 上限（有/无重比对）、
重复派生 read、两条单侧 read、contig 边缘、仅内部、合法零命中与畸形输入、篡改分类被 checker 拒绝。

## T1：固定五株

| 样本 | 导出 records | 审核版 audited_full_span（含冲突） | 其它审核类 | 与 v1 判定不同 | 资源（filled-site 总计） |
|---|---:|---:|---|---:|---|
| PaAR2038 | 762 | 762 | — | 0 | 9.9 min，峰值 7.2 GB（map） |
| SMT4a | 333 | 317（+11） | cigar_fail 2，edge 3 | 1 | 5.7 min，4.7 GB |
| LjRoot223 | 287 | 287 | — | 0 | 0.5 min，1.4 GB |
| R-93128 | 491 | 459 | ambiguous_competing 21，edge 11 | 4 | 1.1 min，2.2 GB |
| R-91497 | 491 | 491 | — | 0 | 1.9 min，2.2 GB |

- **AR20-38 冻结输入回归**：同一 manifest 以 `--methods PILER RS --subset all` 重跑，149 家族、762 成员的
  native/inferred/candidate 坐标、边界来源、侧翼及 X/up/down/triplet 序列逐条完全一致
  （[regression.json](results/t0/ar2038/regression.json)）。RS_R1 11 个 IS5 成员每个 135–167 个审核支持单位。
- 五株 `scripts/check_audited_evidence.py` 全部 pass（由逐 read 表独立复算单位与分类）；
  `scripts/check_callout_bundle.py` live 模式全部 `live_exports_and_v1_tables_consistent`。
- 五株均无 read 达到 `-N 50` secondary 上限，`-N 200 -p 0.5` 重比对路径只在挑战测试中验证。
- SMT4a 的 11 个带冲突成员全部来自 primary 比对：1 个 1,360 bp tnp 位点（4 条重叠 records，
  2 条 MAPQ 60 reads 缺 900–951 bp X，对 30 个完整单位）与 7 个 81–324 bp 短片段。PCR 文库，不作多态结论。
- R-93128 的 4 个变化均为 v1 仅 1 条 read 判 full_span、审核版看到该处只有 secondary 定位（重复区域）。
- 所有支持单位标 `independence_unverified`；AR20-38 与 SMT4a 为 PCR 文库，标 `pcr_library_independence_unverified`。
  ONT-only 组装株标 `assembly_built_from_these_reads`，混合组装株标 `hybrid_assembly_used_these_reads`。
- 统一报告：[members.tsv](results/report/members.tsv)（每 record 的长度层、边界、完整侧翼、ACGT triplet、
  exact alias、部分重叠、v1/审核支持、tnp）、[families.tsv](results/report/families.tsv)、[samples.json](results/report/samples.json)。

### 人工核查

[spot_checks.tsv](results/report/spot_checks.tsv)：53 个成员（AR20-38 RS_R1 全部 11 个；每株前 5 个优先候选；
首个短片段、edge/censored、tnp-unsupported、冲突/歧义），从 assembly 重新切 X 与两侧 500 bp 与导出逐字比较
（0 不符），并核对边界依据与逐 read 证据；3 个类别 none observed。全部判定 consistent。

**核查中发现并修复的系统错误**：第一次运行把 secondary（MAPQ 0）比对上“跳过 X”的 read 计为冲突，
造成 AR20-38 相邻 IS5 拷贝（RS_R1_m6/m7，相距约 6.8 kb）及 LjRoot223_RS_R9_m1 的假冲突。
修正为：X 缺失证据只来自 primary 比对，secondary 记为 competing（commit 556b375，新增挑战测试），
五株全部重跑；第一次结果保留于集群 `attempt1_secondary_conflict_bug/`。

## 运行记录与基础设施问题

- 每阶段命令、退出码、wall、峰值 RSS：各样本 `stages.tsv`；汇总 [resources.tsv](results/report/resources.tsv)。
- 计算节点首次调用 BLAST+ 冷启动约 200 s（scratch 上的 conda 环境加载 93 个共享库），超过 doctor 固定 15 s 探测；
  T0 第一次因此失败，加入工具预热阶段后通过（未改 package 超时）。
- lr6 若干节点 `/global/scratch` 挂起或 I/O 错误：n0292（EIO）、n0082 与 n0031（写入阻塞，日志与输出均不可见，
  探针作业 2 分钟超时）、此前 n0348、n0039。这些节点均未被调度系统 drain。
  filled_site.sbatch 现先以 60 s 超时探测 scratch（失败以 75 退出），日志写到 /global/home；重跑均排除这些节点。
- T0 作业记录的 source commit 为 406dae8；作业启动约 20 s 后提交的 a8856f6 只改了 stage 日志格式、预热和两个未在 T0 使用的脚本。

## 下一批 10 株（提议，先最多 2 个并发）

[next_batch_proposal.tsv](next_batch_proposal.tsv)：Paraburkholderia 剩余合格样本 9 个——CT11-68、DD10（Complete）；
FW306-02-B、LjRoot267、R-92578（高连续性 draft）；fungorum 克隆成员 R-92583、R-92561、R-91495、R-92973
（保留各自样本身份，检查局部插入差异）。按来源与数据资格选取，未参考候选产量。
第 10 个名额固定规则：Burkholderiaceae 中 BioSample 精确匹配 Complete Genome 的 genomic ONT WGS、
环境或植物相关培养分离株、非 PCR 文库，取最小 assembly accession；需先审计 ENA 上 654 个 runs，未执行。
DD10 的 MISAG 模板问题须在运行前核实。

资源估计（按本轮实测，非保证）：每株 callout 2 CPU 约 5–20 min、<0.5 GB；filled-site 16 CPU ≤10 min、峰值 ≤8 GB；
reads 下载 0.1–2.2 GB，ENA 约 0.5–1 MB/s，为主要耗时（≤1 h/株）。2 个并发时 10 株约 5–6 h wall；
scratch 占用每株 reads ≤2.2 GB＋PAF ≤0.2 GB＋callout ≤0.1 GB，合计约 25 GB。

## 早前快照检查（保留）

PR #1 早前版本中对四株归档小表运行的 `--snapshot` 检查见 [snapshot_checks.json](snapshot_checks.json)；
它只检查表间一致性，不能替代上面的 live 与审核结果。
