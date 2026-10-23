# Call-out scale v1：固定验收范围与本次检查

主规范：[CALLOUT_SCALE_PROTOCOL.md](../../docs/CALLOUT_SCALE_PROTOCOL.md)。
集群执行交接：[CLUSTER_TASK.md](CLUSTER_TASK.md)。
本目录是单株候选导出与 read 证据的扩容验收，不是新的 discovery 算法。

## 本次实际执行

- 从 PR #1 的 `ee229861df8d27b6c265ca04c3693662394daa83` 取得源码及归档小结果。
- Linux / Python 3.12，本地安装 multisite-insertion 0.1.0、Biopython 1.88、NumPy 2.5.4。
- `python -m unittest discover -s tests -v`：26 个测试通过；原有 17 个、新增 9 个。
  新检查覆盖 coordinate/evidence 版本错配、缺行/重复行、错误 subset、reads 样本错配、
  censored/N-containing triplet 不升级、AR20-38 IS5 小表回归，以及 live 合法零命中和输入 hash 篡改。
- [snapshot_checks.json](snapshot_checks.json)：四株已归档结果全部通过小表一致性检查。
  该检查包括成员 ID、BED/read 坐标、家族归属、长度、annotation 标签、汇总计数和下载校验记录。

| 样本 | 导出 records | ≥400 bp records | 其中 v1 full≥2 records |
|---|---:|---:|---:|
| PaAR2038 | 762 | 526 | 526 |
| DdD110 | 164 | 64 | 64 |
| AmHOT1A3 | 188 | 134 | 134 |
| RoMoAcy1 | 646 | 360 | 356 |

以上都是 records，含跨工具/边界重复，不是独立 insertion 或背景数。
后三株 PacBio 只用于已有格式回归；不加入本轮 Nanopore 五株发布验收面板。
PaAR2038_RS_R1 的 11 条、1171–1173 bp、v1 full 145–181 已从表核对。
这不代表重新验证了精确边界、tnp CDS 包含关系或 55 对侧翼比对。

## 尚未完成

- 本地没有 pals/piler2/RepeatScout/blastn/minimap2，未执行新的真实原生 discovery 或 mapping。
- 没有集群 round-1 的原始 run/PAF/FASTQ；未独立复核四株新结果。
- 本次没有实现审核版 CIGAR/竞争位置检查，也没有运行跨株 empty-site 验证。
- 因此 [release_check.json](release_check.json) 的 `ready_for_10_sample_batch` 为 false。
  它是明确的当前验收状态，不是对已有生物学结果的否定。

测试清单固定在 [test_panel.tsv](test_panel.tsv)。不要把它直接传给旧 samples.tsv 解析器，
两者字段不同。所有重跑使用新输出目录，保留当前 v1 和 round-1。

小结果 checker 可用 `--snapshot` 在本机复现；完整 run 用默认 live 模式。
具体命令、输入布局、检查范围和扩容通过条件见主规范。
