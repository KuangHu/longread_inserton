# 相邻重复位置的竞争比对：扩容前补充复核

2026-10-23；检查的原始提交：`b772dc75e80e7926eaaff02a4debeba0aa008177`。
原有 43 个测试本地全部复现通过；补充一个反例后发现此前测试未覆盖的错误。

## 可复现错误

`audit_filled_site.py::competing()` 原来在同一 contig、同向且比对起点相距不超过
`candidate_length + 200` 时，直接视为 same_place，跳过竞争判断。
相邻插入拷贝可以满足这个条件，但仍是同一段 read 的不同定位。

反例为纯计算输入：X=[5000,6200)，长度 1200 bp；read r1 的 primary 为
[4000,7500)，secondary 为同一 contig 上 [5300,8800)，均覆盖相同 query 区间、AS=7000。
r2 使用稍不同的长度/端点，避免被支持单位去重。两个替代位置均偏移 1300 bp，
落入原来 1400 bp 的免检查范围。

修复前：`audited_full_span`、2 full units、0 competing，独立汇总 checker 也通过。
应为：`ambiguous_competing`，因为两条 reads 都有等分的不同位置。
checker 原本只是从已发出的 read 状态复算计数/分类，不会从原始 PAF 独立复算竞争位置，
因此不能用它的 PASS 排除这类上游分类错误。

## 修复与验证

- 仅在 target、strand、target start/end、query start/end 和 CIGAR **全部相同**时，
  才把另一条记录视为重复的同一路径。附近但偏移的定位必须接受原 AS/query-overlap 竞争检查。
- 不改 AS ≥95%、query overlap ≥50%、CIGAR 覆盖、支持数、原生发现或边界规则。
- 新增同 contig 相邻等分定位反例，以及完全重复路径不误判为新竞争位置的对照。
- schema 格式仍为 audited_filled_site_v1；新结果明确记录
  `classifier_revision=nearby_competition_v2`，checker 报告也携带该版本。
- 本地 `PYTHONPATH=src python -m unittest discover -s tests -v`：**45 tests PASS**。

修复后的五株真实数据尚未重算，不能据此断言有多少真实候选受影响。
原始 read PAF 在集群；仓库逐 read 证据表不足以恢复所有跨位置竞争比对。

## 修复前产量应如何读

这里的 priority 是 X≥400 bp、双侧各500 bp完整、triplet 仅 ACGT、边界为 deposited/provisional。
下面从已归档 members.tsv 重新统计，支持数为**旧审核器**结果，不能充当修复后验证。

| 样本 | priority records | 旧审核支持 records | 其中 conflict | 支持 records 的完全相同坐标归组 |
|---|---:|---:|---:|---:|
| PaAR2038 | 520 | 520 | 0 | 420 |
| SMT4a | 199 | 199 | 4 | 134 |
| LjRoot223 | 146 | 146 | 0 | 106 |
| R-93128 | 174 | 161 | 0 | 128 |
| R-91497 | 120 | 120 | 0 | 102 |

这些计数有跨软件部分重叠，不能解释为独立 insertion 数。
多数边界仍为 provisional_homology_median；read full-span 不证明边界精确或转座活性。
对外用“候选结构获 read 支持”，不用“fully confirmed insertion”。

## Cluster agent 的下一步（限定范围）

1. 取本修复提交，保持五株 assembly/native/boundary/reads 不变；保存旧结果，不覆盖。
2. 用已有完整 PAF 和相同 reference/member_boundaries，向新的审核输出目录运行 classify。
   确认 PAF 含原 secondary/AS/cg 信息且 hash 与原记录一致。若输入完整无变，无需重新下载、
   原生发现或重新 mapping；达到 secondary 上限时仍执行既有 recheck 路径。
3. 五株运行 check_audited_evidence，更新统一 report、resources 与 release_check。
   每株 summary/checker 必须记录 nearby_competition_v2；不能仅修改旧 JSON 的版本字段。
4. 交付所有改变状态/支持数的成员差异表，人工核查全部改变分类的成员及 AR20-38 RS_R1。
   若无真实成员改变，也保存零差异表，不能省略重算。
5. 全部检查通过、变化有解释后恢复 ready_for_10_sample_batch=true。
   此修复不改变 datatype、发现路线或生物学产量门槛，不引入 empty-site 验证依赖。

不需要重开 pipeline 设计，也不需要重新运行 PILER/RepeatScout；只修复并复核 read 竞争判定。
