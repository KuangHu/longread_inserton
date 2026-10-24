# Cluster agent：收敛 insertion call-out，为批量运行做验收

用户当前优先级是 pipeline / datatype / test run / check result，然后 scale up。
执行规范：[CALLOUT_SCALE_PROTOCOL.md](../../docs/CALLOUT_SCALE_PROTOCOL.md)。
本任务优先于继续扩展跨株 empty-site、低频非参考插入或寻找更多候选物种。

1. 检查工作树，保留尚未提交的 round-1 成果。提交 curl 修复、比较脚本/job、四株
   package summaries/member_boundaries/validation、read evidence、输入/工具 hash、失败重跑记录和中文说明。
   重复记录与真实 locus 分开；此前汇报的 12/2 是跨株结果，不能替代单株产量。
   raw reads、PAF/BAM、大中间文件不进入 Git。
2. 固定 `test_panel.tsv` 五株，核对实际 FNA/注释/reads、BioSample、版本和 library metadata。
   复用已校验大文件；不重下同一版本、不混用 GCA/GCF 序列、不改工程分组来追阳性。
   test_panel 是验收清单，**不是旧 common.sh 的 samples.tsv 格式**；从已核实的旧 samples 行
   生成执行清单，保存映射。不得直接传给旧 sbatch 脚本。
3. 不改 PILER/RS 默认与边界逻辑，明确使用 --subset all。
   将只做 filled-site read check 的安全 wrapper 从旧 03_reads.sbatch 拆出，要求新输出目录、
   显式输入、set -euo pipefail、完整命令和退出码；不能 rm -rf 旧 evidence。
   当前轮无需非参考插入提取/回贴，也无需另一株 empty-site reads。
4. 实现协议第 4 节的审核版 read 检查，另存 evidence_audited/，保留 v1 逐成员差异。
   必须有逐 read ID、原始/派生关系状态、CIGAR 数值、竞争比对、冲突与失败原因；
   达 secondary 上限的候选补查，不能只看 MAPQ。先写明确真值的挑战测试。
5. 运行 package tests、原生 demo、AR20-38 逐成员回归，再固定五株 filled-site 检查。
   对每株现有 v1 输出运行 `scripts/check_callout_bundle.py` live 模式，归档报告。
   该脚本仅检查已有 v1 表与 exports，不会自动读取/认证新的 audited 输出。
   实现审核版时增加独立 schema/checker，不能把 audited 结果伪装成 v1。
6. 生成统一每株/每家族报告：所有 native records；<400/≥400；边界状态；完整两侧数；
   无模糊碱基 triplet 数；v1/审核支持；tnp-supported/unsupported/unavailable；
   exact-coordinate aliases、非完全重合歧义数。不要将跨工具支持相加或称历史事件数。
   X、up、down、triplet FASTA 以同 member ID 关联；no-hit、失败、unresolved 均不能丢。
7. 按协议固定规则人工检查并记录 member ID、assembly 切片、边界证据、reads 与结论。
   报告 wall time、CPU/RSS/disk 和每个阶段退出码。若五株合格，提出固定的下一批 10 株清单
   与资源估算，先最多 2 个并发；不要直接发起全库/百株任务。
8. 输出 release_check.json，逐项列 T0/T1 的 pass/fail/not_run、证据文件路径和阻塞原因。
   仅全部必要技术项通过才能标 ready_for_10_sample_batch=true。
   非参考事件数、空位点验证和候选产量不作为 release gate。
   推到 PR #1，附本轮已完成与未完成事项；不合并 main。

不要再次重新设计路线或要求用户重复确认已授权的核查/修复/提交。
任务完成条件是可复现的五株验收包和明确的扩容判定，而不是又一份只有 future tense 的计划。
