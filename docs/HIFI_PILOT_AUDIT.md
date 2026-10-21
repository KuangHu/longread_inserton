# 真实 CCS pilot 审计：面向完整插入元件的有效产出

审计日期：2026-10-21。原始 calling、恢复结果和阈值均未修改；本轮只读审计全部 276 条 INS/DEL 候选记录，并检查几个具体失败案例。

## 结论

**此前“恢复 20 个结果”不能作为 IS 元件挖掘可用性的依据。原结果中，≥1 kb 且有双状态 read 支持的只有 3 个，≥2 kb 的为 0。3 个严格通过结果仅 139、377、412 bp，对本项目寻找完整编码型插入元件的目标没有实用说服力。**

这里的 500 bp、1 kb、2 kb 是审计分层，不是 IS 的分类定义。长度合适不等于 IS；短片段也不等于不存在任何移动元件来源。尚未做功能、结构或移动性鉴定。

## 候选和恢复结果的长度分布

| VCF 报告长度 | INS/DEL 候选记录 | 恢复序列表示 | 恢复且具备双状态支持 |
|---|---:|---:|---:|
| 50–499 bp | 114 | 14 | 12 |
| 500–999 bp | 35 | 1 | 1 |
| 1,000–1,999 bp | 52 | 4 | 3 |
| 2,000–4,999 bp | 57 | 1 | 0 |
| 5,000–9,999 bp | 18 | 0 | 0 |
| 合计 | 276 | 20 | 16 |

候选记录不是独立 IS 数量：同一区域可能存在多个等位型、重叠候选或复杂结构。

≥1 kb 有 127 条记录（69 INS、58 DEL）；仅 5 条恢复了序列，全部为重复等价断点的 ambiguous 分类。其余 122 条：111 条无合格完整 CIGAR 加双侧锚定；10 条同一精确断点支持不足；1 条插入序列不一致。

5 条已恢复的长结果：

| 长度 | 类型 | occupied / empty reads | 现有解释 |
|---:|---|---:|---|
| 1,025 bp | DEL | 10 / 5 | 双状态支持，断点歧义 |
| 1,120 bp | DEL | 24 / 7 | 双状态支持，断点歧义 |
| 1,161 bp | INS | 10 / 7 | 双状态支持，断点歧义 |
| 1,468 bp | DEL | 0 / 16 | 缺 occupied 支持 |
| 2,807 bp | DEL | 1 / 7 | occupied 支持不足 |

所以，“约 93% 候选被过滤”不等于“约 93% 真阳性被漏掉”。原先把这个比例直接称为主要漏检瓶颈，证据不够。

## 完整 CIGAR 审计

在 122 条未恢复的长候选中：

- 31 条能找到至少 3 个 primary、MAPQ ≥20、长度相容、距报告断点 ≤100 bp 的完整 CIGAR gap。
- 加上当前 100 bp 双侧锚定及每侧 identity ≥95% 后，只剩 3 条仍有至少 3 个支持 read 名。
- 仅把 identity 门槛改为 90% 做诊断计数，也只增至 6 条。该计数没有完成精确断点分组、另一状态验证或 consensus；不是新增有效结果。

因此，问题并非全是“没有读取 supplementary alignment”，也不能靠统一降低 identity 门槛解决。部分 flanks 存在明显分歧或额外 gap；要区分菌株差异、复杂等位型和错误比对，需要按局部序列重新构建证据。

此审计阶段允许统计多个相容 gap；原恢复器还要求一条 read 只有一个合格命中。因此这些计数是诊断证据量，不是原调用器结果的替代。

## 精查案例

### 1. 约 1.95 kb：明确值得修复的精确分组问题

`Sniffles2.INS.3SB9`，Prochlorococcus-A-LRMED-150，contig C11。

- 3 条 read 均通过原来的双侧锚定要求。
- 插入长度：1,950、1,946、1,953 bp。
- 断点：95,195、95,193、95,195（0-based interbase）。
- 原恢复器按精确坐标拆成 2+1，导致最大组不足 3 条。
- 三段插入序列全局两两 identity：96.88%、97.29%、97.03%。
- 在名义断点 ±100 bp 范围发现 38 条局部 reference-like reads。

这是最清楚的边界归一化/等位型聚类改进机会。但仍未正式合并断点、重建并验证完整两种等位型，也未鉴定为 IS；本轮不新增 accepted call。

### 2. 约 3.1 kb：不能把邻近断点直接合并

`Sniffles2.INS.FSB9`，同一 MAG，contig C11。

- 4 条 read 可通过两侧锚定，断点仅相差数 bp。
- 插入段长度分别 2,393、2,478、2,991、3,109 bp。
- 插入序列两两全局 identity 只有 52.2%–54.0%。

这不是把精确断点改成 ±5 bp 就能可靠恢复的单一插入。它至少提示非常不同的插入序列/等位型或复杂比对，不能合成一个混合 consensus。审计不能单独判定是群落真实多态还是错误归属。

### 3. 约 1.27 kb：序列不一致的过滤有实质理由

`Sniffles2.INS.10SBD`，Prochlorococcus-A-LRMED-150，contig C4。

6 个通过锚定的 gap 长度为 1,263–1,382 bp。两两插入序列 identity 约 82.3%–99.1%；虽然有一对相似度约 99.1%，不能据此认定已有至少 3 条支持的统一插入等位型。这个失败不能简单归咎于“过于保守”；需要等位型聚类。

### 4. 约 2.05 kb：split 支持不等于可以救回完整配对

`Sniffles2.DEL.4SB6`，Pelagibacter-LRMED-234，contig C1。

- 6 个 caller-support read 名具有方向、近似位置与长度相容的 split 几何证据。
- 其中 3 条共享参考区间 [1311024,1313074)，两段 query 有 6 bp 重叠。
- 对跨接头实际 query 片段重新检查：200 bp 的参考空位点表示对应 194 bp query，identity 为 91.5%；左锚点 91%，右锚点 98%。
- 按原标准只找到 1 条 occupied reference-like read。

这不是一个仅增加 supplementary 支持就能通过原标准的长事件。仍有参考分歧、断点重叠与另一状态证据不足的问题。

## split-alignment 审计的总体结果

122 条未恢复长候选中，24 条有至少 3 个 caller-support read 名满足同 contig、同方向、MAPQ ≥20、双段长度 ≥100 bp，且 split 几何与候选大致相容的条件。

这只是诊断标准：断点容差 500 bp，长度容差与原流程一致（25%，至少 5 bp），允许少量 query/reference overlap。没有据此提取或新增 accepted 序列。

24 条中，仅 1 条有至少 3 个 read 的“另一侧未对齐 gap”绝对值 ≤10 bp；即上述 2.05 kb 例子。其余很多同时在 reference 和 query 上留有未对齐区段。它们需要局部重比对，才能区分纯插入/删除、断点修剪、替换或复杂结构；不能直接把两个参考 flank 拼接成已证实的真实 empty allele。未对齐片段本身也不足以证明发生了替换。

## 对后续工作的判断

当前原型不能直接作为完整 IS 挖掘工具。应先针对有证据的失效方式做有限改造：

1. 用上述 1.95 kb 案例验证边界归一化及序列一致性约束下的 read 分组。
2. 对 split 候选做局部双等位型重建；必须允许得出“这是复杂替换，不符合简单 insertion/empty 模型”的结论。
3. 保留多种插入等位型，避免把不相似序列合成 consensus。
4. 产出指标改为按长度分层的、两种状态得到支持的完整序列对；短 indel 不再作为本项目有效产出的主指标。

现在没有证据支持“修一下就能得到上百个 IS”。不应继续按本轮海水数据和现有流程扩大下载。即使修复明确漏掉的案例，也需要在来源匹配更明确、目标菌覆盖与读长更合适的数据上验证产出；本轮参考仍是来源未逐 MAG 核实的同研究 MAG 池。

## 审计文件和复现

[逐候选证据表](../validation/real_hifi_pilot/audit/results/candidate_audit.tsv)；[长度统计](../validation/real_hifi_pilot/audit/results/audit_summary.json)；[案例序列一致性](../validation/real_hifi_pilot/audit/results/example_cases.json)；[split 接头复核](../validation/real_hifi_pilot/audit/results/split_example_local_check.json)。

在原 pilot 数据目录仍可用、原 pipeline 依赖已安装的情况下：

```bash
python validation/real_hifi_pilot/audit/audit_candidates.py /path/to/hifi_marine_pilot
python validation/real_hifi_pilot/audit/examine_cases.py /path/to/hifi_marine_pilot
python validation/real_hifi_pilot/audit/check_split_example.py /path/to/hifi_marine_pilot
```

原始 BAM、VCF、records.jsonl 未修改。详细逐 read 审计 JSON 保留在本地 `hifi_marine_pilot/audit/`，可用以上脚本重建；Git 保存脚本、逐候选计数和案例证据。
