# 主线纠偏：相似 core、不同上下文

2026-10-27。此文件优先于旧扩量计划中的 gap-first 发现顺序。

## 目标与当前状态

```text
Assembly A: upstream A -- core X  -- downstream A
Assembly B: upstream B -- core X' -- downstream B
Assembly C: upstream C -- core X''-- downstream C
```

core 同源或相同，左右上下文代表不同候选位置；不要求对应空位点。
主要目标仍是完整 core <2500 bp、1–2 完整编码单元。一个样本一个拷贝也允许参与跨样本发现。
没有 reads、没有 gap、没有已知功能注释都不是发现入口的排除条件。
gap 的用途是辅助推断 core 边界、解释 presence/absence、提供额外读段检查对象。

历史 Stage 1 / Stage 2 已完成的 gap-derived 目录和 read 检查继续保存，作为辅助分支。
历史 `assembly_catalog_ready_for_long_reads` 只描述 gap 目录交付，**不代表本主线已完成**。
现有 16,414 main loci 和 26,643 蛋白家族不能作为本主线全库发现的产出。

2026-10-27 已增加可运行的 `experiments/core_context/core_context.py`，直接从联合 FNA 的重复 graph blocks 取实例，
不使用 gap 清单。真实 H16/H850 smoke 重新构图并通过 verify，找回旧 R33 六位置（近似边界）。
这是入口修复和真实数据回归；GFF3 编码复核、蛋白同源独立发现、精细上下文分组和全库扩量仍待完成。

## Agent 执行顺序

1. 阅读新入口 README、smoke 汇总和测试，复现核心反例：所有样本都带 core，侧翼不同，没有任何空位点，仍能产生实例。
   相同 core/相同上下文必须作为对照；反向、多拷贝、末端、低复杂度、环状跨原点和远缘宿主需有明确状态。
2. 复用已核验 FNA/GFF3/FAA，不重复下载；按实测距离组织小型可比子集。原始 contig 和 sample 映射不得丢失。
3. 保留两条独立于 gap 的候选入口：
   - 联合图重复 blocks/经验证的 block 路径；
   - 全部 FAA/GFF3 中的同源蛋白与邻接结构（不是仅旧 gap 候选蛋白），单基因及双基因分别组织。
   对照已有 PILER/RepeatScout 单株结果可作为附加来源，不将单株多拷贝设为准入条件。
4. 蛋白入口复用 MMseqs2 的最终代表成员审计，保留两个蛋白的顺序、方向和间隔。
   一个共同蛋白不等于整个双基因 insertion 相同。蛋白家族、DNA core、结构亚型分开。
5. 将 core 实例直接定位回全部被选 FNA，保留所有合格位置和多个拷贝；不得只保留最佳命中。
   搜索范围、参数、未覆盖的样本和命中上限要明示。远缘同源发现和近缘 DNA 复核分开记录。
6. 接入已有 GFF3/coding 单元检查（strand/phase/partial/pseudo/transl_table/transl_except/circular）；
   编码范围与完整 core 边界分开。用跨实例保守区和上下文转折提出边界，gap 只补充边界证据。
   不确定、缺注释、边缘截短者保留待审查，不能进入 1–2 complete CDS 主清单。
7. 统一 core 方向，分别分析 U/D，并保留扩展窗口与邻近宿主基因。
   同一正交位置重复观察、不同上下文候选、单侧不同、无法判断分别报告。
   无侧翼比对不是不同位点的充分证据；对邻近基因同源/共线关系复核后再给 site groups，仍不叫独立事件。
8. 出具 U/core/D 案例图和可追溯表，至少有相似 core/不同两侧、相似 core/相同两侧、单株多拷贝、
   单侧争议及无 gap 来源的实例。所有种子来源和辅助 gap 链接分别存储。
9. 原型审阅后才按 batch 扩量。读段阶段检查各实例 U–core–D 连续性，随后检索新 read contexts；
   空位点验证可选，不阻塞主线。不要为了继续旧 gap 主线而启动 500 GB 下载。

## 本轮验收与停止点

- 新入口实际运行且没有 gap/empty-site 依赖；实例保留来源、方向、序列和坐标。
- 真实 smoke 和反例测试可复现；边界、annotation、orthology 状态诚实呈现。
- 后续达到 1–2 gene 与上下文审查要求前，不标记新主线生产就绪。
- 计划、代码、测试、参数、输入哈希、结果表均提交到现有 PR；历史结果不删除、不改写为新方法结果。
- 目前停止状态为 `core_context_graph_prototype_smoke_passed`，不是全量发现完成。
