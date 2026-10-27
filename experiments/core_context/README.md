# Core 相似、上下文不同：恢复独立发现入口

2026-10-27。此目录是纠偏后的 assembly 发现原型。它直接读取联合 FNA 构建的 PanGraph，
枚举重复 block 的每个位置，再比较定向 U/core/D；**不读取 gap calls，不要求空位点，不要求 FASTQ。**
它没有替换或删除已有 gap、S4 或 read-check 结果。

## 找回的旧 smoke test

本地 `cross_assembly_pilot`（2026-10-24）保留了 PanGraph 1.4.0、minimap2 2.31、H16/H850 FNA、
Nanopore 前缀及已审阅的 R33 六位置对照。旧 pilot 的 70 个 seeds 来自 gap，
但它还实现过 seed → reads → assembly 上下文扩展；后来的 fna_pool 规模化流程只接入了 gap 发现和已知位点验证。
PanGraph 当时只建立并验证了图，没有用图的重复位置作为 caller 入口。原单株 package 的 PILER/RepeatScout 路线仍保留。

本次重新构建 PanGraph 并执行 `--verify`，仅用已有 H16/H850 联合 FNA：

| 检查 | 结果 |
|---|---:|
| 输入样本 / contigs | 2 / 7 |
| 50 <= block 实例长度 <2500 的重复 block groups | 264 |
| 实例 / 实例对 | 766 / 1,413 |
| 低侧翼相似性且左右邻接 block 都不同的候选对 | 348 |
| R33 对应 graph block 实例 | 8 |
| 旧 smoke 的六个位置找回 | 6/6，边界容差 20 bp，非精确边界验证 |
| 新增下载 / 使用 gap 清单 / 使用 reads | 0 / 0 / 0 |

R33 只在发现完成后用于回归核对，没有输入发现程序。八个实例中七个有完整 500 bp 侧翼，
另一个位于 contig 末端，core 也可能截短，保留为 edge-limited。额外两个实例不能称为新验证的 sites。
完整侧翼的七个实例产生 21 个不同邻接/低侧翼相似性的诊断候选对，其正交关系尚未审查。
小表、汇总和示意图见 `smoke_h16_h850/`。大 FNA、图 JSON 和实例序列未上传。

## 执行入口

依赖 Python >=3.10、Biopython；构图工具在本次真实测试中为 PanGraph 1.4.0。
准备 TSV，列为 `sample`、`fasta`。每个样本的 FNA 保留独立 contigs，不串接碱基。
相对路径相对于 manifest 所在目录；每次输出到新目录，避免误用残留结果。

```bash
python experiments/core_context/core_context.py build --manifest samples.tsv --out graph_bundle --pangraph /path/to/pangraph
python experiments/core_context/core_context.py call --bundle graph_bundle --out context_calls
python -m unittest discover -s experiments/core_context -p 'test_*.py' -v
```

重新运行已经存在的本地 pilot（无需网络）：

```bash
python experiments/core_context/reproduce_local_smoke.py --pilot /path/to/cross_assembly_pilot --out new_smoke_output
```

`build` 写出 joint FNA、contig 来源映射、图和带哈希/版本/命令的 provenance。
`call` 校验来源后枚举重复 block；保留所有合格长度实例，包括同一 assembly 内多拷贝。
输出 `occurrences.tsv/json`、`contexts.fna`、`context_pairs.tsv`、`families.tsv`、排除项和汇总。
输出均为 0-based 半开坐标；展示图使用 1-based inclusive。core 及上下游统一到 graph block 方向。

比较规则是未校准的原型诊断：core 全局 identity >=0.9；两侧均 >=0.9 为相似上下文候选；
两侧均 <0.7 且两个定向邻接 block 均不同，标为不同上下文候选。
**低侧翼 identity 本身不足以确认不同正交位置。** 不满足条件的保留 unresolved。
每个 group 默认最多审查 20,000 对，截断明确标记；实例不截断。本次真实 smoke 没有发生截断。

## 不得跳过的后续验收

- Graph block 不是已确认的完整 insertion；边界可能切开基因、重复或跨多个生物学单元。
- 本原型尚未接入 1–2 complete CDS、translation/GFF3 复核；所有实例 `annotation_status=not_reviewed`。
  不得把这 264 个 groups 当作合格 insertion families 或追加到旧 main 清单。
- 下一步需要接入已有 S4 的编码检查，并增加基于全部 FAA/GFF3 的蛋白同源入口，保留单基因/双基因结构。
  相邻 graph blocks 的合并/边界扩展也需单独验证，不能为满足 2500 bp 截短候选。
- Context 比较需要扩展窗口和邻近基因同源/共线证据；暂未计算 site groups。
- 环状或跨原点 graph paths 明确延后处理；未标环状不能证明 contig 线性。
- 此处没有 read validation、转座活性、独立事件、敏感性或新颖性结论。
- `production_ready=false`。先完成下述计划里的 prototype review，再扩量，不能直接对 22,759 assemblies 建一个大图。

主计划：[PLAN_CORE_CONTEXT_DISCOVERY.md](../../docs/PLAN_CORE_CONTEXT_DISCOVERY.md)。
