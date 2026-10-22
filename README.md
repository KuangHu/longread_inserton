# multisite-insertion

从**组装 FASTA** 中寻找“内部相似、外部不同”的多位点重复区段，并导出候选
`upstream + insertion + downstream`。无需 perfect empty site。

这是 `KuangHu/longread_inserton` 的新路线。旧的 reads → SV → 双状态恢复原型已退出
当前代码和安装入口；历史保留在
[旧版本 5ad842b](https://github.com/KuangHu/longread_inserton/tree/5ad842b0c8261485c27fa685de9fc538ff49cd23)。
仓库名称沿用，但本 package **不读取 FASTQ/BAM、不执行长读段 SV calling**。

## 固定流程

1. 每个样本独立运行 **PALS → PILER-DF** 和维护版 **RepeatScout**，使用原生默认发现参数。
2. 原样保存所有软件输出，解析全部家族和成员；不增加长度、拷贝数、复杂度或注释淘汰规则。
3. 发现后添加已有 CDS 注释和 BLASTX 蛋白证据，保留全部家族的注释表。
4. 默认对 **tnp-supported 家族**做上下文比对、候选边界推断和两侧各 500 bp 提取。
   `--subset all` 可分析全部原生家族；这不会改变原生发现结果。
5. 输出同家族对应侧翼比较、跨软件坐标重叠和每条记录的依据/未解决状态。

**原生发现与我们的后处理明确分开。** 边界推断不是 PILER/RepeatScout 的精确插入判定。
tnp-supported 是家族层面的支持，不要求每个成员都有完整 tnp；也不证明转座活性。
没有 TIR/TSD 必需条件，没有空位点门槛，不合并跨软件记录为“独立事件数”。

## 安装

Linux/WSL，Python ≥3.10。先安装原生依赖，见 [INSTALL.md](docs/INSTALL.md)。

```bash
git clone https://github.com/KuangHu/longread_inserton.git
cd longread_inserton
python -m pip install .
multisite-insertion doctor
```

Python package 只依赖 Biopython。PALS、PILER、MUSCLE 3、RepeatScout 和 NCBI BLAST+
是独立可执行程序；`pip install` 不会安装它们。未向 PyPI 发布。

## 运行

创建 `samples.json`，相对路径相对于 manifest 所在目录：

```json
{
  "schema_version": 1,
  "samples": [
    {
      "id": "isolate1",
      "fasta": "data/isolate1.fasta",
      "annotations": ["data/chromosome.embl", "data/plasmid.embl"],
      "annotation_format": "embl"
    }
  ]
}
```

一个样本的染色体和质粒放在同一个 FASTA 中；不同 isolate 分为不同 sample。
**不要把近缘菌的 MAG 混成一个 sample**：跨 MAG 的保守宿主基因也会被发现为重复。
注释文件可为 EMBL 或 GenBank；每个已提供注释的序列会与 FASTA 逐碱基核对。
注释 record ID 必须对应 FASTA ID 或其中竖线分隔的 accession。

```bash
multisite-insertion run --manifest samples.json --out runs/isolate1
multisite-insertion validate --run runs/isolate1
```

只看软件原始发现结果：

```bash
multisite-insertion discover --manifest samples.json --out runs/native_only
```

没有注释也可以发现和推断：省略 `annotations` 并使用 `run --subset all`。
此时不会产生 tnp 支持；默认 `tnp-supported` 子集可能为空，但全部原生家族仍保留。
两种方法可通过 `--methods PILER` 或 `--methods RS` 单独运行。
所有命令要求新目录或空目录，避免旧结果混入；失败记录为 `failed`，不冒充零命中。

## 复用已有原生输出

在 sample 中添加 `native`；不指定的方法才会重新运行。导入目录应对应**同一份 FASTA**。
程序核对 contig/坐标并记录导入文件哈希，不能仅凭坐标证明旧输出来自哪一份序列。

```json
"native": {
  "PILER": "previous/isolate1_PILER_default",
  "RS": "previous/isolate1_RS_default"
}
```

PILER 目录需含 `families.gff` 和 `library.fasta`；RS 目录需含 `ranges.tsv` 和
`families.fasta`。完整目录原样复制到新运行的 `native/`。继续使用同一 `run` 命令。

## 主要输出

| 路径 | 内容 |
|---|---|
| `native/<sample>/<method>/` | 原生输出和可追踪的导入记录；后处理不修改 |
| `results/all_native_families.json` | 所有原生家族、成员、有效性与已有注释重叠 |
| `results/all_families_annotated.json` | 全家族的蛋白证据和 tnp-supported 标记 |
| `results/member_boundaries.json` | 可提取成员的原生、推断、导出坐标与证据 |
| `results/record_inventory.json` | 所选子集全部成员，包括未导出序列的越界记录 |
| `results/invalid_native_records.json` | 越界原生范围原样保留，不裁剪成假序列 |
| `sequences/candidate_insertion.fasta` | 候选区段；必须结合头部的边界状态解释 |
| `sequences/upstream_500.fasta`、`downstream_500.fasta` | 同 member ID 的两侧序列；contig 末端可不足 500 bp |
| `sequences/triplet.fasta` | 三部分按家族方向拼接；负链已反向互补 |
| `sequences/genomic_*.fasta` | 基因组正向的左侧、核心、右侧版本 |
| `results/family_summary.*`、`flank_pair_comparison.json` | 全部可提取家族及逐对侧翼证据 |
| `results/cross_tool_locus_overlaps.json` | 所选家族成员的跨软件原生范围重叠，不自动去重 |
| `results/candidates.bed`、`member_coordinates.txt` | 逐条坐标；BED 为 0-based 半开区间 |
| `run.json`、`provenance.json`、`native_manifest.json` | 状态、命令、输入/二进制/原生文件哈希 |

FASTA 头部候选范围是 **1-based 闭区间**；JSON 的 `start0/end0` 为 **0-based 半开区间**。
up/down 采用原生家族方向，**不一定是 tnp 转录方向**。空侧翼仍保留 FASTA 条目，BLAST
使用独立的非空侧翼输入。无显著侧翼匹配不等于已证明侧翼完全不相关。

边界状态、算法和局限见 [METHOD.md](docs/METHOD.md)。

## 已做的验证

三个环境分离株（AltDE、UT26S、RHA1）的既有原生输出用于回归验证：

| 方法 | 全部原生家族 | tnp 支持家族 | 所选原生记录 |
|---|---:|---:|---:|
| PILER | 23 | 17 | 72 |
| RepeatScout | 238 | 52 | 228 |

打包后 **300 条坐标记录、8 份 FASTA 与此前结果一致**；家族、侧翼及跨软件重叠统计一致。
这是软件重构回归，不是独立的生物学准确率评估。详情和可复现检查见
[benchmarks/README.md](benchmarks/README.md)。

单元测试与无需下载天然基因组的合成运行示例：

```bash
python -m unittest discover -s tests -v
python examples/make_demo.py --out data/demo
multisite-insertion run --manifest data/demo/manifest.json --out runs/demo
multisite-insertion validate --run runs/demo
```

合成示例只有随机重复序列和人为测试注释，用于检查软件连通性，不能充当天然 IS 阳性证据。

## 范围与限制

本版本固定为这条 assembly-based 路线，不加载旧 `site_recovery`、Arm A、Arm B、图气泡或
empty/filled 恢复模块。完整元件边界、独立事件数、历史插入方向和转座活性均需另行验证。
高度分化的重复、低拷贝元件、组装折叠及断裂仍可能漏检；短/碎片化 RS 输出仍可能只是片段。
当前 ±2000 bp 上下文和 500 bp 侧翼为 pilot 固定设置，**不跨 circular origin 绕回**。
