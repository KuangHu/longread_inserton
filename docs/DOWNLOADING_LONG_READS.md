# 长读段下载规范与实测记录

验证日期：2026-10-21。目的：规范公共数据获取、FASTQ 转换、完整性验证和数据来源记录。
原 `nanopore_insertion_callout` 仓库保持不变。本文和辅助脚本归入独立的 insertion-site-recovery 项目。

## 实际完成的工作

8/8 成功。每个源文件取**前 500 条完整 reads**，输出独立、完整的 gzip FASTQ，共 4,000 条 reads。
土壤与匹配肠道样本分别测试 ONT 和 HiFi，因此六组数据共八个文件。
实际获取远端文件前缀共 36 MiB；重新压缩后的 FASTQ 总计约 17.18 MB。

**这是有限前缀试下载，不是八个完整 run 的下载。**
所有子集已读至 EOF，验证其自身 gzip CRC、FASTQ 四行结构、序列/质量长度和 SHA256。
原始大文件的完整 MD5 **没有**被验证；manifest 明确记录 `original_file_md5_verified=false`。
源文件名前缀不是随机抽样，因此下表不可当作全 run 的质量或读长估计。

| 数据 | Run/源文件 | 子集 reads | 子集中位读长 bp | 子集 N50 bp | 报告质量情况 |
|---|---|---:|---:|---:|---|
| 海水 CCS | SRR28395891 | 500 | 3968.5 | 6038 | 全部碱基统一 Q30，不能用于推断原始测序准确度 |
| 土壤 ONT | ERR15289757 | 500 | 5705 | 8030 | 非统一质量值 |
| 土壤 HiFi | ERR15289804 | 500 | 9711 | 11396 | 7 种离散质量值 |
| 配对肠道 ONT | ERR15285694 | 500 | 1422.5 | 6159 | 非统一质量值 |
| 配对肠道 HiFi | ERR15289675 | 500 | 9001 | 9925 | 7 种离散质量值 |
| 102 MAG 研究的肠道 HiFi | SRR17687125 | 500 | 13566.5 | 14552 | 非统一质量值 |
| 纵向肠道 ONT | ERR14838512 | 500 | 1499.5 | 8169 | 非统一质量值 |
| 堆肥 ONT SUP | PAW97819.sup_v5.0.0.bam | 500 | 6693 | 11306 | BAM→FASTQ，保留原有质量 |

质量值来自所下载文件，是报告值，**不是通过真值比对测得的准确率**。
质量值压缩/分箱和占位值可能使不同来源不可直接比较。

## 发现的两个重要元数据问题

1. ENA 为海水和 102-MAG 肠道数据生成的文件名都带 `_subreads.fastq.gz`。
   海水样本 500 条 header 均含 `CCS`，且对应论文报告 CCS；不能仅凭文件名把它归为 CLR。
   肠道数据的 HiFi 判断也依据对应研究，而不是仅凭 `PACBIO_SMRT`。
   `source.platform` 是研究支持的预期处理路线，不是仅凭 read 内容完成的技术鉴定。
2. 海水子集所有质量字符都是 Q30。后续分析应标记统一质量，不能把它当成每个碱基均具有真实 Q30 的证据。

## 下载通道规范

### A. SRA/ENA 中的原始测序 run：ENA metadata → HTTPS FASTQ

SRA、PacBio 数据目录和论文是发现入口；此批 archive 数据实际统一从 ENA HTTPS 取得。
不需要为了同一个 accession 从多个目录重复下载，也不需要先下载 `.sra` 再转 FASTQ。

每个 run 先查询 ENA `filereport`，保存：

- run、BioSample、测序平台/仪器、library strategy/source；
- 原始 `fastq_ftp`、`fastq_bytes`、`fastq_md5`；
- 元数据查询 URL 和 UTC 时间；
- 论文依据，以及预期的 ONT/HiFi 路线。

本脚本仅接受所选 run 对应一个 FASTQ 的情形。多文件 run 会显式报错，避免静默漏文件。

`probe` 通过 HTTP Range 获取源文件前缀，达到指定完整记录数后停止。
它要求 HTTP 206 和正确的 Content-Range，并检查连续请求的 ETag 一致。
源 gzip 前缀本身不是完整 gzip，**不会作为最终 FASTQ 交付**；只把完整 reads 重写成新的有效 gzip 文件。

### B. ONT Open Data：S3 basecalled BAM → FASTQ

选择 SUP basecalls，不下载庞大的 POD5，也不重复下载 HAC 和 SUP。
选中的单个完整 BAM 为 126,937,137,871 字节，约 127 GB。
这次通过 HTTP Range 获取 BAM 的前 4 MiB，以 pysam 提取 500 条完整 primary reads，转为 FASTQ。

读取截断的 BAM 前缀时，htslib 的 `no BGZF EOF marker` 警告是预期行为：
我们没有把前缀冒充完整 BAM。最终生成的 FASTQ 则经过了完整 EOF/CRC 验证。
转换使用 forward orientation，并同步处理质量方向；secondary/supplementary 记录不重复输出。
FASTQ 不保留 BAM 的修饰碱基等附加标签；本次 insertion 序列流程只需要序列和质量。

如果以后下载完整 basecalled BAM，需先按存储空间和网络预算单独安排，再用正常的 BAM→FASTQ 转换流程。
本工具的 `full` 模式不支持完整 BAM 转换，遇到 BAM 会明确拒绝，避免意外下载 127 GB。

## 复用命令

以下命令从仓库根目录执行。


Python 3.10+；只有 BAM 路线需要 `pysam`（本次为 0.24.1）。

```bash
python -m pip install pysam==0.24.1

# 八个文件各取前 500 条完整 reads，使用全新输出目录。
python tools/downloads/download_probe.py probe --out next_probe --reads 500 --workers 3

# 只检查同一土壤样本的两条路线。
python tools/downloads/download_probe.py probe --out soil_probe \
  --source soil_ont soil_hifi --reads 500

# 下载一个完整 archive FASTQ，并验证 metadata 给出的大小和 MD5。
# 此命令为后续完整下载提供；本次没有对远端完整文件运行它。
python tools/downloads/download_probe.py full --out complete_marine \
  --source marine_ccs --max-full-gb 1
```

完整 FASTQ 下载默认限制为 10 GB，并检查剩余磁盘空间；参数单位为十进制 GB。
下载使用 `.part` 文件，只有大小和 MD5 都匹配后才发布正式文件；失败保留 `.part`，不会宣称成功。
当前完整下载不做断点续传；再次运行会重新下载 `.part`，不会复用未知来源的部分文件。
原始文件 MD5 检验用于传输完整性，不构成真实性或生物学质量鉴定。

`probe` 默认每次获取 4 MiB，最多获取 32 MiB；不足以提取指定 reads 时失败，不静默输出更少的 reads。
Range 获取有最多三次即时尝试，元数据获取失败会记录为 failed。
输出目录必须为空，不自动复用旧结果。

## 配套文件

仓库包含 `tools/downloads/download_probe.py` 和 `tools/downloads/test_download_probe.py`。
以下 `results/` 和 `download.log` 是运行产物，未提交到 Git；实测样本另存于已交付的 `longread_download_pilot.zip`。


- `download_probe.py`：可复用下载/转换/质量摘要程序。
- `test_download_probe.py`：5 项完整性回归测试，已通过；涵盖 Range 被忽略、FASTQ 截断/长度不匹配、完整文件 MD5 错误时不得发布。
- `results/summary.json`：8 个源文件的完整记录。
- `results/<dataset>/manifest.json`：源 URL、元数据、传输范围、输出 SHA256、质量摘要。
- `results/<dataset>/*.fastq.gz`：真实 reads 的 500-read 子集。
- `download.log`：实际下载中的 stderr（包括预期的截断 BAM 警告）。

测试完整下载的 MD5 发布逻辑使用本地受控数据；不等于已下载完整远端 run。

## 本轮没有验证的内容

尚未配对下载 assembly/MAG，也没有进行真实 insertion calling。
500 条 reads 分散在复杂群落中，通常不足以构成目标位点覆盖，因此不能把此子集跑出零事件解释为数据没有 insertion。
下一步应选择一个明确匹配 assembly 的完整样本，确定目标 MAG 覆盖度，再评估 occupied/empty 两种状态的检出。

## 数据来源

- 海水 CCS：<https://doi.org/10.1186/s40793-025-00773-8>
- 同一样本土壤/肠道双平台：<https://doi.org/10.1038/s41467-026-69760-y>
- 肠道 HiFi/102 MAG：<https://doi.org/10.1038/s41467-022-34149-0>
- 纵向肠道 ONT：<https://doi.org/10.1038/s41586-025-09786-2>
- 堆肥：<https://epi2me.nanoporetech.com/compost_mgx_2026.04/>

这里只验证上述选定源文件的访问和转换，并非验证 ENA、SRA、MGnify、JGI 等全部平台的所有下载方式。

## 本地测试

```bash
(cd tools/downloads && python -m unittest -q test_download_probe)
```

