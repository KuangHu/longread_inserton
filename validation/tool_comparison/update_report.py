"""Fill the checked-in report's numerical section from completed run outputs."""
import argparse,csv,json
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('root',type=Path);p.add_argument('report',type=Path);a=p.parse_args()
s=json.loads((a.root/'comparison/summary.json').read_text())
lines=['| 方法/配置 | ≥500 bp INS候选 | 其中native PASS | 有显式INS序列 | ≥3条occupied read序列支持 | DEL候选 |',
       '|---|---:|---:|---:|---:|---:|']
for k,label in [('Sniffles2','Sniffles2原候选基线'),('DeBreak_default','DeBreak默认深度过滤'),('DeBreak_no_depth_cap','DeBreak关闭深度上限（诊断）'),('TLDR','TLDR无家族库/每簇最多20 reads'),('TLDR_focal','TLDR局部诊断，不能与全范围比产量')]:
    if k not in s:continue
    v=s[k]['INS'];d=s[k]['DEL']
    lines.append(f'| {label} | {v["candidates"]} | {v["native_pass"]} | {v["sequence_present"]} | {v["three_occupied_read_support"]} | {d["candidates"] if not k.startswith("TLDR") else "不适用"} |')
lines+=['','INS的3-read序列支持不依赖empty证据，也不代表IS鉴定。TLDR数量包括非PASS候选，需同时查看native filter。','',
        '| 其他路线 | 实际完成的范围 | 本次输出 |','|---|---|---|',
        '| DIVE | 427条focal reads完成；60415条targets reads触及3 GiB限制失败 | 1个25 bp内部anchor；未恢复完整插入 |',
        '| Palidis | 原生IR检测/切片/映射完成；完整Palidis4未运行 | 探索性配对得1个790 bp结构候选，非已确认IS |','',
        '长度分层（INS；TLDR含其原生非PASS记录）：','',
        '| 方法 | 500–999 bp | 1,000–1,999 bp | 2,000–4,999 bp | 5,000–50,000 bp |','|---|---:|---:|---:|---:|']
for k in ['Sniffles2','DeBreak_default','DeBreak_no_depth_cap','TLDR']:
    if k in s:lines.append('| '+k+' | '+' | '.join(map(str,s[k]['INS']['length_bins'].values()))+' |')
run=a.root/'runs/tldr_no_detail/run.json'
if run.exists():
    r=json.loads(run.read_text());lines+=['',f'TLDR全范围重跑状态：`{r["status"]}`，退出码{r["exit_code"]}，耗时{r["wall_seconds"]:.1f}秒；时间不是受控性能benchmark。']
    t=a.root/'runs/tldr_no_detail/calls.table.txt'
    if t.exists():
        raw=list(csv.DictReader(t.open(),delimiter='\t'))
        lines.append(f'原生TLDR表共{len(raw)}条记录（含长度0及ShortIns等失败结果）；上表仅统计500–50,000 bp部分。')
text=a.report.read_text();start='<!-- RUN_SUMMARY -->';end='<!-- END_RUN_SUMMARY -->'
if end in text:
    before,rest=text.split(start,1);_,after=rest.split(end,1);text=before+start+'\n\n'+'\n'.join(lines)+'\n\n'+end+after
else:text=text.replace(start,start+'\n\n'+'\n'.join(lines)+'\n\n'+end)
a.report.write_text(text)
