"""Add a family-level planning view without changing sample identities or priorities."""
import collections
import hashlib
import json
import sqlite3
from pathlib import Path

ROOT = Path(__file__).resolve().parent / 'ena_bacteria_ont_2026-10-24'
SOURCE = ROOT / 'taxonomy_batches'
OUT = ROOT / 'family_batches'


def load(name):
    return [json.loads(line) for line in (SOURCE / name).read_text(encoding='utf-8').splitlines()]


def save(name, rows):
    with (OUT / name).open('w', encoding='utf-8') as stream:
        for row in rows:
            stream.write(json.dumps(row, ensure_ascii=False) + '\n')


def esc(value):
    return str(value).replace('|', '\\|').replace('\n', ' ')


def main():
    OUT.mkdir(exist_ok=True)
    taxa = {r['original_tax_id']: r for r in load('taxon_mapping.jsonl')}
    groups = {g['group_id']: g for g in load('species_groups.jsonl')}
    old_batches = load('batch_catalog.jsonl')
    members = load('run_membership.jsonl')
    con = sqlite3.connect(ROOT / 'run_inventory.sqlite')
    sizes = dict(con.execute('SELECT run_accession,fastq_total_bytes FROM runs'))
    con.close()
    group_rows = []
    batches = {}
    for gid, group in groups.items():
        lineages = [taxa[t] for t in group['original_tax_ids']]
        families = {t['family_tax_id'] for t in lineages}
        assert len(families) == 1, ('Conflicting family lineage', gid)
        family_id = families.pop()
        family_name = lineages[0]['family_name']
        bid = f'family_{family_id}' if family_id else 'unplaced_' + group['batch_id']
        row = dict(group_id=gid, display_name=group['display_name'],
                   genus_batch_id=group['batch_id'], genus_tax_id=group['genus_tax_id'], genus_name=group['genus_name'],
                   family_batch_id=bid, family_tax_id=family_id, family_name=family_name,
                   review_priority=group['review_priority'], fastq_unique_samples=group['fastq_unique_samples'],
                   grouping_basis='NCBI family ancestor' if family_id else 'No family ancestor; retain previous group')
        group_rows.append(row)
        if bid not in batches:
            batches[bid] = dict(family_batch_id=bid, family_tax_id=family_id, family_name=family_name,
                                previous_batches=set(), genus_tax_ids=set(), genus_names=set(), group_ids=[],
                                singleton_groups=0, runs=0, fastq_runs=0, samples=set(), fastq_samples=set(),
                                known_fastq_bytes=0, low_reference_runs=0)
        b = batches[bid]
        b['previous_batches'].add(group['batch_id'])
        if group['genus_tax_id']:
            b['genus_tax_ids'].add(group['genus_tax_id'])
            b['genus_names'].add(group['genus_name'])
        b['group_ids'].append(gid)
        b['singleton_groups'] += group['fastq_unique_samples'] == 1
    lookup = {r['group_id']: r for r in group_rows}
    new_members = []
    for member in members:
        group = lookup[member['group_id']]
        b = batches[group['family_batch_id']]
        sample = member['sample_accession'] or 'RUN:' + member['run_accession']
        b['samples'].add(sample)
        b['runs'] += 1
        if member['availability'] == 'fastq_listed':
            b['fastq_runs'] += 1
            b['fastq_samples'].add(sample)
            b['known_fastq_bytes'] += int(sizes[member['run_accession']] or 0)
        b['low_reference_runs'] += group['review_priority'] == 'low_reference'
        new_members.append(dict(member, family_batch_id=b['family_batch_id'], review_priority=group['review_priority']))
    for b in batches.values():
        for field in ['previous_batches', 'genus_tax_ids', 'genus_names']:
            b[field] = sorted(b[field])
        b['previous_batch_count'] = len(b['previous_batches'])
        b['genus_count'] = len(b['genus_tax_ids'])
        b['group_count'] = len(b['group_ids'])
        b['unique_samples'] = len(b.pop('samples'))
        b['fastq_unique_samples'] = len(b.pop('fastq_samples'))
        b['batch_role'] = 'family_level_planning_container; samples stay separate'
        b['group_ids'].sort()
    result = sorted(batches.values(), key=lambda b: (b['family_name'].casefold(), b['family_batch_id']))
    save('family_catalog.jsonl', result)
    save('group_to_family.jsonl', sorted(group_rows, key=lambda r: (r['family_batch_id'], r['display_name'])))
    save('run_membership.jsonl', new_members)
    summary = dict(previous_batches=len(old_batches), family_batches=sum(bool(b['family_tax_id']) for b in result),
                   unresolved_retained_batches=sum(not bool(b['family_tax_id']) for b in result), total_batches=len(result),
                   runs=len(new_members), fastq_runs=sum(b['fastq_runs'] for b in result),
                   taxonomy_groups=len(groups),
                   multi_previous_batch_families=sum(bool(b['family_tax_id']) and b['previous_batch_count'] > 1 for b in result),
                   singleton_groups_brought_into_shared_batches=sum(b['singleton_groups'] for b in result if b['group_count'] > 1),
                   low_reference_runs=sum(b['low_reference_runs'] for b in result),
                   source_taxonomy_sha256=hashlib.sha256((SOURCE / 'taxon_mapping.jsonl').read_bytes()).hexdigest(),
                   policy='Group by family ancestor, not numeric proximity. No family assignment means retain prior group. Keep original run/group/genus membership and per-species review priorities.',
                   streptomycetaceae=next((b for b in result if b['family_name']=='Streptomycetaceae'), None))
    assert len(new_members) == len({m['run_accession'] for m in new_members}) == 93682
    assert sum(b['runs'] for b in result) == 93682
    assert summary['fastq_runs'] == 89134
    assert len(group_rows) == len({r['group_id'] for r in group_rows}) == len(groups)
    assert summary['low_reference_runs'] == 23564
    (OUT / 'summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    lines = ['# 按科合并的规划 batch', '',
             f"原 {len(old_batches):,} 个属级/未定属 batch → {summary['family_batches']:,} 个科级 batch，另保留 {summary['unresolved_retained_batches']:,} 个没有科级祖先的原分组，共 {len(result):,} 个规划 batch。", '',
             '合并依据是 NCBI taxonomy 中相同的 family 祖先 taxon ID。Taxon ID 是标识符，数值大小或差值不能作为亲缘距离。', '',
             '层级：科级规划 batch → 原属级 batch → 物种/编号类群 → BioSample → run。同一科里不同属可以共享一个规划 batch；原属级信息仍保留。', '',
             '这是批次组织方式，不是把样本、reads 或 assembly 混成一个样本，也不表示科内所有样本都适合直接比较。下载和执行任务可在选择范围后按样本量/体积拆成子批。', '',
             'E. coli 的低优先级只作用于该物种的成员，不传播到整个科。无明确 family 祖先的记录保持原 batch，不根据名称猜测。', '',
             '## 文件', '',
             '- [科级 batch 概览](FAMILY_BATCHES.md)',
             '- [科级清单](family_catalog.jsonl)：科、所含属、原 batch、分类组、去重样本数及容量。',
             '- [分类组到科级 batch 的映射](group_to_family.jsonl)',
             '- [每个 run 的归属](run_membership.jsonl)：保留原属级 batch_id，新增 family_batch_id。',
             '- [核对与示例](summary.json)', '',
             '计数基于 2026-10-24 冻结的 ENA 清单与已保存的 NCBI taxonomy 快照。没有下载新的 reads，全部 93,682 runs 均恰好保留一次。', '',
             '[NCBI Taxonomy 数据模型](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/data-processing/taxonomy-processing/taxonomy/)', '']
    (OUT / 'README.md').write_text('\n'.join(lines), encoding='utf-8')
    lines = ['# 科级规划 batch 一览', '', '按科名排序。无科祖先的记录列为未定；去重样本以 sample accession 为单位，不等于独立菌株数。', '',
             '| 科/待核查分组 | Family Taxon ID | 原 batch 数 | 属数 | 分类组数 | 单样本组数 | FASTQ 样本数 | Runs | FASTQ GB |',
             '|---|---|---:|---:|---:|---:|---:|---:|---:|']
    for b in result:
        values = [b['family_name'] or b['family_batch_id'], b['family_tax_id'] or '', b['previous_batch_count'],
                  b['genus_count'], b['group_count'], b['singleton_groups'], b['fastq_unique_samples'],
                  b['runs'], round(b['known_fastq_bytes']/1e9, 3)]
        lines.append('| ' + ' | '.join(esc(v) for v in values) + ' |')
    (OUT / 'FAMILY_BATCHES.md').write_text('\n'.join(lines) + '\n', encoding='utf-8')
    compact = dict(summary)
    if compact['streptomycetaceae']:
        compact['streptomycetaceae'] = {k:v for k,v in compact['streptomycetaceae'].items() if k != 'group_ids'}
    print(json.dumps(compact, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
