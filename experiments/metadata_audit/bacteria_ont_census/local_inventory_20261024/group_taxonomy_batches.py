"""Administrative taxonomy grouping; no sequence analysis or novelty ranking."""
import collections
import hashlib
import json
import sqlite3
import tarfile
from array import array
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent
BASE = ROOT / 'ena_bacteria_ont_2026-10-24'
OUT = BASE / 'taxonomy_batches'
ARCHIVE = ROOT / 'taxonomy_cache' / 'taxdump.tar.gz'
RANKS = {'species': 1, 'genus': 2, 'family': 3}
REVIEW_BANDS = {
    'very_dense_1000_plus': 'B01_very_dense',
    'dense_100_999': 'B02_dense',
    'medium_10_99': 'B03_medium',
    'sparse_2_9': 'B04_sparse',
    'singleton_1': 'B05_singleton',
    'no_fastq_sample': 'B06_no_fastq',
}


def read_taxonomy(ids):
    parents = array('I', [0]) * 4000000
    ranks = bytearray(4000000)
    merged = {}
    with tarfile.open(ARCHIVE, 'r:gz') as tf:
        deleted = {int(raw.split(b'|')[0]) for raw in tf.extractfile('delnodes.dmp')}
        for raw in tf.extractfile('merged.dmp'):
            parts = raw.split(b'|')
            merged[int(parts[0])] = int(parts[1])
        for raw in tf.extractfile('nodes.dmp'):
            parts = raw.split(b'|')
            taxid, parent = int(parts[0]), int(parts[1])
            if taxid >= len(parents):
                grow = taxid + 100000 - len(parents)
                parents.extend(array('I', [0]) * grow)
                ranks.extend(bytes(grow))
            parents[taxid] = parent
            ranks[taxid] = RANKS.get(parts[2].strip().decode(), 0)
    records = {}
    needed = set()
    for old in sorted(ids):
        current = old
        seen = set()
        while current in merged:
            if current in seen: raise ValueError('Cycle in merged taxids')
            seen.add(current)
            current = merged[current]
        record = dict(original_tax_id=old, current_tax_id=current,
                      taxid_was_merged=current != old,
                      taxid_listed_as_deleted=current in deleted,
                      species_tax_id=None, genus_tax_id=None, family_tax_id=None)
        node = current
        seen = set()
        while node and node < len(parents) and parents[node] and node not in seen:
            seen.add(node)
            rank = ranks[node]
            for key, code in RANKS.items():
                if rank == code and record[key + '_tax_id'] is None:
                    record[key + '_tax_id'] = node
            if node == parents[node]: break
            node = parents[node]
        record['taxonomy_resolved'] = current < len(parents) and bool(parents[current])
        record['grouping_status'] = 'species_resolved' if record['species_tax_id'] else 'no_species_ancestor' if record['taxonomy_resolved'] else 'taxid_unresolved'
        needed.add(current)
        needed.update(record[k + '_tax_id'] for k in RANKS if record[k + '_tax_id'])
        records[old] = record
    names = {}
    with tarfile.open(ARCHIVE, 'r:gz') as tf:
        for raw in tf.extractfile('names.dmp'):
            parts = raw.split(b'|')
            taxid = int(parts[0])
            if taxid in needed and parts[3].strip() == b'scientific name':
                names[taxid] = parts[1].strip().decode('utf-8')
    for r in records.values():
        r['current_scientific_name'] = names.get(r['current_tax_id'], '')
        for rank in RANKS:
            r[rank + '_name'] = names.get(r[rank + '_tax_id'], '')
    return records


def density(n):
    if n == 0: return 'no_fastq_sample'
    if n == 1: return 'singleton_1'
    if n < 10: return 'sparse_2_9'
    if n < 100: return 'medium_10_99'
    if n < 1000: return 'dense_100_999'
    return 'very_dense_1000_plus'


def write_jsonl(name, rows):
    with (OUT / name).open('w', encoding='utf-8') as target:
        for r in rows:
            target.write(json.dumps(r, ensure_ascii=False) + '\n')


def esc(value):
    return str(value).replace('|', '\\|').replace('\n', ' ').replace('\r', ' ')


def main():
    OUT.mkdir(exist_ok=True)
    con = sqlite3.connect(BASE / 'run_inventory.sqlite')
    con.row_factory = sqlite3.Row
    runs = [dict(r) for r in con.execute('SELECT run_accession,tax_id,scientific_name,sample_accession,study_accession,availability,library_strategy,fastq_total_bytes FROM runs')]
    con.close()
    contexts = sqlite3.connect(BASE / 'organism_environment' / 'organism_environment.sqlite')
    contexts.row_factory = sqlite3.Row
    context = {r['run_accession']: dict(r) for r in contexts.execute('SELECT run_accession,environment_metadata_status FROM run_context')}
    contexts.close()
    ids = {int(r['tax_id']) for r in runs}
    print('Resolving', len(ids), 'taxids using official NCBI ranks', flush=True)
    mapping = read_taxonomy(ids)
    write_jsonl('taxon_mapping.jsonl', mapping.values())
    species = {}
    memberships = []
    for r in runs:
        m = mapping[int(r['tax_id'])]
        # Genus-only/unclassified records retain their own taxon; never merge all sp. labels together.
        group_id = 'species_' + str(m['species_tax_id']) if m['species_tax_id'] else 'unresolved_taxon_' + str(m['current_tax_id'])
        batch_id = 'genus_' + str(m['genus_tax_id']) if m['genus_tax_id'] else 'unplaced_taxon_' + str(m['current_tax_id'])
        member = dict(run_accession=r['run_accession'], sample_accession=r['sample_accession'],
                      original_tax_id=int(r['tax_id']), original_scientific_name=r['scientific_name'],
                      group_id=group_id, batch_id=batch_id,
                      grouping_status=m['grouping_status'], availability=r['availability'])
        memberships.append(member)
        if group_id not in species:
            species[group_id] = dict(group_id=group_id, batch_id=batch_id,
                species_tax_id=m['species_tax_id'], species_name=m['species_name'],
                display_name=m['species_name'] or m['current_scientific_name'] or r['scientific_name'],
                genus_tax_id=m['genus_tax_id'], genus_name=m['genus_name'], family_name=m['family_name'],
                grouping_status=m['grouping_status'], runs=0, fastq_runs=0, known_fastq_bytes=0,
                samples=set(), fastq_samples=set(), wgs_fastq_samples=set(), studies=set(),
                original_tax_ids=set(), source_samples=set(), host_only_samples=set())
        group = species[group_id]
        group['runs'] += 1
        sample = r['sample_accession'] or 'RUN:' + r['run_accession']
        group['samples'].add(sample)
        group['studies'].add(r['study_accession'])
        group['original_tax_ids'].add(int(r['tax_id']))
        if r['availability'] == 'fastq_listed':
            group['fastq_runs'] += 1
            group['fastq_samples'].add(sample)
            group['known_fastq_bytes'] += int(r['fastq_total_bytes'] or 0)
            if r['library_strategy'] == 'WGS': group['wgs_fastq_samples'].add(sample)
        status = context[r['run_accession']]['environment_metadata_status']
        if status == 'environment_or_isolation_recorded': group['source_samples'].add(sample)
        if status == 'host_context_only': group['host_only_samples'].add(sample)
    batches = {}
    for group in species.values():
        batch_id = group['batch_id']
        if batch_id not in batches:
            batches[batch_id] = dict(batch_id=batch_id, genus_tax_id=group['genus_tax_id'],
                                    genus_name=group['genus_name'], species_group_count=0,
                                    unresolved_group_count=0, runs=0, samples=set(), fastq_samples=set(),
                                    known_fastq_bytes=0, group_ids=[])
        b = batches[batch_id]
        b['species_group_count'] += bool(group['species_tax_id'])
        b['unresolved_group_count'] += not bool(group['species_tax_id'])
        b['runs'] += group['runs']
        b['samples'].update(group['samples'])
        b['fastq_samples'].update(group['fastq_samples'])
        b['known_fastq_bytes'] += group['known_fastq_bytes']
        b['group_ids'].append(group['group_id'])
    for group in species.values():
        group['unique_samples'] = len(group.pop('samples'))
        group['fastq_unique_samples'] = len(group.pop('fastq_samples'))
        group['wgs_fastq_unique_samples'] = len(group.pop('wgs_fastq_samples'))
        group['distinct_studies'] = len(group.pop('studies') - {''})
        group['source_unique_samples'] = len(group.pop('source_samples'))
        group['host_only_unique_samples'] = len(group.pop('host_only_samples'))
        group['original_tax_ids'] = sorted(group['original_tax_ids'])
        group['data_density_band'] = density(group['fastq_unique_samples'])
        group['review_batch_id'] = REVIEW_BANDS[group['data_density_band']] if group['species_tax_id'] else 'B07_unresolved_species'
        group['literature_saturation_status'] = 'not_assessed'
        group['assembly_link_status'] = 'not_assessed'
        group['user_example_ecoli'] = group['species_tax_id'] == 562
        group['review_priority'] = 'low_reference' if group['user_example_ecoli'] else 'not_assessed'
        group['priority_evidence'] = 'User confirmed on 2026-10-24: lower priority for E. coli; retain as reference' if group['user_example_ecoli'] else ''
        group['selection_status'] = 'not_selected'
    for b in batches.values():
        b['unique_samples'] = len(b.pop('samples'))
        b['fastq_unique_samples'] = len(b.pop('fastq_samples'))
        b['group_ids'].sort()
        b['batch_type'] = 'logical_genus_inventory_not_download_job'
        b['single_fastq_sample_groups'] = sum(species[gid]['fastq_unique_samples'] == 1 for gid in b['group_ids'])
    groups = sorted(species.values(), key=lambda g: (g['genus_name'].casefold(), g['display_name'].casefold(), g['group_id']))
    batch_list = sorted(batches.values(), key=lambda b: (b['genus_name'].casefold(), b['batch_id']))
    write_jsonl('species_groups.jsonl', groups)
    write_jsonl('batch_catalog.jsonl', batch_list)
    write_jsonl('run_membership.jsonl', sorted(memberships, key=lambda m: m['run_accession']))
    review_batches = []
    for batch_id in list(REVIEW_BANDS.values()) + ['B07_unresolved_species']:
        selected = [g for g in groups if g['review_batch_id'] == batch_id]
        review_batches.append(dict(review_batch_id=batch_id, group_count=len(selected),
                                   runs=sum(g['runs'] for g in selected),
                                   fastq_runs=sum(g['fastq_runs'] for g in selected),
                                   group_ids=[g['group_id'] for g in selected],
                                   purpose='Metadata review stratum, not novelty rank or download job'))
    (OUT / 'review_batches.json').write_text(json.dumps(review_batches, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    bands = {}
    for band in ['very_dense_1000_plus', 'dense_100_999', 'medium_10_99', 'sparse_2_9', 'singleton_1', 'no_fastq_sample']:
        selected = [g for g in groups if g['species_tax_id'] and g['data_density_band'] == band]
        bands[band] = dict(species_groups=len(selected), runs=sum(g['runs'] for g in selected),
                           fastq_samples_sum_per_species=sum(g['fastq_unique_samples'] for g in selected))
    stats = dict(run_count=len(runs), original_tax_ids=len(ids),
                 species_groups=sum(bool(g['species_tax_id']) for g in groups),
                 unresolved_groups=sum(not bool(g['species_tax_id']) for g in groups),
                 logical_batches=len(batches), genus_batches=sum(bool(b['genus_tax_id']) for b in batches.values()),
                 unplaced_batches=sum(not bool(b['genus_tax_id']) for b in batches.values()),
                 merged_taxids=sum(m['taxid_was_merged'] for m in mapping.values()),
                 unresolved_taxids=sum(not m['taxonomy_resolved'] for m in mapping.values()),
                 deleted_taxids=sum(m['taxid_listed_as_deleted'] for m in mapping.values()),
                 density_bands=bands,
                 ecoli=next((g for g in groups if g['species_tax_id'] == 562), None),
                 generated_utc=datetime.now(timezone.utc).isoformat(),
                 taxonomy_archive_sha256=hashlib.sha256(ARCHIVE.read_bytes()).hexdigest(),
                 taxonomy_source='https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz',
                 evidence='NCBI taxonomy ranks and ENA inventory counts; no literature saturation measurement')
    assert sum(g['runs'] for g in groups) == len(runs) == 93682
    assert sum(b['runs'] for b in batches.values()) == len(runs)
    assert sum(g['fastq_runs'] for g in groups) == 89134
    assert len({m['run_accession'] for m in memberships}) == len(runs)
    (OUT / 'summary.json').write_text(json.dumps(stats, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    lines = ['# 物种归并与逻辑 batch', '',
             '2026-10-24 ENA 细菌 ONT 候选池；用 NCBI taxonomy 的 species / genus 层级归并。', '',
             f"原 {len(ids):,} 个 taxon ID → {stats['species_groups']:,} 个 species 层级分组，另有 {stats['unresolved_groups']:,} 个无法归到 species 层级的保留组。", '',
             f"生成 {stats['genus_batches']:,} 个 genus 逻辑 batch，另有 {stats['unplaced_batches']:,} 个无 genus 的单独 batch。全部 {len(runs):,} 个 runs 均保留且只归属一个组。", '',
             '## 实际分批的组织方式', '',
             '以 genus batch 为组织单元，多个 species / sp. 编号共享一个 batch。只有一个 run 或样本的 species group 不需要单独成批。物种归并表中的每行是分类单位，不是一个处理任务。', '',
             '例如 Streptomyces sp. NBC_01197 归入 genus_1883（Streptomyces），与其他 Streptomyces 类群一起规划处理。组内每个样本的原始 run、reads 和未来 assembly 对应关系仍独立保留；同 batch 不代表拼接 reads 或合成一个生物样本。', '',
             '数据密度仅作筛选标签，不是分批边界；同一 genus batch 可以同时包含 singleton、稀疏和高密度物种。较大的 genus batch 在决定下载范围后再按体积或样本量拆为执行子批。', '',
             '## 数据密度标签（不是处理 batch）', '',
             '以下标签按有 FASTQ 的去重 sample accession 数划分，阈值只是清单统计规则，不代表研究热度或新颖性。旧文件中的 review_batch_id / review_batches.json 是这些审阅标签的历史名称，处理归属应读取 batch_id / batch_catalog.jsonl。', '',
             '| 审阅标签（旧 ID） | 每个 species 的 FASTQ 样本数 | Species 分组数 |', '|---|---|---:|']
    for band, val in bands.items(): lines.append(f"| {REVIEW_BANDS[band]} | {band} | {val['species_groups']:,} |")
    lines += [f"| B07_unresolved_species | 物种层级待核查 | {stats['unresolved_groups']:,} |"]
    lines += ['', '## 文件入口', '',
              '- [可阅读的物种归并表](SPECIES_GROUPS.md)：名称、属、样本数、研究项目数、数据密度及原 taxon 数。',
              '- [可阅读的 batch 概览](BATCHES.md)：每个属的物种组数、样本数和登记体积；点击 Batch 可查看成员。',
              '- [全部 batch 的成员明细](BATCH_MEMBERS.md)：包含只有一个样本的 species group，统一按属组织。',
              '- [七个数据密度审阅标签（旧文件名）](review_batches.json)：这些统计分层不决定处理 batch，也不是新颖性排名。',
              '- [物种分组数据](species_groups.jsonl)；[batch 清单](batch_catalog.jsonl)；[逐 run 归属](run_membership.jsonl)；[分类映射证据](taxon_mapping.jsonl)。',
              '- [汇总与快照来源](summary.json)。', '',
              '## 怎样理解“没被挖透”', '',
              '本轮只完成分类归并、数据密度和来源覆盖统计。数据少可能是研究少，也可能是上传少、缺乏组装、分类名称变化或测序平台偏好，不能直接当作新颖性。', '',
              '后续应结合具体研究主题的文献、既有系统筛查、样本独立性、环境覆盖及同一 BioSample 的 assembly 可用性，再判断是否值得优先研究。没有为全表计算“研究饱和度”，没有据此选定样本。', '',
              '用户于 2026-10-24 确认：降低 E. coli 等研究充分类群的优先级。E. coli（species taxon 562）已标为 low_reference（低优先级参考组），保留全部记录。其他类群的 review_priority 暂为 not_assessed，需结合文献逐项判断，不按数据量自动降级。不把整个 Escherichia 属或 Shigella 自动并入。', '',
              '同一物种下菌株级 taxon 会合并；无 species 祖先的记录保留各自 taxon，不把所有“sp.”或“uncultured”误合成同一物种。', '',
              f"本次有 {stats['deleted_taxids']} 个原 taxon ID 在当前 NCBI 快照中已删除；这些记录保留原名称并标为待核查，没有凭名称强行匹配。", '',
              'species-rank 中仍可能有未正式命名的 sp. 类群；species 层级并不保证正式双名或鉴定可靠。', '',
              'Studies 数为不同 study accession 的数量，不意味着独立实验或独立样本。相同 sample 跨物种出现时，各物种内去重，但跨组统计不可直接当成全局唯一菌株数。', '',
              'batch 是用于浏览和后续规划的属级集合，不是可以立即提交的下载/计算任务。超大属后续可以按所选样本量或体积拆分。原始文件和来源字段未修改。', '',
              '[NCBI taxonomy 层级与名称说明](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/data-processing/taxonomy-processing/taxonomy/)', '']
    (OUT / 'README.md').write_text('\n'.join(lines), encoding='utf-8')
    lines = ['# Species 分组表', '', '每行是分类分组，不是一个独立 batch。同一 Batch ID 的物种可以一起规划处理；样本身份与文件不合并。数据密度不是研究饱和度；样本数按 sample accession 去重。', '',
             '| 名称 | Genus | 分组 ID | Batch ID | FASTQ 样本 | WGS FASTQ 样本 | 全部 runs | Studies | 原 taxon 数 | 数据密度 | 审阅优先级 |',
             '|---|---|---|---|---:|---:|---:|---:|---:|---|---|']
    for g in groups:
        vals = [g['display_name'], g['genus_name'], g['group_id'], g['batch_id'], g['fastq_unique_samples'], g['wgs_fastq_unique_samples'], g['runs'], g['distinct_studies'], len(g['original_tax_ids']), g['data_density_band'], g['review_priority']]
        lines.append('| ' + ' | '.join(esc(v) for v in vals) + ' |')
    (OUT / 'SPECIES_GROUPS.md').write_text('\n'.join(lines) + '\n', encoding='utf-8')
    lines = ['# 属级逻辑 batch', '', '按名称排序；不是优先级排序，不是下载任务。没有 genus 的记录单列。', '',
             '| Batch | Genus | Species 组 | 未定 species 组 | 单 FASTQ 样本组 | FASTQ 样本 | Runs | FASTQ GB |',
             '|---|---|---:|---:|---:|---:|---:|---:|']
    for b in batch_list:
        vals = [f"[{b['batch_id']}](BATCH_MEMBERS.md#{b['batch_id']})", b['genus_name'], b['species_group_count'], b['unresolved_group_count'], b['single_fastq_sample_groups'], b['fastq_unique_samples'], b['runs'], round(b['known_fastq_bytes']/1e9, 3)]
        lines.append('| ' + ' | '.join(esc(v) for v in vals) + ' |')
    (OUT / 'BATCHES.md').write_text('\n'.join(lines) + '\n', encoding='utf-8')
    lines = ['# 按属组织的 batch 成员', '', '每个 batch 可包含多个物种/编号类群；样本和原始文件保持独立。单样本类群不单独成批。', '']
    for b in batch_list:
        lines += ['## ' + b['batch_id'], '',
                  f"{b['genus_name'] or '属未定'}：{len(b['group_ids']):,} 个分类组，{b['fastq_unique_samples']:,} 个有 FASTQ 的去重样本，{b['runs']:,} 个 runs。", '',
                  '| 名称 | 分类组 | FASTQ 样本数 | Runs | 审阅优先级 |', '|---|---|---:|---:|---|']
        for gid in sorted(b['group_ids'], key=lambda k: species[k]['display_name'].casefold()):
            g = species[gid]
            lines.append('| ' + ' | '.join(esc(g[k]) for k in ['display_name', 'group_id', 'fastq_unique_samples', 'runs', 'review_priority']) + ' |')
        lines.append('')
    (OUT / 'BATCH_MEMBERS.md').write_text('\n'.join(lines), encoding='utf-8')
    print(json.dumps(stats, ensure_ascii=False, indent=2), flush=True)


if __name__ == '__main__':
    main()
