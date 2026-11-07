"""Annotate the frozen run inventory with submitter-reported organism/source metadata."""
import argparse
import collections
import hashlib
import json
import sqlite3
from datetime import datetime, timezone
from pathlib import Path

from build_ena_census import ROOT as BASE, QUERY, fetch

OUT = BASE / 'organism_environment'
ENV = ['isolation_source', 'environment_biome', 'environment_feature', 'environment_material',
       'broad_scale_environmental_context', 'local_environmental_context', 'environmental_medium']
HOST = ['host', 'host_scientific_name', 'host_tax_id', 'host_body_site']
FIELDS = ['run_accession', 'sample_accession', 'secondary_sample_accession', 'tax_id',
          'scientific_name', 'strain', 'sample_alias', 'sample_title', 'sample_description',
          *ENV, *HOST, 'country', 'location', 'collection_date']
MISSING = {'', 'na', 'n/a', 'n.a.', 'none', 'null', 'unknown', 'not known', 'not provided',
           'not available', 'not applicable', 'not collected', 'not recorded', 'unspecified',
           'missing', 'not determined', 'unavailable', '-', '?', 'nan',
           'not aplicable', 'not avaliable', 'not colleted', 'not collect',
           'not relevant', 'unknown source', 'unknown2', 'unknown4',
           'not provided submitted under migs 2.1'}


def meaningful(value):
    text = str(value or '').strip().lower()
    return text not in MISSING and not text.startswith(('missing:', 'not applicable:'))


def download():
    OUT.mkdir(exist_ok=True)
    started = datetime.now(timezone.utc).isoformat()
    url = fetch('search', dict(result='read_run', query=QUERY, format='json',
                               fields=','.join(FIELDS), limit=0), OUT / 'ena_context_raw.json')
    (OUT / 'provenance.json').write_text(json.dumps(dict(
        started_utc=started, completed_utc=datetime.now(timezone.utc).isoformat(),
        query=QUERY, url=url, fields=FIELDS, source='ENA public read_run indexed sample metadata',
        scope='Annotates frozen run list; newly indexed runs are reported separately, not silently added.'),
        indent=2) + '\n', encoding='utf-8')


def md_escape(value):
    return str(value).replace('|', '\\|').replace('\n', ' ').replace('\r', ' ')


def build():
    rows = json.loads((OUT / 'ena_context_raw.json').read_text(encoding='utf-8'))
    mapping = {r['run_accession']: r for r in rows}
    assert len(mapping) == len(rows), 'Duplicate context run IDs'
    base_conn = sqlite3.connect(BASE / 'run_inventory.sqlite')
    base_conn.row_factory = sqlite3.Row
    frozen = {r['run_accession']: dict(r) for r in base_conn.execute(
        'SELECT run_accession,sample_accession,tax_id,scientific_name,availability FROM runs')}
    base_conn.close()
    missing_ids = sorted(set(frozen) - set(mapping))
    extra_ids = sorted(set(mapping) - set(frozen))
    conflicts = []
    counts = collections.Counter()
    coverage = {f: collections.Counter() for f in FIELDS if f not in ['run_accession']}
    field_samples = {f: set() for f in coverage}
    organisms = {}
    value_agg = {}
    db = sqlite3.connect(OUT / 'organism_environment.sqlite')
    db.execute('DROP TABLE IF EXISTS run_context')
    derived_fields = ['availability', 'context_join_status', 'environment_metadata_status',
                      'environment_source_fields', 'host_metadata_present', 'geography_metadata_present',
                      'context_summary', 'meaningful_fields_json']
    db.execute('CREATE TABLE run_context (' + ','.join('"' + f + '" TEXT' for f in FIELDS + derived_fields) + ', PRIMARY KEY(run_accession))')
    with (OUT / 'run_organism_environment.jsonl').open('w', encoding='utf-8') as stream:
        for run, old in sorted(frozen.items()):
            raw = mapping.get(run)
            join_status = 'matched' if raw else 'missing_from_context_export'
            # A run/sample mismatch must not attach the context of a different sample.
            if raw and raw.get('sample_accession') != old['sample_accession']:
                conflicts.append(dict(run_accession=run, field='sample_accession', old=old['sample_accession'], new=raw.get('sample_accession')))
                join_status = 'sample_conflict_context_not_attached'
                raw = None
            if raw:
                for field in ['tax_id', 'scientific_name']:
                    if raw.get(field) != old[field]:
                        conflicts.append(dict(run_accession=run, field=field, old=old[field], new=raw.get(field)))
            record = {f: (raw.get(f, '') if raw else old.get(f, '')) for f in FIELDS}
            present = {f: meaningful(record[f]) for f in FIELDS}
            env_fields = [f for f in ENV if present[f]]
            host_fields = [f for f in HOST if present[f]]
            env_status = ('environment_or_isolation_recorded' if env_fields else
                          'host_context_only' if host_fields else 'no_structured_environment_or_host')
            geo = present['country'] or present['location']
            # These are literal evidence summaries, not inferred ecological labels.
            summary = '; '.join(f + '=' + record[f] for f in ENV + HOST if present[f])
            derived = dict(availability=old['availability'], context_join_status=join_status,
                           environment_metadata_status=env_status, environment_source_fields=';'.join(env_fields),
                           host_metadata_present=bool(host_fields), geography_metadata_present=geo,
                           context_summary=summary, meaningful_fields_json=json.dumps(present))
            record.update(derived)
            stream.write(json.dumps(record, ensure_ascii=False) + '\n')
            vals = [record[f] for f in FIELDS + derived_fields]
            db.execute('INSERT INTO run_context VALUES (' + ','.join('?' for _ in vals) + ')', vals)
            counts[env_status] += 1
            counts['host_context_recorded'] += bool(host_fields)
            counts['geography_recorded'] += geo
            sample = record['sample_accession'] or ('RUN:' + run)
            for f in coverage:
                if present[f]:
                    coverage[f]['meaningful_runs'] += 1
                    field_samples[f].add(sample)
                elif str(record[f] or '').strip():
                    coverage[f]['placeholder_runs'] += 1
                else:
                    coverage[f]['blank_runs'] += 1
            key = (record['tax_id'], record['scientific_name'])
            if key not in organisms:
                organisms[key] = dict(tax_id=key[0], scientific_name=key[1], runs=0, samples=set(),
                                      fastq_runs=0, source_metadata_runs=0, host_only_runs=0)
            item = organisms[key]
            item['runs'] += 1
            item['samples'].add(sample)
            item['fastq_runs'] += old['availability'] == 'fastq_listed'
            item['source_metadata_runs'] += bool(env_fields)
            item['host_only_runs'] += env_status == 'host_context_only'
            for f in ENV + HOST + ['country']:
                if not present[f]: continue
                val = record[f].strip()
                k = (f, val)
                if k not in value_agg:
                    value_agg[k] = dict(field=f, value=val, runs=0, samples=set())
                value_agg[k]['runs'] += 1
                value_agg[k]['samples'].add(sample)
    db.execute('CREATE INDEX context_taxon ON run_context(tax_id)')
    db.execute('CREATE INDEX context_sample ON run_context(sample_accession)')
    db.execute('CREATE INDEX context_availability ON run_context(availability)')
    db.execute('CREATE INDEX context_status ON run_context(environment_metadata_status)')
    db.commit()
    assert db.execute('PRAGMA integrity_check').fetchone()[0] == 'ok'
    assert db.execute('SELECT COUNT(*) FROM run_context').fetchone()[0] == len(frozen)
    db.close()
    orgs = sorted(organisms.values(), key=lambda r: (-r['runs'], r['tax_id'], r['scientific_name']))
    values = sorted(value_agg.values(), key=lambda r: (r['field'], -r['runs'], r['value']))
    for entries, filename in [(orgs, 'organisms.jsonl'), (values, 'reported_source_values.jsonl')]:
        with (OUT / filename).open('w', encoding='utf-8') as target:
            for entry in entries:
                entry['unique_samples'] = len(entry.pop('samples'))
                target.write(json.dumps(entry, ensure_ascii=False) + '\n')
    coverage_out = {f: dict(meaningful_runs=c['meaningful_runs'], placeholder_runs=c['placeholder_runs'],
                            blank_runs=c['blank_runs'], unique_samples=len(field_samples[f])) for f, c in coverage.items()}
    audit = dict(frozen_runs=len(frozen), context_export_runs=len(rows), missing_run_ids=missing_ids,
                 new_run_ids_not_added=extra_ids, metadata_changes=conflicts,
                 distinct_taxon_name_pairs=len(orgs), distinct_tax_ids=len({r['tax_id'] for r in orgs if r['tax_id']}),
                 counts=dict(counts), field_coverage=coverage_out,
                 raw_sha256=hashlib.sha256((OUT / 'ena_context_raw.json').read_bytes()).hexdigest(),
                 environment_inference='None: raw source and host fields are preserved. Taxonomic names and geography do not imply habitat.',
                 placeholder_policy=sorted(MISSING) + ['prefix missing:', 'prefix not applicable:'],
                 unique_sample_definition='sample_accession, or run-specific fallback if blank')
    (OUT / 'summary.json').write_text(json.dumps(audit, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')
    org_lines = ['# Organism 全量汇总', '',
                 '按 ENA 提交的 taxon ID / scientific name 汇总，按 run 数降序。不是经过统一 species-rank 归并的物种表；不代表样本选择优先级。', '',
                 '| Scientific name | Taxon ID | Runs | 去重样本数 | 有 FASTQ 的 runs | 有环境/分离来源的 runs |',
                 '|---|---|---:|---:|---:|---:|']
    for r in orgs:
        org_lines.append('| ' + ' | '.join(md_escape(r[k]) for k in ['scientific_name', 'tax_id', 'runs', 'unique_samples', 'fastq_runs', 'source_metadata_runs']) + ' |')
    (OUT / 'ORGANISMS.md').write_text('\n'.join(org_lines) + '\n', encoding='utf-8')
    source_lines = ['# 已报告的环境与来源原文', '',
                    '各字段展示 run 数最多的 30 个原始取值。大小写、拼写和同义词未合并；相同 run 可出现在多个字段中。完整值表见 reported_source_values.jsonl。', '',
                    '旧版与新版环境字段可能是别名，不能把它们的计数相加。宿主、地理来源与生态环境分开记录。', '']
    for field in ['isolation_source', 'environment_biome', 'environment_feature', 'environment_material', 'host_scientific_name', 'country']:
        source_lines += ['## ' + field, '', '| 原文 | Runs | 去重样本数 |', '|---|---:|---:|']
        for r in [v for v in values if v['field'] == field][:30]:
            source_lines.append(f"| {md_escape(r['value'])} | {r['runs']:,} | {r['unique_samples']:,} |")
        source_lines.append('')
    (OUT / 'SOURCES.md').write_text('\n'.join(source_lines), encoding='utf-8')
    lines = ['# Nanopore organism 与环境来源目录', '', '基于 2026-10-24 冻结的 ENA 细菌 ONT 候选池，逐 run 补充 ENA 索引中的样本元数据。', '',
             f'共 {len(frozen):,} runs，{len(orgs):,} 个 taxon ID / scientific name 组合。这里不是物种数：名称可能停留在属级，或细化到菌株级。', '',
             '## 环境信息覆盖', '', '| 信息类别 | Runs |', '|---|---:|']
    labels = [('environment_or_isolation_recorded', '至少一个分离来源或结构化环境字段有有效值'),
              ('host_context_only', '仅有宿主/宿主部位上下文'),
              ('no_structured_environment_or_host', '上述结构化字段均缺失或为占位值')]
    for k, label in labels: lines.append(f'| {label} | {counts[k]:,} |')
    lines += ['', '前三类互斥且覆盖全部 runs。字段缺失表示本次 ENA 索引未提供有效值，不代表原始 BioSample XML、论文或补充材料中也没有。', '',
              '## 逐字段覆盖', '', '| 字段 | 有效值 runs | 占位值 runs | 空白 runs |', '|---|---:|---:|---:|']
    for f in ['scientific_name', 'strain', *ENV, *HOST, 'country', 'location', 'collection_date']:
        c = coverage_out[f]
        lines.append(f"| {f} | {c['meaningful_runs']:,} | {c['placeholder_runs']:,} | {c['blank_runs']:,} |")
    lines += ['', '## 如何查看每条记录', '',
              '- [逐 run 完整注释](run_organism_environment.jsonl)：每行一个 run，保留 organism、taxon、strain、BioSample、环境原文、宿主、地点、采集日期和下载可用状态。',
              '- [本地检索数据库](organism_environment.sqlite)：`run_context` 表，可按 organism、环境字段、是否有 FASTQ 等组合检索。',
              '- [organism 汇总](organisms.jsonl)：每个 taxon/name 的 run 数、去重样本数及来源信息覆盖。',
              '- [可直接阅读的 organism 全量表](ORGANISMS.md)；[环境与来源原文概览](SOURCES.md)。',
              '- [来源原文汇总](reported_source_values.jsonl)：按字段和值汇总 run 数与去重样本数，未把相似词自动合并。',
              '- [完整统计和对齐核对](summary.json)；[来源请求](provenance.json)。', '',
              '## 解读原则', '',
              'organism 是提交者标注的测序对象；host 是样本来源宿主，两者不混用。country/location 是地理来源，不直接等同于生态环境。', '',
              '保留 environmental biome、feature、material 以及新版 broad/local context、medium 字段；不从物种名称猜测其环境。原始字段保留 not provided、missing 等占位字符串，但统计时将这些值视为缺失。', '',
              '`context_summary` 只拼接有值字段及其原文，不做自动生态分类。sample title/description 保留供核查，未从自由文本猜测环境。', '',
              '同一 sample 可以有多个 runs；来源汇总同时给出去重样本数。不同字段的统计可以重叠，不可相加为总样本数。', '',
              '本轮没有查询原始样本 XML、做物种鉴定、判定宿主健康状态或筛选适合的样本；没有下载 reads。', '',
              '## 索引对齐核对', '',
              f'- 冻结列表缺失的 context runs：{len(missing_ids):,}。',
              f'- 新增而未纳入冻结底表的 runs：{len(extra_ids):,}。',
              f'- sample / taxon / scientific name 变更记录：{len(conflicts):,}。', '',
              '[ENA 样本元数据返回字段定义](https://www.ebi.ac.uk/ena/portal/api/returnFields?result=read_run&format=json)', '']
    (OUT / 'README.md').write_text('\n'.join(lines), encoding='utf-8')
    print(json.dumps({k: v for k, v in audit.items() if k not in ['placeholder_policy', 'raw_sha256']}, ensure_ascii=False, indent=2), flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--fetch', action='store_true')
    args = p.parse_args()
    if args.fetch: download()
    build()
