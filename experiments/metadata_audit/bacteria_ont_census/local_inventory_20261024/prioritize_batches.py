"""Rank metadata-review batches using sample sufficiency and generic literature attention."""
import argparse
import collections
import concurrent.futures
import hashlib
import json
import math
import re
import sqlite3
import subprocess
import shutil
import tarfile
import time
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlencode, quote

BASE = Path(__file__).resolve().parent / 'ena_bacteria_ont_2026-10-24'
OUT = BASE / 'batch_priority'
CACHE = OUT / 'literature_cache'
MIN_SAMPLES = 100  # Explicitly confirmed by user, 2026-10-24.
MIN_STUDIES = 2
# Curated nomenclatural aliases absent from the downloaded NCBI names snapshot.
VERIFIED_ALIASES = {
    2100: ('Mycoplasma hyorhinis', 'https://lpsn.dsmz.de/text/taxonomic-suspension'),
    29562: ('Mycoplasma ovipneumoniae', 'https://lpsn.dsmz.de/text/taxonomic-suspension'),
    92402: ('Mycoplasma equirhinis', 'https://lpsn.dsmz.de/species/metamycoplasma-equirhinis'),
    56142: ('Mycoplasma faucium', 'https://lpsn.dsmz.de/species/metamycoplasma-faucium'),
}


def lines(path):
    return [json.loads(x) for x in path.read_text(encoding='utf-8').splitlines()]


def search(query):
    key = hashlib.sha256(query.encode()).hexdigest()
    path = CACHE / (key + '.json')
    if path.exists():
        record = json.loads(path.read_text(encoding='utf-8'))
        if record.get('status') == 'ok' and record['query'] == query:
            return record
    url = 'https://www.ebi.ac.uk/europepmc/webservices/rest/search?' + urlencode(
        dict(query=query, format='json', resultType='idlist', pageSize=1))
    error = ''
    for attempt in range(3):
        try:
            proc = subprocess.run([shutil.which('curl.exe') or shutil.which('curl') or 'curl', '--silent', '--show-error', '--fail',
                                   '--connect-timeout', '12', '--max-time', '40', url],
                                  capture_output=True, text=True, encoding='utf-8', timeout=45)
            if proc.returncode: raise RuntimeError(proc.stderr[:300])
            raw = json.loads(proc.stdout)
            if not isinstance(raw.get('hitCount'), int): raise ValueError('API omitted integer hitCount')
            if raw.get('request', {}).get('queryString') != query: raise ValueError('Query echo mismatch')
            record = dict(query=query, url=url, fetched_utc=datetime.now(timezone.utc).isoformat(),
                          status='ok', hit_count=raw['hitCount'], response=raw)
            path.write_text(json.dumps(record, ensure_ascii=False) + '\n', encoding='utf-8')
            return record
        except Exception as exc:
            error = str(exc)
            if attempt < 2: time.sleep(1 + attempt)
    record = dict(query=query, url=url, fetched_utc=datetime.now(timezone.utc).isoformat(), status='error', error=error)
    path.write_text(json.dumps(record, ensure_ascii=False) + '\n', encoding='utf-8')
    return record


def query_for(names):
    terms = ['TITLE_ABS:"' + name.replace('"', '') + '"' for name in sorted(set(names)) if name]
    if not terms: return None
    return 'SRC:MED AND (' + ' OR '.join(terms) + ')'


def original_genus(name):
    parts = name.replace('[', '').replace(']', '').split()
    if not parts: return None
    if parts[0].lower() in ['uncultured', 'unclassified', 'bacterium', 'bacteria', 'metagenome']: return None
    if parts[0] == 'Candidatus': return ' '.join(parts[:2]) if len(parts) > 1 else None
    return parts[0] if parts[0][0].isupper() and parts[0].isalpha() else None


def build(fetch_literature=False):
    OUT.mkdir(exist_ok=True)
    CACHE.mkdir(exist_ok=True)
    families = {b['family_batch_id']: b for b in lines(BASE / 'family_batches/family_catalog.jsonl')}
    members = lines(BASE / 'family_batches/run_membership.jsonl')
    groups = {g['group_id']:g for g in lines(BASE / 'taxonomy_batches/species_groups.jsonl')}
    maps = {m['original_tax_id']:m for m in lines(BASE / 'taxonomy_batches/taxon_mapping.jsonl')}
    conn = sqlite3.connect(BASE / 'run_inventory.sqlite')
    conn.row_factory = sqlite3.Row
    runs = {r['run_accession']:dict(r) for r in conn.execute('SELECT run_accession,library_strategy,library_source,study_accession,fastq_total_bytes,base_count,scientific_name FROM runs')}
    conn.close()
    data = {bid:dict(all_wgs_samples=set(), samples=set(), studies=set(), sample_studies=collections.defaultdict(set),
                    sample_genus=collections.defaultdict(set), sample_bases=collections.Counter(), genomic_samples=set(),
                    wgs_runs=0, bytes=0, group_ids=set(), genera={}) for bid in families}
    for m in members:
        r = runs[m['run_accession']]
        if m['availability'] != 'fastq_listed' or r['library_strategy'] != 'WGS': continue
        bid = m['family_batch_id']
        d = data[bid]
        sample = m['sample_accession'] or ('RUN:' + m['run_accession'])
        d['all_wgs_samples'].add(sample)
        if m['review_priority'] == 'low_reference': continue
        d['samples'].add(sample)
        if r['study_accession']:
            d['studies'].add(r['study_accession'])
            d['sample_studies'][sample].add(r['study_accession'])
        d['wgs_runs'] += 1
        d['bytes'] += int(r['fastq_total_bytes'] or 0)
        d['group_ids'].add(m['group_id'])
        if r['base_count'].isdigit(): d['sample_bases'][sample] += int(r['base_count'])
        if r['library_source'] == 'GENOMIC': d['genomic_samples'].add(sample)
        g = groups[m['group_id']]
        gid = g['genus_tax_id']
        unit = 'genus_' + str(gid) if gid else 'family_fallback_' + str(families[bid]['family_tax_id'])
        d['sample_genus'][sample].add(unit)
        if unit not in d['genera']:
            d['genera'][unit] = dict(unit_id=unit, genus_tax_id=gid,
                                   label=g['genus_name'] or families[bid]['family_name'], names=set(), original_names=set())
        u = d['genera'][unit]
        if u['label']: u['names'].add(u['label'])
        name = original_genus(r['scientific_name'])
        if name and gid:
            u['original_names'].add(name)
    qualified = [bid for bid,d in data.items() if len(d['samples']) >= MIN_SAMPLES and len(d['studies']) >= MIN_STUDIES and families[bid]['family_tax_id']]
    units = {}
    for bid in qualified:
        for unit,u in data[bid]['genera'].items():
            if unit not in units: units[unit] = dict(u, names=set(u['names']), original_names=set(u['original_names']))
            else:
                units[unit]['names'].update(u['names'])
                units[unit]['original_names'].update(u['original_names'])
    # Official scientific/synonym/equivalent genus names mitigate simple renaming bias.
    wanted = {u['genus_tax_id']:u for u in units.values() if u['genus_tax_id']}
    species_units = {}
    for bid in qualified:
        for gid in data[bid]['group_ids']:
            g = groups[gid]
            unit = 'genus_' + str(g['genus_tax_id'])
            if g['species_tax_id'] and unit in units:
                species_units[g['species_tax_id']] = units[unit]
    for u in units.values(): u['species_alias_terms'] = set()
    archive = BASE.parent / 'taxonomy_cache/taxdump.tar.gz'
    with tarfile.open(archive, 'r:gz') as tf:
        for raw in tf.extractfile('names.dmp'):
            p = raw.split(b'|')
            taxid = int(p[0])
            if taxid in wanted and p[3].strip() in [b'scientific name', b'synonym', b'equivalent name']:
                name = p[1].strip().decode('utf-8')
                if name and len(name) >= 4: wanted[taxid]['names'].add(name)
            if taxid in species_units and p[3].strip() in [b'synonym', b'equivalent name', b'scientific name']:
                name = p[1].strip().decode('utf-8')
                prior_genus = original_genus(name)
                u = species_units[taxid]
                if prior_genus and prior_genus != u['label'] and ' ' in name:
                    u['species_alias_terms'].add(name)
                    u['names'].add(name)
                    # NCBI occasionally includes naming authorities inside synonym strings.
                    canonical = re.match(r'^(?:Candidatus )?[A-Z][a-z]+ [a-z][a-z-]+', name.strip('"'))
                    if canonical and canonical.group().split()[-1] not in ['sp', 'spp', 'serotype', 'bacterium']:
                        u['species_alias_terms'].add(canonical.group())
                        u['names'].add(canonical.group())
    for taxid, (alias, source) in VERIFIED_ALIASES.items():
        if taxid in species_units:
            species_units[taxid]['names'].add(alias)
            species_units[taxid]['species_alias_terms'].add(alias)
    (OUT / 'curated_aliases.json').write_text(json.dumps({str(k):dict(name=v[0],source=v[1]) for k,v in VERIFIED_ALIASES.items()}, indent=2) + '\n', encoding='utf-8')
    for u in units.values():
        u['names'] = sorted(u['names'])
        u['original_names'] = sorted(u['original_names'])
        u['species_alias_terms'] = sorted(u['species_alias_terms'])
        u['query'] = query_for(u['names'])
    (OUT / 'query_plan.json').write_text(json.dumps(dict(qualified_families=len(qualified), units=list(units.values())), ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    print('Data-qualified families:', len(qualified), '; literature queries:', len(units), flush=True)
    results = {}
    if fetch_literature:
        with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
            jobs = {pool.submit(search, u['query']): unit for unit,u in units.items() if u['query']}
            for index, future in enumerate(concurrent.futures.as_completed(jobs), 1):
                results[jobs[future]] = future.result()
                if index % 25 == 0 or index == len(jobs): print('Literature queries complete:', index, '/', len(jobs), flush=True)
    else:
        for unit,u in units.items():
            if not u['query']: continue
            path = CACHE / (hashlib.sha256(u['query'].encode()).hexdigest() + '.json')
            if path.exists(): results[unit] = json.loads(path.read_text(encoding='utf-8'))
    audits = []
    rows = []
    for bid,b in families.items():
        d = data[bid]
        count = len(d['samples'])
        studies = len(d['studies'])
        study_counts = collections.Counter(study for studyset in d['sample_studies'].values() for study in studyset)
        top_project_share = max(study_counts.values(), default=0) / count if count else None
        row = dict(family_batch_id=bid, family_name=b['family_name'], family_tax_id=b['family_tax_id'],
                   all_wgs_fastq_samples=len(d['all_wgs_samples']), eligible_wgs_fastq_samples=count,
                   eligible_genomic_wgs_samples=len(d['genomic_samples']), distinct_studies=studies,
                   eligible_wgs_runs=d['wgs_runs'], eligible_fastq_bytes=d['bytes'],
                   eligible_species_groups=len(d['group_ids']), largest_study_sample_share=top_project_share,
                   data_gate_passed=count >= MIN_SAMPLES and studies >= MIN_STUDIES,
                   attention_score=None, literature_coverage=0, attention_review='not_assessed',
                   assembly_pairing_status='not_assessed', assembly_verified_sample_count=None,
                   min_samples=MIN_SAMPLES, min_studies=MIN_STUDIES)
        row['review_flags'] = []
        if count >= MIN_SAMPLES and len(d['genomic_samples']) < MIN_SAMPLES:
            row['review_flags'].append('fewer_than_100_GENOMIC_WGS_samples')
        if top_project_share is not None and top_project_share > .8:
            row['review_flags'].append('over_80pct_samples_in_one_study')
        weights = collections.Counter()
        for sample, unit_set in d['sample_genus'].items():
            for unit in unit_set: weights[unit] += 1 / len(unit_set)
        known_weight = 0
        logsum = 0
        if bid in qualified:
            for unit, weight in weights.items():
                res = results.get(unit, {})
                ok = res.get('status') == 'ok'
                if ok:
                    known_weight += weight
                    logsum += weight * math.log10(1 + res['hit_count'])
                audits.append(dict(family_batch_id=bid, unit_id=unit, name=units[unit]['label'],
                                   terms=units[unit]['names'], sample_weight=weight,
                                   pubmed_title_abstract_hits=res.get('hit_count'),
                                   query=units[unit]['query'], url=res.get('url'), status=res.get('status', 'missing'),
                                   fetched_utc=res.get('fetched_utc')))
            row['literature_coverage'] = known_weight / count if count else 0
            if abs(known_weight - count) < 1e-6:
                row['attention_score'] = logsum / count
                row['attention_review'] = 'generic_literature_proxy_not_mining_saturation'
        if not row['data_gate_passed']: row['priority'] = 'hold_insufficient_data'
        elif not b['family_tax_id']: row['priority'] = 'hold_taxonomy_review'
        elif row['attention_score'] is None: row['priority'] = 'hold_literature_missing'
        else: row['priority'] = 'pending_tier'
        rows.append(row)
    rankable = sorted([r for r in rows if r['priority']=='pending_tier'], key=lambda r:(r['attention_score'], -r['eligible_wgs_fastq_samples'], r['family_batch_id']))
    scores = sorted(r['attention_score'] for r in rankable)
    # Rank-based quartiles are explicitly a review heuristic, not a novelty statistic.
    low_end = math.ceil(len(rankable)/4)
    high_start = len(rankable) - math.ceil(len(rankable)/4)
    for i,r in enumerate(rankable):
        r['review_rank'] = i+1
        r['priority'] = 'P1_lower_attention_review' if i < low_end else 'P3_higher_attention_reference' if i >= high_start else 'P2_middle_attention_review'
    ordered = rankable + sorted([r for r in rows if r not in rankable], key=lambda r:(r['priority'], -r['eligible_wgs_fastq_samples'], r['family_batch_id']))
    for filename, records in [('batch_priorities.jsonl', ordered), ('literature_evidence.jsonl', audits)]:
        with (OUT / filename).open('w', encoding='utf-8') as stream:
            for r in records: stream.write(json.dumps(r, ensure_ascii=False) + '\n')
    summary = dict(batch_count=len(rows), min_samples=MIN_SAMPLES, min_studies=MIN_STUDIES,
                   priority_counts=dict(collections.Counter(r['priority'] for r in rows)),
                   literature_queries=len(units), successful_queries=sum(r.get('status')=='ok' for r in results.values()),
                   thresholds='Relative quartiles of sample-weighted log10(1 + PubMed title/abstract hits) for qualifying families',
                   lower_attention_group_count=low_end, generated_utc=datetime.now(timezone.utc).isoformat(),
                   exclusion='Species already marked low_reference, currently E. coli, omitted from eligibility counts and sample weights; source inventory retained',
                   caveat='Generic literature attention is not demonstrated saturation of insertion/mobile-element mining. Assembly pairing is not yet verified.',
                   first_review_batches=[{k:r[k] for k in ['family_batch_id','family_name','eligible_wgs_fastq_samples','distinct_studies','largest_study_sample_share','attention_score']} for r in rankable[:10]])
    (OUT / 'summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    labels = {'P1_lower_attention_review':'P1 先核查', 'P2_middle_attention_review':'P2 后续核查',
              'P3_higher_attention_reference':'P3 低优先级参考', 'hold_insufficient_data':'暂缓：数据不足',
              'hold_taxonomy_review':'暂缓：分类待核查', 'hold_literature_missing':'暂缓：文献查询未完成'}
    table = ['# Batch 初步审阅优先级', '', '优先级用于决定先核查哪批元数据与文献，不是已经证明的新颖性排名或可立即下载的最终名单。', '',
             '门槛：扣除 low_reference 成员后至少 100 个有 WGS FASTQ 的去重 sample accession，至少 2 个不同 study accession。E. coli 的原记录保留。', '',
             '文献指标为样本加权的 log10(1+PubMed 命中数)，越低表示此口径下关注度相对较低；不是“未被挖掘”的概率。', '',
             '| 优先级 | 科/待核查分组 | 合格 WGS 样本 | GENOMIC WGS 样本 | 项目数 | 最大项目占比 | 文献关注度指标 | FASTQ GB |',
             '|---|---|---:|---:|---:|---:|---:|---:|']
    for r in ordered:
        vals = [labels[r['priority']], r['family_name'] or r['family_batch_id'], r['eligible_wgs_fastq_samples'], r['eligible_genomic_wgs_samples'],
                r['distinct_studies'], f"{r['largest_study_sample_share']:.1%}" if r['largest_study_sample_share'] is not None else '',
                f"{r['attention_score']:.3f}" if r['attention_score'] is not None else '', f"{r['eligible_fastq_bytes']/1e9:.2f}"]
        table.append('| ' + ' | '.join(str(v).replace('|','\\|') for v in vals) + ' |')
    (OUT / 'PRIORITIES.md').write_text('\n'.join(table) + '\n', encoding='utf-8')
    shortlist = ['# P1：优先核查的 batch', '',
                 '本表是通过至少 100 个 WGS FASTQ 样本、至少 2 个项目门槛后，通用文献关注度相对最低四分之一的 batch。旧名已纳入检索，但仍不是插入元件研究新颖性的证明。', '',
                 '| 科级 batch | WGS 样本 | 其中 GENOMIC | 项目数 | 最大项目占比 | 需核查事项 |',
                 '|---|---:|---:|---:|---:|---|']
    flag_labels = {'fewer_than_100_GENOMIC_WGS_samples':'GENOMIC 样本不足 100，核查其他文库来源',
                   'over_80pct_samples_in_one_study':'超过 80% 样本集中于一个项目'}
    for r in ordered:
        if r['priority'] != 'P1_lower_attention_review': continue
        notes = '；'.join(flag_labels[f] for f in r['review_flags']) or '待核查 assembly 配对及主题文献'
        shortlist.append(f"| {r['family_name']} | {r['eligible_wgs_fastq_samples']} | {r['eligible_genomic_wgs_samples']} | {r['distinct_studies']} | {r['largest_study_sample_share']:.1%} | {notes} |")
    shortlist += ['', '全部 batch、P2/P3 和暂缓项见 [完整优先级表](PRIORITIES.md)；规则与证据见 [说明](README.md)。', '']
    (OUT / 'SHORTLIST.md').write_text('\n'.join(shortlist), encoding='utf-8')
    readme = ['# 数据量与文献关注度：初步优先级', '',
              '用户确认最少 100 个样本。此处按有 WGS FASTQ 的 sample accession 去重；另保留至少 2 个来源项目作为检查条件，项目数量并不证明实验独立。', '',
              '## 规则', '',
              '1. 从每个 batch 的计数中扣除已标记 low_reference 的物种（目前仅 E. coli），但不删除底表记录。',
              '2. 达不到 100 个 WGS FASTQ 样本或 2 个不同项目的 batch 暂缓；没有科归属的单列待核查。',
              '3. 对数据达标的科，查询其实际有合格样本的各属在 Europe PMC 中 PubMed 来源（SRC:MED）的标题/摘要命中文献数。检索当前属名、NCBI 属级同义名，并加入本批物种在旧属下的完整同义名称；保留完整查询供审计。旧属拆分时不把整个旧属的文献都归给新属。',
              '4. 计算每科各属 log10(1+命中文献数) 的样本加权平均。同一 sample 跨属时均分其权重。避免一个少量出现的热门属直接代表整科，也避免新科名文献少造成虚假低关注度。',
              '5. 仅在数据达标且文献查询完整的科之间，按关注度从低到高排列：最低四分之一 P1、居中一半 P2、最高四分之一 P3。并列值按样本数降序。这是初始审阅规则，不是用户给定的生物学阈值。', '',
              '## 如何解读', '',
              'P1 表示相对较低的通用文献关注度、且达到数据门槛，值得先核查；不等于“该类群没被挖过”。插入元件相关研究或系统筛查的饱和程度本轮未验证。', '',
              'PubMed 名称检索会受历史改名、歧义、数据库收录和摘要措辞影响。物种旧名有助于减少改名造成的漏检；只在旧属名下发表但未提及所检索物种的文献仍可能漏检，属级同义名仍可能导致高估。全部词表和查询均可检查。', '',
              '样本数不等于独立菌株数；WGS 元数据不保证样本为纯培养，因此另列 GENOMIC WGS 数量。最大项目占比反映采样集中程度，只展示不额外改变优先级。', '',
              'reads 长度、覆盖度、测序质量和同一 BioSample 的 assembly 配对尚未验证；100 个 WGS 样本不代表已有 100 对 reads + assembly。', '',
              '科级 batch 是组织容器；保留属、物种和样本归属，不混合不同样本的文件，也不以优先级改变原样本记录。', '',
              '## 输出', '',
              '- [P1 优先核查清单](SHORTLIST.md)：同时提示文库来源和项目集中程度。',
              '- [全部 batch 优先级](PRIORITIES.md)',
              '- [机器可读结果](batch_priorities.jsonl)',
              '- [文献查询与加权证据](literature_evidence.jsonl)',
              '- [查询词及旧名](query_plan.json)',
              '- [经名称数据库核实的补充旧名](curated_aliases.json)：另将 NCBI 同义名称中的命名人/年份去除后检索，避免整句匹配漏检。',
              '- [汇总](summary.json)', '',
              '[Europe PMC 官方 API 说明](https://europepmc.org/RestfulWebService)', '']
    (OUT / 'README.md').write_text('\n'.join(readme), encoding='utf-8')
    assert len(rows) == 518
    assert all(r['eligible_genomic_wgs_samples'] <= r['eligible_wgs_fastq_samples'] <= r['all_wgs_fastq_samples'] for r in rows)
    assert all(r['eligible_wgs_fastq_samples'] >= MIN_SAMPLES and r['distinct_studies'] >= MIN_STUDIES and r['literature_coverage'] > .999999 for r in rankable)
    print(json.dumps(summary, ensure_ascii=False, indent=2), flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--fetch-literature', action='store_true')
    args = p.parse_args()
    build(args.fetch_literature)
