"""Review tables and U/core/D case figures for one graph-entry subset run.

Reads the scratch bundle/calls of core_context.py plus the subset manifest (FNA/GFF3 paths)
and writes small, committable outputs. GFF3 genes are drawn and overlapped only as a boundary
diagnostic; this is not the S4 coding check, and annotation_status stays not_reviewed.
Case pairs are chosen by a fixed rule, not by eye: per relation, the pair from the group with
the most samples (then most occurrences, then id) whose two occurrences are in different samples
if any exist, then the smallest occurrence ids.
"""
import argparse
import csv
import gzip
import json
import re
import shutil
from collections import Counter, defaultdict
from pathlib import Path

RELATIONS = ['similar_flanks_candidate_same_context', 'different_adjacency_low_flank_similarity_candidate',
             'context_unresolved', 'core_similarity_below_threshold', 'edge_or_missing_context']
SHOWN = RELATIONS[:3]
WINDOW = 4000      # bp of genes drawn on each side of the core
END_TOL = 20       # a core end within this many bp of a CDS end counts as at the CDS end


def read_tsv(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter='\t'))


def write_tsv(path, rows, columns):
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, columns, delimiter='\t', lineterminator='\n', extrasaction='ignore')
        w.writeheader(); w.writerows(rows)


def read_cds(gff):
    out = defaultdict(list)
    with open(gff) as f:
        for line in f:
            if line.startswith('#'):
                continue
            c = line.rstrip('\n').split('\t')
            if len(c) < 9 or c[2] not in ('CDS', 'rRNA', 'tRNA'):
                continue
            a = dict(kv.split('=', 1) for kv in c[8].split(';') if '=' in kv)
            name = a.get('gene') or a.get('locus_tag') or ''
            product = re.sub('%2C', ',', a.get('product', ''))
            pseudo = a.get('pseudo') == 'true' or 'partial' in a
            out[c[0]].append((int(c[3]) - 1, int(c[4]), c[6], c[2], name, product, pseudo))
    return {k: sorted(set(v)) for k, v in out.items()}


def boundary_flags(x, cds, contig_len):
    s, e = int(x['start0']), int(x['end0'])
    over = [g for g in cds.get(x['contig'], []) if g[0] < e and g[1] > s]
    inside = [g for g in over if g[0] >= s - END_TOL and g[1] <= e + END_TOL]
    # a feature is cut when a core end falls inside it, away from its ends
    cut = [g for g in over if g not in inside and
           any(g[0] + END_TOL < p < g[1] - END_TOL for p in (s, e))]
    cds_inside = [g for g in inside if g[3] == 'CDS']
    return dict(features_overlapping=len(over), complete_features_inside=len(inside),
                cds_inside=len(cds_inside), cds_inside_pseudo_or_partial=sum(g[6] for g in cds_inside),
                features_cut_by_core_end=len(cut),
                at_contig_end=s < 500 or contig_len - e < 500,
                inside_products='|'.join(sorted({g[5] for g in cds_inside}))[:300])


def svg_case(cases, occ, cds, contig_len, title, subtitle):
    W, row_h, x0, scale = 1180, 70, 330, 1180 - 330 - 20
    y = 70
    parts = []
    colors = dict(core='#edc06f', up='#b6d1df', down='#b4d4b9')
    for rel, pair in cases:
        parts.append(f'<text x="20" y="{y}" font-size="14" font-weight="bold">{rel}</text>')
        if pair is None:
            parts.append(f'<text x="20" y="{y+22}" font-size="12">no pair with this relation</text>'); y += 50; continue
        ids = (f"core {float(pair['core_identity']):.3f}  up {float(pair['upstream_identity'] or 0):.3f}  "
               f"down {float(pair['downstream_identity'] or 0):.3f}  group {pair['family_id']}")
        parts.append(f'<text x="20" y="{y+18}" font-size="12">{ids}</text>')
        y += 28
        for key in ('occurrence_a', 'occurrence_b'):
            x = occ[pair[key]]
            s, e, strand = int(x['start0']), int(x['end0']), x['strand']
            lo, hi = max(0, s - WINDOW), min(contig_len[(x['sample'], x['contig'])], e + WINDOW)
            span = 2 * WINDOW + (e - s)
            px = lambda p: x0 + scale * ((p - lo) if strand == '+' else (hi - p)) / span
            parts.append(f'<text x="20" y="{y+16}" font-size="11">{x["sample"]} {x["contig"]}:{s+1}-{e} ({strand})</text>')
            parts.append(f'<text x="20" y="{y+30}" font-size="10" fill="#555">{x["occurrence_id"]} len {e-s}'
                         f'{"" if x["flank_complete"] == "True" else "  edge-limited"}</text>')
            for a, b, c in ((max(lo, s - 500), s, 'up'), (s, e, 'core'), (e, min(hi, e + 500), 'down')):
                if strand == '-':
                    c = {'up': 'down', 'down': 'up'}.get(c, c)
                l, r = sorted((px(a), px(b)))
                parts.append(f'<rect x="{l:.1f}" y="{y+4}" width="{max(r-l,1):.1f}" height="12" fill="{colors[c]}"/>')
            for g in cds.get(x['contig'], []):
                if g[1] <= lo or g[0] >= hi:
                    continue
                l, r = sorted((px(max(g[0], lo)), px(min(g[1], hi))))
                fwd = (g[2] == '+') == (strand == '+')
                fill = '#999' if g[6] else ('#c0504d' if g[3] != 'CDS' else '#4f6d8f')
                tip = 5 if r - l > 8 else 0
                pts = (f'{l:.1f},{y+20} {r-tip:.1f},{y+20} {r:.1f},{y+26} {r-tip:.1f},{y+32} {l:.1f},{y+32}' if fwd else
                       f'{l+tip:.1f},{y+20} {r:.1f},{y+20} {r:.1f},{y+32} {l+tip:.1f},{y+32} {l:.1f},{y+26}')
                label = (g[4] or g[5])[:40].replace('&', '&amp;').replace('<', '&lt;')
                parts.append(f'<polygon points="{pts}" fill="{fill}" opacity="0.85"><title>{label}</title></polygon>')
            y += row_h - 30
        y += 22
    head = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{y+20}">',
            '<rect width="100%" height="100%" fill="white"/>',
            f'<text x="20" y="26" font-size="18">{title}</text>',
            f'<text x="20" y="46" font-size="12">{subtitle}</text>']
    return '\n'.join(head + parts + ['</svg>']) + '\n'


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run', type=Path, required=True, help='scratch subset dir with bundle/ and calls/')
    p.add_argument('--manifest', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    calls, bundle = a.run / 'calls', a.run / 'bundle'
    manifest = read_tsv(a.manifest)
    samples = {m['sample']: m for m in manifest}
    occs = read_tsv(calls / 'occurrences.tsv')
    occ = {x['occurrence_id']: x for x in occs}
    pairs = read_tsv(calls / 'context_pairs.tsv')
    families = read_tsv(calls / 'families.tsv')
    contig_len = {(r['sample'], r['contig']): int(r['length']) for r in read_tsv(bundle / 'contigs.tsv')}
    cds = {}
    for s, m in samples.items():
        for k, v in read_cds(m['gff3']).items():
            cds[k] = v

    graph = json.loads((bundle / 'graph.json').read_text())
    per_block = Counter(str(n['block_id']) for n in graph['nodes'].values())
    block_len = {}
    for n in graph['nodes'].values():
        s, e = n['position']
        block_len.setdefault(str(n['block_id']), []).append(e - s)
    repeated = {b for b, c in per_block.items() if c >= 2}
    def length_class(b):
        ls = block_len[b]
        if min(ls) >= 2500: return 'all_ge_2500'
        if max(ls) < 50: return 'all_lt_50'
        return 'some_in_range'
    repeat_len_class = Counter(length_class(b) for b in repeated)
    occ_json = {x['occurrence_id']: x for x in json.loads((calls / 'occurrences.json').read_text())}

    flags = []
    for x in occs:
        f = boundary_flags(x, cds, contig_len[(x['sample'], x['contig'])])
        anchors = occ_json[x['occurrence_id']]['anchors']
        f['adjacent_repeated_blocks'] = sum(1 for n in anchors if n is not None and n[0] in repeated)
        f.update(occurrence_id=x['occurrence_id'], family_id=x['family_id'], sample=x['sample'], core_length=x['core_length'],
                 flank_complete=x['flank_complete'])
        flags.append(f)
    write_tsv(a.out / 'boundary_flags.tsv', flags, ['occurrence_id', 'family_id', 'sample', 'core_length', 'flank_complete',
        'at_contig_end', 'adjacent_repeated_blocks', 'features_overlapping', 'complete_features_inside', 'cds_inside',
        'cds_inside_pseudo_or_partial', 'features_cut_by_core_end', 'inside_products'])

    by_fam = {r['family_id']: r for r in families}
    cases, case_rows = [], []
    for rel in RELATIONS:
        xs = [q for q in pairs if q['relation'] == rel]
        if not xs:
            cases.append((rel, None)); continue
        cross = [q for q in xs if occ[q['occurrence_a']]['sample'] != occ[q['occurrence_b']]['sample']]
        pool = cross or xs
        pick = min(pool, key=lambda q: (-int(by_fam[q['family_id']]['samples']), -int(by_fam[q['family_id']]['occurrences']),
                                        q['family_id'], q['occurrence_a'], q['occurrence_b']))
        cases.append((rel, pick))
        case_rows.append(dict(pick, cross_sample=bool(cross), relation_pairs=len(xs)))
    write_tsv(a.out / 'case_pairs.tsv', case_rows, ['relation', 'family_id', 'occurrence_a', 'occurrence_b', 'core_identity',
        'upstream_identity', 'downstream_identity', 'cross_sample', 'relation_pairs'])
    name = a.manifest.name.replace('.samples.tsv', '')
    (a.out / 'cases.svg').write_text(svg_case([c for c in cases if c[0] in SHOWN], occ, cds, contig_len,
        f'{name}: repeated graph blocks, U/core/D review cases',
        'Assembly only; graph block boundaries provisional; GFF3 features drawn for review (blue CDS, red RNA, grey pseudo/partial), '
        'not a coding check. 500 bp flanks coloured; core in orange.'))

    summary = json.loads((calls / 'summary.json').read_text())
    prov = json.loads((bundle / 'graph_provenance.json').read_text())
    times = {}
    for k in ('build', 'call'):
        t = (a.run / f'{k}.time').read_text()
        wall = re.search(r'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): (\S+)', t).group(1)
        rss = int(re.search(r'Maximum resident set size \(kbytes\): (\d+)', t).group(1))
        times[k] = dict(wall=wall, max_rss_gb=round(rss / 1048576, 2))
    multi = Counter((x['family_id'], x['sample']) for x in occs)
    cross_sample_fams = sum(1 for r in families if int(r['samples']) >= 2)
    # Copy structure separates multi-copy cores from blocks present once per genome (mostly shared backbone)
    multicopy = {f for (f, _), c in multi.items() if c >= 2}
    by_structure = defaultdict(Counter)
    for q in pairs:
        same = occ[q['occurrence_a']]['sample'] == occ[q['occurrence_b']]['sample']
        key = ('multicopy_group' if q['family_id'] in multicopy else 'single_copy_per_sample_group') + \
              (':same_sample_pair' if same else ':cross_sample_pair')
        by_structure[key][q['relation']] += 1
    diff_groups = {q['family_id'] for q in pairs if q['relation'] == RELATIONS[1]}
    review = dict(subset=name, version=summary['version'], evidence_source='assembly_only', production_ready=False,
        samples=len(manifest), contigs=len(contig_len), total_bp=sum(contig_len.values()),
        pangraph=prov['tool_version'], binary_sha256=prov.get('binary_sha256'), graph_verify_passed=prov['graph_verify_passed'],
        graph_sha256=prov['graph_sha256'], build_seconds=prov.get('build_seconds'), cost=times,
        graph_blocks=len(per_block), repeated_blocks_any_length=len(repeated), repeated_block_length_classes=dict(repeat_len_class),
        repeated_block_groups_in_range=summary['families'], groups_spanning_2plus_samples=cross_sample_fams,
        groups_single_sample_only=summary['families'] - cross_sample_fams,
        occurrences=summary['occurrences'], occurrences_second_copy_same_sample=sum(c - 1 for c in multi.values()),
        context_pairs=summary['context_pairs'], relations={r: summary['relations'].get(r, 0) for r in RELATIONS},
        groups_multicopy_in_some_sample=len(multicopy),
        groups_with_different_context_candidate=len(diff_groups),
        groups_with_different_context_candidate_multicopy=len(diff_groups & multicopy),
        relations_by_copy_structure={k: {r: v.get(r, 0) for r in RELATIONS} for k, v in sorted(by_structure.items())},
        truncated_pair_reviews=summary['truncated_pair_reviews'], excluded_circular_nodes=summary['excluded_nodes'],
        boundary=dict(
            edge_limited_occurrences=sum(f['flank_complete'] != 'True' for f in flags),
            occurrences_within_500bp_of_contig_end=sum(f['at_contig_end'] for f in flags),
            occurrences_next_to_another_repeated_block=sum(f['adjacent_repeated_blocks'] > 0 for f in flags),
            occurrences_with_a_feature_cut_by_core_end=sum(f['features_cut_by_core_end'] > 0 for f in flags),
            occurrences_with_no_feature_overlap=sum(f['features_overlapping'] == 0 for f in flags),
            occurrences_with_1_2_complete_cds_inside_and_no_cut=sum(1 <= f['cds_inside'] <= 2 and f['features_cut_by_core_end'] == 0
                                                                    for f in flags),
            note='GFF3 overlap preview only; S4 coding checks not applied; graph blocks are not complete insertions'),
        case_pair_rule=__doc__.split('Case pairs')[1].strip(),
        annotation_complete=False, site_groups_computed=False)
    (a.out / 'review_summary.json').write_text(json.dumps(review, indent=2) + '\n')
    for f in ('summary.json', 'families.tsv', 'occurrences.tsv', 'excluded_nodes.json'):
        shutil.copy(calls / f, a.out / f)
    shutil.copy(bundle / 'graph_provenance.json', a.out / 'graph_provenance.json')
    shutil.copy(bundle / 'contigs.tsv', a.out / 'contigs.tsv')
    with open(calls / 'context_pairs.tsv', 'rb') as s, gzip.GzipFile(a.out / 'context_pairs.tsv.gz', 'wb', mtime=0) as d:
        shutil.copyfileobj(s, d)
    print(json.dumps(review, indent=2))


if __name__ == '__main__':
    main()
