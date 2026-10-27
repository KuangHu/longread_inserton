"""Review of one protein-entry subset run: product-keyword hint, overlap with the graph entry, gene-order case figure.

Product keywords are a reading aid only (not a mobile-element classification). Case cores are fixed by rule, in the
target_lt2500 scope: single-gene cores with >= 1 different_context_candidate pair, ranked by context classes linked by those pairs, then
samples, then core key; the top three are drawn with one occurrence per context class (up to four), plus the
similar-context core with the most samples as a reference. Neighbour genes are coloured by homology group.
"""
import argparse
import csv
import gzip
import hashlib
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

KEYWORDS = [('transposase_or_IS', r'transposase|\bIS\d|insertion element|IstB'),
            ('integrase_recombinase_resolvase', r'integrase|recombinase|resolvase|invertase'),
            ('phage_like', r'phage|capsid|portal|terminase'), ('intron_rt', r'intron|maturase'),
            ('hypothetical', r'hypothetical'), ('translation_or_ribosomal', r'elongation factor|ribosomal')]


def hint(product):
    for k, r in KEYWORDS:
        if re.search(r, product, re.I):
            return k
    return 'other'


def read_tsv(path):
    op = gzip.open if str(path).endswith('.gz') else open
    with op(path, 'rt') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def colour(label):
    if not label or label == 'nolabel':
        return '#dddddd'
    h = int(hashlib.sha256(label.encode()).hexdigest()[:6], 16)
    return 'hsl(%d,55%%,%d%%)' % (h % 360, 55 + (h >> 9) % 20)


def draw(cases, title, subtitle):
    W, y, box, gap = 1180, 70, 92, 6
    parts = []
    for heading, rows in cases:
        parts.append('<text x="20" y="%d" font-size="13" font-weight="bold">%s</text>' % (y, esc(heading)))
        y += 12
        for o in rows:
            up = o['up_labels'].split(',')[::-1] if o['up_labels'] else []
            down = o['down_labels'].split(',') if o['down_labels'] else []
            parts.append('<text x="20" y="%d" font-size="11">%s %s:%s-%s (%s) class %s</text>' % (
                y + 18, o['sample'], o['contig'], int(o['start0']) + 1, o['end0'], o['strand'], o['context_class_target_lt2500'].split('|')[-1]))
            x = 360
            for i, lab in enumerate(up + ['CORE'] + down):
                rev = lab.endswith('(rev)')
                base = lab.replace('(rev)', '')
                fill = '#edc06f' if lab == 'CORE' else colour(base)
                pts = ('%d,%d %d,%d %d,%d %d,%d %d,%d' % (x, y + 6, x + box - 10, y + 6, x + box, y + 15, x + box - 10, y + 24, x, y + 24)
                       if not rev else '%d,%d %d,%d %d,%d %d,%d %d,%d' % (x + 10, y + 6, x + box, y + 6, x + box, y + 24, x + 10, y + 24, x, y + 15))
                parts.append('<polygon points="%s" fill="%s" stroke="#666" stroke-width="0.5"><title>%s</title></polygon>' % (pts, fill, esc(base)))
                parts.append('<text x="%d" y="%d" font-size="8">%s</text>' % (x + 12, y + 18, esc(base[-9:] if lab != 'CORE' else 'core')))
                x += box + gap
            y += 30
        y += 16
    head = ['<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d">' % (W, y + 10),
            '<rect width="100%" height="100%" fill="white"/>',
            '<text x="20" y="26" font-size="18">%s</text>' % esc(title),
            '<text x="20" y="46" font-size="12">%s</text>' % esc(subtitle)]
    return '\n'.join(head + parts + ['</svg>']) + '\n'


def esc(s):
    return s.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--subset-dir', type=Path, required=True, help='committed subset dir with protein_entry/ and graph outputs')
    a = p.parse_args()
    d, pe = a.subset_dir, a.subset_dir / 'protein_entry'
    cores = read_tsv(pe / 'protein_cores.tsv')
    occs = read_tsv(pe / 'protein_occurrences.tsv.gz')
    by_core = defaultdict(list)
    for o in occs:
        by_core[('all', o['structure'], o['core_key'])].append(o)
        if o['in_target_length'] == 'True':
            by_core[('target_lt2500', o['structure'], o['core_key'])].append(o)

    # graph-entry overlap: graph occurrences in groups with a different-context candidate pair
    gocc = {r['occurrence_id']: r for r in read_tsv(d / 'occurrences.tsv')}
    gdiff = {r['family_id'] for r in read_tsv(d / 'context_pairs.tsv.gz')
             if r['relation'] == 'different_adjacency_low_flank_similarity_candidate'}
    giv = defaultdict(list)
    for r in gocc.values():
        giv[(r['sample'], r['contig'])].append((int(r['start0']), int(r['end0']), r['family_id'] in gdiff))

    def graph_state(o):
        s, e = int(o['start0']), int(o['end0'])
        hits = [x for x in giv[(o['sample'], o['contig'])] if x[0] < e and x[1] > s]
        return 'graph_diff_group' if any(h[2] for h in hits) else ('graph_block_other' if hits else 'no_graph_block_50_2500')

    rows, table = [], Counter()
    for c in cores:
        diff = int(c['different_context_candidate']) > 0
        states = Counter(graph_state(o) for o in by_core[(c['scope'], c['structure'], c['core_key'])])
        c['product_hint'] = hint(c['products'])
        c['graph_overlap'] = ','.join('%s=%d' % kv for kv in sorted(states.items()))
        if diff:
            table[(c['scope'], c['structure'], c['product_hint'])] += 1
            rows.append(c)
    cols = ['scope', 'structure', 'core_key', 'occurrences', 'samples', 'max_copies_in_one_sample', 'context_classes_provisional',
            'context_classes_linked_by_different_candidate', 'different_context_candidate', 'similar_context',
            'one_side_shared_unresolved', 'insufficient_context_labels', 'span_range', 'spacing_range', 'product_hint', 'graph_overlap', 'products']
    rows.sort(key=lambda c: (c['scope'], c['structure'], -int(c['context_classes_linked_by_different_candidate']), -int(c['samples']), c['core_key']))
    with open(pe / 'different_context_cores.tsv', 'w', newline='') as f:
        w = csv.DictWriter(f, cols, delimiter='\t', lineterminator='\n', extrasaction='ignore')
        w.writeheader(); w.writerows(rows)

    overlap = Counter()
    for c in rows:
        if c['structure'] != 'single_gene':
            continue
        got = 'graph_diff_group' in c['graph_overlap']
        overlap[(c['scope'], 'any_occurrence_in_graph_diff_group' if got else 'none_in_graph_diff_group', c['product_hint'] == 'transposase_or_IS')] += 1

    single = [c for c in rows if c['structure'] == 'single_gene' and c['scope'] == 'target_lt2500']
    cases = []
    for c in single[:3]:
        xs, seen = [], set()
        cls = 'context_class_target_lt2500'
        for o in sorted(by_core[('target_lt2500', 'single_gene', c['core_key'])], key=lambda o: (o[cls], o['sample'], o['occurrence_id'])):
            if o[cls] not in seen and len(xs) < 4:
                seen.add(o[cls]); xs.append(o)
        cases.append(('different-context candidate %s: %s (%s occ, %s samples, %s classes)' % (
            c['core_key'], c['products'][:60], c['occurrences'], c['samples'], c['context_classes_provisional']), xs))
    ref = [c for c in cores if c['scope'] == 'target_lt2500' and c['structure'] == 'single_gene' and int(c['similar_context'])
           and int(c['similar_context']) == int(c['pairs_checked'])]
    if ref:
        r = min(ref, key=lambda c: (-int(c['samples']), c['core_key']))
        cases.append(('similar-context reference %s: %s' % (r['core_key'], r['products'][:60]),
                      sorted(by_core[('target_lt2500', 'single_gene', r['core_key'])], key=lambda o: o['sample'])[:3]))
    name = d.name
    (pe / 'cases_protein.svg').write_text(draw(cases, '%s: protein entry, target <2500 bp, oriented neighbour genes (K=3)' % name,
        'Assembly only. Arrows = genes in core orientation, coloured by homology group (same colour = homologous); '
        'grey = no protein label. Not a boundary or orthology call.'))
    review = dict(subset=name, product_hint_by_scope_structure={'%s:%s:%s' % k: v for k, v in sorted(table.items())},
                  single_gene_graph_overlap={'%s:%s:%s' % (sc, k, 'IS_hint' if t else 'non_IS_hint'): v
                                             for (sc, k, t), v in sorted(overlap.items())},
                  case_rule=__doc__.split('Case cores')[1].strip(), cases=[c[0] for c in cases])
    (pe / 'protein_review.json').write_text(json.dumps(review, indent=2) + '\n')
    print(json.dumps(review, indent=2))


if __name__ == '__main__':
    main()
