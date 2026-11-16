"""Review layer on top of S6 v0.1 (site_groups.py outputs are read, never changed). Evidence stays assembly_only;
production_ready=false; site groups and position units are not independent insertion events.

1. Boundary denominators. Every S6 occurrence is one of: both sides resolved (candidate conserved ends), reviewed
   with an unknown side, or not reviewed (no different-context partner in core_review). Stable ends are further split
   by how many valid partners and how many distinct position units those partners come from; one valid partner only
   means one comparison supports the end.

2. Different coordinates vs different host context. A non-overlapping same-assembly pair always has different
   coordinates; its host context (S6 host windows, N=5 per side) is classified separately:
     context_different_both_sides   both sides 'different' (the target U-core-D pattern)
     context_repeated_both_sides    both sides 'same' (repeated_flank: segment copied with the core, or a repeat)
     context_repeated_one_side      one side 'same', the other 'different'
     context_unknown                any side unknown
   A flank DNA cross-check (blastn of the outer flank, INNER..W bp beyond each core end, between the two copies) is
   reported next to the class and never changes it.

3. Relations between position units. A unit is a site group or a single-observation site (unresolved occurrences are
   not units). Two units of one core are compared through every cross pair of members:
     ambiguous     some cross pair is same_position_candidate (only possible as a multiple-correspondence link)
     conflicting   some cross pair has conflicting_sides
     different     every cross pair is a different-position verdict
     unknown       anything else
   A 'different' relation also carries the context class of its cross pairs. Group counts are not counts of
   confirmed distinct positions: per core we report the largest set of units that are pairwise 'different' (and
   pairwise different with both-side context different), never the plain unit count.

4. Representative cases (rules: S6_SITE_GROUPS_v0.1.md, review layer), fixed rules over the existing assemblies: context different, repeated context, conserved
   paralog positions, multiple correspondence (see pick()).
"""
import argparse
import json
from collections import Counter, defaultdict
from pathlib import Path

import core_review as cr
import protein_context as pc
from review_protein import colour, esc

VERSION = 's6_review_v0.1'
W = 3000           # flank length per side for the DNA cross-check
INNER = 500        # bp next to each core end left out of the cross-check (element ends may extend past the CDS)
MIN_OUTER = 1000   # shorter outer flanks are reported as short_flank
DIFFERENT = {'different_position_same_assembly', 'different_position_candidate'}
SUPPORTED = {'site_group', 'single_observation_site'}


def context_class(up, down):
    s = {up, down}
    if s == {'different'}:
        return 'context_different_both_sides'
    if s == {'same'}:
        return 'context_repeated_both_sides'
    if s == {'same', 'different'}:
        return 'context_repeated_one_side'
    return 'context_unknown'


def unit_relation(cross):
    """cross: pair rows between the members of two units. Returns (relation, context)."""
    vs = [r['verdict'] for r in cross]
    if 'same_position_candidate' in vs:
        return 'ambiguous', ''
    if 'conflicting_sides' in vs:
        return 'conflicting', ''
    if not vs or not all(v in DIFFERENT for v in vs):
        return 'unknown', ''
    cs = {context_class(r['up_side'], r['down_side']) for r in cross}
    if cs == {'context_different_both_sides'}:
        return 'different', 'context_different_both_sides'
    if cs & {'context_repeated_both_sides', 'context_repeated_one_side'}:
        return 'different', 'context_repeated_in_some_pair'
    return 'different', 'context_unknown_in_some_pair'


def max_clique(nodes, adj):
    """Size of the largest set of nodes that are pairwise adjacent (Bron-Kerbosch with pivot)."""
    best = 0

    def bk(r, p, x):
        nonlocal best
        if not p and not x:
            best = max(best, r)
            return
        if r + len(p) <= best:
            return
        u = max(p | x, key=lambda v: len(adj[v] & p))
        for v in list(p - adj[u]):
            bk(r + 1, p & adj[v], x & adj[v])
            p = p - {v}; x = x | {v}
    bk(0, set(nodes), set())
    return best


def covered(intervals, lo, hi):
    """bp of [lo, hi) covered by the union of [s, e) intervals."""
    n, cur = 0, lo
    for s, e in sorted(intervals):
        s, e = max(s, cur), min(e, hi)
        if e > s:
            n += e - s; cur = e
    return n


def outer(side, n):
    """The outer part of a flank of length n in core orientation: up flanks end at the core, down flanks start there."""
    return (0, max(0, n - INNER)) if side == 'up' else (min(n, INNER), n)


def flank_class(cov, length):
    if length < MIN_OUTER:
        return 'short_flank'
    return 'shared' if cov >= 0.5 else ('not_shared' if cov <= 0.1 else 'partial')


def flank_check(fx, fy):
    """fx, fy: {'up': seq, 'down': seq} in core orientation. Per side: the smaller outer-flank fraction covered by
    blastn HSPs to the other copy's flank on the same side, the shorter outer length, and a descriptive class.
    SUPERSEDED by flank_dna.outer_flank_check (v2): this v0.1 check aligns the whole flanks, so an outer window can
    count as shared while it matches only the other copy's excluded inner window. Kept unchanged so the committed
    dna_* columns of the S6 review outputs stay reproducible; read flank_dna_v2.tsv.gz instead."""
    out = {}
    for side in ('up', 'down'):
        x, y = fx[side], fy[side]
        (xl, xh), (yl, yh) = outer(side, len(x)), outer(side, len(y))
        length = min(xh - xl, yh - yl)
        cov = 0.0
        if length > 0:
            hits = cr.blast('blastn', {'x': x, 'y': y}, 'blastn')
            cx = covered([(h['qs'] - 1, h['qe']) for h in hits.get(('x', 'y'), [])], xl, xh) / (xh - xl)
            cy = covered([(h['qs'] - 1, h['qe']) for h in hits.get(('y', 'x'), [])], yl, yh) / (yh - yl)
            cov = min(cx, cy)
        out[side] = (round(cov, 3), length, flank_class(cov, length))
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--subset-dir', type=Path, required=True)
    p.add_argument('--work', type=Path, required=True, help='protein work dir (genes.tsv.gz)')
    a = p.parse_args()
    s6 = a.subset_dir / 'protein_entry' / 's6'
    out = s6 / 'review'
    out.mkdir(exist_ok=True)
    manifest = a.subset_dir.parent / (a.subset_dir.name + '.samples.tsv')

    occs = {o['occurrence_id']: o for o in pc.read_tsv(s6 / 's6_occurrences.tsv.gz')}
    pairs = pc.read_tsv(s6 / 's6_pairs.tsv.gz')
    pair = {frozenset((r['occurrence_a'], r['occurrence_b'])): r for r in pairs}
    cores = {(c['structure'], c['core_key']): c for c in pc.read_tsv(s6 / 's6_cores.tsv')}
    stable = pc.read_tsv(s6 / 'stable_end_candidates.tsv')
    genes = {g['gene_uid']: g for g in pc.read_tsv(a.work / 'genes.tsv.gz')}
    circular = {(g['sample'], g['contig']): g['circular_contig'] == 'True' for g in genes.values()}

    # flanks in core orientation, for the DNA cross-check
    fasta = {r['sample']: r['fasta'] for r in pc.read_tsv(manifest)}
    for smp in sorted({o['sample'] for o in occs.values()}):
        seqs = pc.read_fna(fasta[smp])
        for o in occs.values():
            if o['sample'] == smp:
                sub, c0, c1 = cr.window(seqs[o['contig']], int(o['start0']), int(o['end0']), o['strand'] == '+', W,
                                        circular[(smp, o['contig'])])
                o['flank'] = dict(up=sub[:c0], down=sub[c1:])

    by_core = defaultdict(list)
    for oid, o in sorted(occs.items()):
        by_core[(o['structure'], o['core_key'])].append(oid)
    unit_of = {x: (o['site_group'] or 'single:' + x) for x, o in occs.items() if o['s6_status'] in SUPPORTED}

    # 2. same-assembly pairs: coordinates vs context
    sa_rows = []
    for r in pairs:
        if r['verdict'] != 'different_position_same_assembly':
            continue
        x, y = occs[r['occurrence_a']], occs[r['occurrence_b']]
        fc = flank_check(x['flank'], y['flank'])
        sa_rows.append(dict(structure=r['structure'], core_key=r['core_key'], sample=r['sample_a'],
            occurrence_a=r['occurrence_a'], occurrence_b=r['occurrence_b'],
            status_a=x['s6_status'], status_b=y['s6_status'],
            both_supported_positions=x['s6_status'] in SUPPORTED and y['s6_status'] in SUPPORTED,
            up_side=r['up_side'], up_support=r['up_support'], down_side=r['down_side'], down_support=r['down_support'],
            context_class=context_class(r['up_side'], r['down_side']), repeated_label=r['repeated_label'],
            dna_outer_cov_up=fc['up'][0], dna_outer_len_up=fc['up'][1], dna_up=fc['up'][2],
            dna_outer_cov_down=fc['down'][0], dna_outer_len_down=fc['down'][1], dna_down=fc['down'][2]))
    pc.write_tsv(out / 's6_same_assembly_pairs.tsv.gz', sa_rows, list(sa_rows[0]))
    sa = {frozenset((r['occurrence_a'], r['occurrence_b'])): r for r in sa_rows}

    # 3. position units and their relations
    unit_rows, rel_rows, core_rows = [], [], []
    for key, ids in sorted(by_core.items()):
        members = defaultdict(list)
        for x in ids:
            if x in unit_of:
                members[unit_of[x]].append(x)
        units = sorted(members)
        for u in units:
            unit_rows.append(dict(unit=u, structure=key[0], core_key=key[1],
                                  type='site_group' if u.startswith('SG|') else 'single_observation_site',
                                  members=len(members[u]), assemblies=';'.join(sorted({occs[x]['sample'] for x in members[u]})),
                                  occurrences=';'.join(members[u])))
        adj_diff, adj_ctx = defaultdict(set), defaultdict(set)
        rels = Counter()
        for i, u in enumerate(units):
            for v in units[i + 1:]:
                cross = [pair[frozenset((x, y))] for x in members[u] for y in members[v]]
                rel, ctx = unit_relation(cross)
                rels[rel] += 1
                if rel == 'different':
                    adj_diff[u].add(v); adj_diff[v].add(u)
                    if ctx == 'context_different_both_sides':
                        adj_ctx[u].add(v); adj_ctx[v].add(u)
                rel_rows.append(dict(structure=key[0], core_key=key[1], unit_a=u, unit_b=v, relation=rel, context=ctx,
                                     cross_pairs=len(cross),
                                     cross_verdicts=','.join('%s=%d' % kv for kv in sorted(Counter(r['verdict'] for r in cross).items())),
                                     same_assembly_cross_pairs=sum(r['sample_a'] == r['sample_b'] for r in cross),
                                     same_assembly_cross_context=','.join('%s=%d' % kv for kv in sorted(Counter(
                                         context_class(r['up_side'], r['down_side']) for r in cross
                                         if r['sample_a'] == r['sample_b']).items()))))
        # same-assembly classes, all pairs and per assembly with >= 2 supported positions
        cls = Counter(r['context_class'] for k, r in sa.items() if (r['structure'], r['core_key']) == key)
        per_asm = defaultdict(set)
        for k, r in sa.items():
            if (r['structure'], r['core_key']) == key and r['both_supported_positions']:
                per_asm[r['sample']].add(r['context_class'])
        best = ''
        if any('context_different_both_sides' in s for s in per_asm.values()):
            best = 'context_different_both_sides'
        elif any(s & {'context_repeated_both_sides', 'context_repeated_one_side'} for s in per_asm.values()):
            best = 'context_repeated'
        elif per_asm:
            best = 'context_unknown'
        c = cores[key]
        core_rows.append(dict(structure=key[0], core_key=key[1], products=c['products'], product_hint=c['product_hint'],
            occurrences=c['occurrences'], assemblies=c['assemblies'], copies_per_assembly=c['copies_per_assembly'],
            site_groups=c['multi_observation_site_groups'], single_observation_sites=c['single_observation_sites'],
            unresolved=c['unresolved'], position_units=len(units),
            unit_pairs_different=rels['different'], unit_pairs_ambiguous=rels['ambiguous'],
            unit_pairs_conflicting=rels['conflicting'], unit_pairs_unknown=rels['unknown'],
            max_pairwise_different_units=max_clique(units, adj_diff) if units else 0,
            max_pairwise_context_different_units=max_clique(units, adj_ctx) if units else 0,
            same_assembly_pairs=sum(cls.values()),
            **{'sa_' + k: cls[k] for k in ('context_different_both_sides', 'context_repeated_both_sides',
                                           'context_repeated_one_side', 'context_unknown')},
            assemblies_2plus_supported_positions=c['assemblies_with_2plus_supported_positions'],
            assemblies_2plus_supported_split=';'.join('%s:%s' % (s, '+'.join(sorted(v))) for s, v in sorted(per_asm.items())),
            core_2plus_supported_class=best))
    pc.write_tsv(out / 's6_position_units.tsv', unit_rows, list(unit_rows[0]))
    pc.write_tsv(out / 's6_unit_relations.tsv.gz', rel_rows, list(rel_rows[0]))
    pc.write_tsv(out / 's6_core_positions.tsv', core_rows, list(core_rows[0]))

    # 1. boundary denominators and stable-end partner diversity
    summary = dict(version=VERSION, subset=a.subset_dir.name, evidence_source='assembly_only', production_ready=False,
                   site_groups_are_independent_insertion_events=False,
                   parameters=dict(flank_bp=W, inner_bp_excluded=INNER, min_outer_bp=MIN_OUTER,
                                   dna_class='shared >= 0.5, not_shared <= 0.1 of the outer flank covered; descriptive only'),
                   structures={})
    for st in ('single_gene', 'two_gene'):
        os_ = [o for o in occs.values() if o['structure'] == st]
        se = [r for r in stable if r['structure'] == st]
        sar = [r for r in sa_rows if r['structure'] == st]
        cs = [c for c in core_rows if c['structure'] == st]
        rs = [r for r in rel_rows if r['structure'] == st]

        def bnd(o):
            if o['boundary_up'] == 'not_reviewed_no_different_partner':
                return 'not_reviewed_no_different_partner'
            return 'both_sides_resolved' if o['boundary_up'] == o['boundary_down'] == 'resolved' else 'reviewed_side_unknown'
        summary['structures'][st] = dict(
            occurrences=len(os_), boundary=dict(Counter(bnd(o) for o in os_)),
            stable_end_partners=dict(
                one_valid_partner=sum(int(r['valid_partners_both']) == 1 for r in se),
                two_plus_partners_one_assembly=sum(int(r['valid_partners_both']) > 1 and int(r['valid_partner_assemblies']) == 1 for r in se),
                two_plus_partner_assemblies=sum(int(r['valid_partner_assemblies']) > 1 for r in se),
                partners_in_two_plus_position_units=sum(int(r['valid_partner_site_groups']) > 1 for r in se),
                partners_in_no_position_unit=sum(int(r['valid_partner_site_groups']) == 0 for r in se)),
            same_assembly_pairs=dict(Counter(r['context_class'] for r in sar)),
            same_assembly_pairs_both_supported=dict(Counter(r['context_class'] for r in sar if r['both_supported_positions'])),
            same_assembly_dna_vs_class={k: dict(Counter('%s/%s' % (r['dna_up'], r['dna_down']) for r in sar if r['context_class'] == k))
                                        for k in sorted({r['context_class'] for r in sar})},
            cores_with_assembly_2plus_supported_positions=dict(Counter(c['core_2plus_supported_class'] for c in cs
                                                                       if c['core_2plus_supported_class'])),
            position_units=sum(c['position_units'] for c in cs),
            unit_pairs=dict(Counter(r['relation'] for r in rs)),
            unit_pairs_different_context=dict(Counter(r['context'] for r in rs if r['relation'] == 'different')),
            group_group_pairs=dict(Counter(r['relation'] for r in rs if r['unit_a'].startswith('SG|') and r['unit_b'].startswith('SG|'))),
            sum_max_pairwise_different_units=sum(c['max_pairwise_different_units'] for c in cs),
            sum_max_pairwise_context_different_units=sum(c['max_pairwise_context_different_units'] for c in cs),
            cores_with_2plus_pairwise_context_different_units=sum(c['max_pairwise_context_different_units'] >= 2 for c in cs))

    cases = pick(sa_rows, rel_rows, core_rows, pairs, occs, unit_of)
    case_rows = []
    for c in cases:
        for i, x in enumerate(c['rows']):
            for y in c['rows'][i + 1:]:
                r = pair[frozenset((x, y))]
                k = frozenset((x, y))
                fc = (lambda s: dict(up=(s['dna_outer_cov_up'], s['dna_up']), down=(s['dna_outer_cov_down'], s['dna_down'])))(sa[k]) \
                    if k in sa else {sd: (v[0], v[2]) for sd, v in flank_check(occs[x]['flank'], occs[y]['flank']).items()}
                case_rows.append(dict(case=c['case'], structure=c['key'][0], core_key=c['key'][1], occurrence_a=x, occurrence_b=y,
                    sample_a=occs[x]['sample'], sample_b=occs[y]['sample'], unit_a=unit_of.get(x, 'unresolved'),
                    unit_b=unit_of.get(y, 'unresolved'), verdict=r['verdict'], up_side=r['up_side'], up_support=r['up_support'],
                    down_side=r['down_side'], down_support=r['down_support'],
                    ambiguous_multiple_correspondence=r['ambiguous_multiple_correspondence'],
                    context_class=context_class(r['up_side'], r['down_side']) if r['sample_a'] == r['sample_b']
                    else 'cross_assembly_host_sides:%s/%s' % (r['up_side'], r['down_side']),
                    dna_outer_cov_up=fc['up'][0], dna_up=fc['up'][1], dna_outer_cov_down=fc['down'][0], dna_down=fc['down'][1]))
    pc.write_tsv(out / 'review_cases.tsv', case_rows, list(case_rows[0]) if case_rows else ['case'])
    summary['cases'] = [dict(case=c['case'], core_key=c['key'][1], occurrences=c['rows']) for c in cases]
    summary['case_rule'] = pick.__doc__
    (out / 's6_review_summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    (out / 'cases_s6_review.svg').write_text(draw(cases, case_rows, occs, cores, genes, unit_of, a.subset_dir.name))
    print(json.dumps(summary, indent=2))


def pick(sa_rows, rel_rows, core_rows, pairs, occs, unit_of):
    """Fixed rules, single-gene cores, one case each when available:
    context_different: a same-assembly pair of supported positions with both sides different; prefer DNA not_shared on
      both sides, then a mobile-keyword core, then more core occurrences.
    repeated_context: a same-assembly pair with both sides same; prefer DNA shared on both sides, then larger support.
    conserved_paralog_positions: two site groups of a non-mobile core related 'different' with both-side context
      different; the pair maximising the assemblies shared by both groups, then the smaller group size. Shows the two
      copies in up to two shared assemblies.
    multiple_correspondence: an occurrence with ambiguous links to >= 2 copies in one other assembly; most links first.
    Ties break on identifiers."""
    hint = {(c['structure'], c['core_key']): c['product_hint'] for c in core_rows}
    nocc = {(c['structure'], c['core_key']): int(c['occurrences']) for c in core_rows}
    out = []
    sg = [r for r in sa_rows if r['structure'] == 'single_gene']
    xs = [r for r in sg if r['both_supported_positions'] and r['context_class'] == 'context_different_both_sides']
    if xs:
        r = min(xs, key=lambda r: (not (r['dna_up'] == r['dna_down'] == 'not_shared'),
                                   hint[('single_gene', r['core_key'])] not in cr.MOBILE_HINTS,
                                   -nocc[('single_gene', r['core_key'])], r['occurrence_a'], r['occurrence_b']))
        out.append(dict(case='context_different', key=('single_gene', r['core_key']), rows=[r['occurrence_a'], r['occurrence_b']]))
    xs = [r for r in sg if r['context_class'] == 'context_repeated_both_sides']
    if xs:
        r = min(xs, key=lambda r: (not (r['dna_up'] == r['dna_down'] == 'shared'),
                                   -(int(r['up_support']) + int(r['down_support'])), r['occurrence_a'], r['occurrence_b']))
        out.append(dict(case='repeated_context', key=('single_gene', r['core_key']), rows=[r['occurrence_a'], r['occurrence_b']]))
    best = None
    for r in rel_rows:
        if r['structure'] != 'single_gene' or not (r['unit_a'].startswith('SG|') and r['unit_b'].startswith('SG|')) \
                or r['relation'] != 'different' or r['context'] != 'context_different_both_sides' \
                or hint[('single_gene', r['core_key'])] in cr.MOBILE_HINTS:
            continue
        ma = {occs[x]['sample']: x for x, u in unit_of.items() if u == r['unit_a']}
        mb = {occs[x]['sample']: x for x, u in unit_of.items() if u == r['unit_b']}
        shared = sorted(set(ma) & set(mb))
        score = (-len(shared), -min(len(ma), len(mb)), r['unit_a'], r['unit_b'])
        if shared and (best is None or score < best[0]):
            best = (score, r, [z for s in shared[:2] for z in (ma[s], mb[s])])
    if best:
        out.append(dict(case='conserved_paralog_positions', key=('single_gene', best[1]['core_key']), rows=best[2]))
    links = defaultdict(lambda: defaultdict(list))
    for r in pairs:
        if r['structure'] == 'single_gene' and r['ambiguous_multiple_correspondence'] == 'True':
            for x, y in ((r['occurrence_a'], r['occurrence_b']), (r['occurrence_b'], r['occurrence_a'])):
                links[x][occs[y]['sample']].append(y)
    xs = [(x, sorted(ys)) for x, d in links.items() for ys in d.values() if len(ys) >= 2]
    if xs:
        x, ys = min(xs, key=lambda t: (-len(t[1]), t[0]))
        out.append(dict(case='multiple_correspondence', key=('single_gene', occs[x]['core_key']), rows=[x] + ys[:3]))
    return out


def draw(cases, case_rows, occs, cores, genes, unit_of, name):
    W_, y, box, gap = 1500, 70, 72, 4
    parts = []
    for c in cases:
        k = cores[c['key']]
        parts.append('<text x="14" y="%d" font-size="12" font-weight="bold">%s | %s %s: %s | %s occ / %s assemblies, copies %s</text>' % (
            y, esc(c['case']), esc(c['key'][0]), esc(c['key'][1]), esc(k['products'][:50]), k['occurrences'], k['assemblies'],
            k['copies_per_assembly']))
        y += 8
        for x in c['rows']:
            o = occs[x]
            parts.append('<text x="14" y="%d" font-size="10">%s %s:%s-%s (%s)</text>' % (
                y + 14, o['sample'], o['contig'], int(o['start0']) + 1, o['end0'], o['strand']))
            parts.append('<text x="14" y="%d" font-size="9" fill="#555">%s</text>' % (
                y + 25, esc(unit_of.get(x, 'unresolved: ' + o['unresolved_reason']).replace('SG|s|', 'group ').replace('single:', 'single-obs '))))
            up = list(zip(o['host_up'].split(','), o['host_up_genes'].split(',')))[::-1] if o['host_up'] else []
            down = list(zip(o['host_down'].split(','), o['host_down_genes'].split(','))) if o['host_down'] else []
            xpos = 300 + (5 - len(up)) * (box + gap)
            items = up + [('SKIP', o['skipped_mobile_up'])] + [('CORE', '')] + [('SKIP', o['skipped_mobile_down'])] + down
            for lab, uid in items:
                if lab == 'SKIP':
                    if uid != '0':
                        parts.append('<rect x="%d" y="%d" width="22" height="18" fill="none" stroke="#c0392b" stroke-dasharray="3,2"/>'
                                     '<text x="%d" y="%d" font-size="8" fill="#c0392b">%sm</text>' % (xpos, y + 6, xpos + 3, y + 18, uid))
                    xpos += 26
                    continue
                core = lab == 'CORE'
                rev = lab.endswith('(rev)')
                lab = lab.replace('(rev)', '')
                fill = '#edc06f' if core else colour(lab)
                pts = ('%d,%d %d,%d %d,%d %d,%d %d,%d' % (xpos, y + 6, xpos + box - 9, y + 6, xpos + box, y + 15, xpos + box - 9, y + 24, xpos, y + 24)
                       if not rev else '%d,%d %d,%d %d,%d %d,%d %d,%d' % (xpos + 9, y + 6, xpos + box, y + 6, xpos + box, y + 24, xpos + 9, y + 24, xpos, y + 15))
                prod = 'core' if core else genes[uid]['product']
                parts.append('<polygon points="%s" fill="%s" stroke="#666" stroke-width="0.5"><title>%s</title></polygon>' % (
                    pts, fill, esc('core' if core else '%s | %s' % (lab, prod))))
                parts.append('<text x="%d" y="%d" font-size="7">%s</text>' % (xpos + 10, y + 18, esc(prod[:12])))
                xpos += box + gap
            y += 32
        for r in [r for r in case_rows if r['case'] == c['case']]:
            parts.append('<text x="300" y="%d" font-size="9" fill="#333">%s vs %s: %s | host up %s (%s) down %s (%s) | %s | '
                         'outer flank DNA up %s (%s) down %s (%s)%s</text>' % (
                y + 10, r['occurrence_a'].rsplit('|', 1)[1], r['occurrence_b'].rsplit('|', 1)[1], r['verdict'], r['up_side'],
                r['up_support'], r['down_side'], r['down_support'], r['context_class'], r['dna_up'], r['dna_outer_cov_up'],
                r['dna_down'], r['dna_outer_cov_down'], ' | ambiguous link' if r['ambiguous_multiple_correspondence'] == 'True' else ''))
            y += 12
        y += 22
    head = ['<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d">' % (W_, y + 10),
            '<rect width="100%" height="100%" fill="white"/>',
            '<text x="14" y="24" font-size="16">%s: S6 v0.1 review cases (target scope, assembly only)</text>' % esc(name),
            '<text x="14" y="42" font-size="10">Host windows (5 genes per side after skipping mobile-keyword genes, red box = '
            'number skipped) in core orientation, coloured by homology group. Outer flank DNA = blastn coverage of '
            '%d-%d bp beyond each core end; a cross-check only.</text>' % (INNER, W),
            '<text x="14" y="56" font-size="10">Position units (site groups, single-observation sites) are not independent '
            'insertion events.</text>']
    return '\n'.join(head + parts + ['</svg>']) + '\n'


if __name__ == '__main__':
    main()
