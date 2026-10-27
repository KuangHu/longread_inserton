"""S6 v0.1 prototype: site groups and the stable-end list for the reviewed protein-entry cores of one subset.

A site group is repeated observation of one core architecture (structure + core key: single- vs two-gene and gene
composition are never mixed) at a corresponding genome position. Evidence stays assembly_only; a site group is not an
independent insertion event, and neither gaps nor empty sites are used or required. Rules: S6_SITE_GROUPS_v0.1.md.

Occurrences: every target_lt2500 occurrence of every core reviewed by core_review.py (>= 1 different-context pair).
Windows per side in core orientation: near = the first K=3 genes as they are; host = the first N=5 genes after
skipping genes whose product matches a mobile keyword (review_protein.hint), scanning at most MAX_SCAN genes. A
family being a candidate core is never a reason to skip. Labels are homology groups (CDS) or RNA type + product.

Side relation (host windows): labels occurring more than once in either window are ambiguous and never match;
matches must keep the same relative orientation, and the matched labels must lie in the same order in both windows
(longest increasing chain). same = >= MIN_SUPPORT ordered, oriented matches; different = both windows full (N genes,
all labelled) and no shared label; otherwise unknown (unknown is never different).

Pair verdict:
  different_position_same_assembly   same assembly, non-overlapping (both sides 'same' -> repeated_flank flag)
  same_position_candidate            different assemblies, both sides same
  different_position_candidate       different assemblies, both sides different
  conflicting_sides                  one side same, the other different
  unresolved                         anything else (missing labels, short context, shared but unordered labels)
An occurrence with same_position_candidate links to >= 2 occurrences of one other assembly keeps those links as
ambiguous_multiple_correspondence (not used for grouping).

Grouping (no transitive merging): connected components of the remaining same_position_candidate links are accepted
as a site group only if every pair inside is same_position_candidate and no assembly appears twice; otherwise every
member stays unresolved (conflict_in_component). An occurrence without same-position links whose every pair is a
different-position verdict is a single_observation_site. Everything else is unresolved with a reason; unresolved
occurrences are never counted as sites.

Stable-end list (from core_review outputs): coding span, shared-sequence estimate, per-side status, valid partners
per side, their assemblies and site groups, single-partner flag, strict < 2500 bp re-check on the estimated span.
Occurrences with an unknown side go to a separate boundary review list; nothing is dropped.
"""
import argparse
import json
from argparse import Namespace
from collections import Counter, defaultdict
from pathlib import Path

import core_review as cr
import protein_context as pc
from review_protein import colour, esc, hint

VERSION = 's6_site_groups_v0.1'
N = 5
MIN_SUPPORT = 2
MAX_SCAN = 15
MOBILE_HINTS = cr.MOBILE_HINTS


def lis(seq):
    """Length of the longest strictly increasing subsequence."""
    tails = []
    for x in seq:
        lo, hi = 0, len(tails)
        while lo < hi:
            mid = (lo + hi) // 2
            if tails[mid] < x:
                lo = mid + 1
            else:
                hi = mid
        tails[lo:lo + 1] = [x]
    return len(tails)


def side_relation(x, y, n=N, need=MIN_SUPPORT):
    """x, y: host windows [(label, same direction as core, gene_uid)], nearest first.
    Returns (same | different | unknown, ordered oriented matches, repeated label seen)."""
    lx, ly = Counter(l for l, _, _ in x if l), Counter(l for l, _, _ in y if l)
    rep = {l for l, c in lx.items() if c > 1} | {l for l, c in ly.items() if c > 1}
    py = {l: (i, d) for i, (l, d, _) in enumerate(y) if l and l not in rep}
    pts = [py[l][0] for l, d, _ in x if l and l not in rep and l in py and py[l][1] == d]
    support = lis(pts)
    if support >= need:
        return 'same', support, bool(rep)
    full = len(x) == n == len(y) and all(l for l, _, _ in x + y)
    if full and not set(lx) & set(ly):
        return 'different', support, bool(rep)
    return 'unknown', support, bool(rep)


def overlaps(a, b):
    return a['sample'] == b['sample'] and a['contig'] == b['contig'] and \
        int(a['start0']) < int(b['end0']) and int(b['start0']) < int(a['end0'])


def pair_verdict(a, b):
    up = side_relation(a['host']['up'], b['host']['up'])
    down = side_relation(a['host']['down'], b['host']['down'])
    sides = {up[0], down[0]}
    if a['sample'] == b['sample']:
        if overlaps(a, b):
            v = 'unresolved'
        else:
            v = 'different_position_same_assembly'
    elif sides == {'same'}:
        v = 'same_position_candidate'
    elif sides == {'different'}:
        v = 'different_position_candidate'
    elif sides == {'same', 'different'}:
        v = 'conflicting_sides'
    else:
        v = 'unresolved'
    return dict(verdict=v, up_side=up[0], up_support=up[1], down_side=down[0], down_support=down[1],
                repeated_label=up[2] or down[2], repeated_flank=a['sample'] == b['sample'] and sides == {'same'})


DIFFERENT = {'different_position_same_assembly', 'different_position_candidate'}


def group(occ_ids, sample_of, verdicts):
    """occ_ids of one core; verdicts {frozenset(a, b): verdict dict}. Returns ({occ: (status, group index or reason)},
    groups as sorted member lists, ambiguous link set)."""
    same = defaultdict(set)
    for k, v in verdicts.items():
        if v['verdict'] == 'same_position_candidate':
            a, b = tuple(k)
            same[a].add(b); same[b].add(a)
    ambiguous = set()
    for a in occ_ids:
        per = defaultdict(list)
        for b in same[a]:
            per[sample_of[b]].append(b)
        for xs in per.values():
            if len(xs) > 1:
                ambiguous |= {frozenset((a, b)) for b in xs}
    links = defaultdict(set)
    for a in occ_ids:
        for b in same[a]:
            if frozenset((a, b)) not in ambiguous:
                links[a].add(b)
    out, groups, seen = {}, [], set()
    for a in sorted(occ_ids):
        if a in seen or not links[a]:
            continue
        comp, stack = set(), [a]
        while stack:
            x = stack.pop()
            if x not in comp:
                comp.add(x); stack.extend(links[x] - comp)
        seen |= comp
        members = sorted(comp)
        clique = all(verdicts[frozenset((x, y))]['verdict'] == 'same_position_candidate' and
                     frozenset((x, y)) not in ambiguous for i, x in enumerate(members) for y in members[i + 1:])
        if clique and len({sample_of[x] for x in members}) == len(members):
            groups.append(members)
            for x in members:
                out[x] = ('site_group', len(groups) - 1)
        else:
            for x in members:
                out[x] = ('unresolved', 'conflict_in_component')
    for a in sorted(occ_ids):
        if a in out:
            continue
        rel = [verdicts[frozenset((a, b))]['verdict'] for b in occ_ids if b != a]
        if any(frozenset((a, b)) in ambiguous for b in occ_ids if b != a):
            out[a] = ('unresolved', 'ambiguous_multiple_correspondence')
        elif rel and all(r in DIFFERENT for r in rel):
            out[a] = ('single_observation_site', None)
        elif 'conflicting_sides' in rel:
            out[a] = ('unresolved', 'conflicting_sides')
        else:
            out[a] = ('unresolved', 'unresolved_pairs')
    return out, groups, ambiguous


def windows(o, gs, labels, mobile):
    idx = [int(u.rsplit('|', 1)[1]) for u in o['genes'].split(';')]
    fwd = o['strand'] == '+'
    nu, nd = cr.scan(gs, min(idx), max(idx), fwd, labels, lambda g: '', max_scan=pc.K)
    su, sd = cr.scan(gs, min(idx), max(idx), fwd, labels, mobile, max_scan=MAX_SCAN, keep=N)
    o['near'] = dict(up=[x[:3] for x in nu], down=[x[:3] for x in nd])
    o['scan'] = dict(up=su, down=sd)
    o['host'] = dict(up=[x[:3] for x in su if not x[3]][:N], down=[x[:3] for x in sd if not x[3]][:N])


def fmt(xs):
    return ','.join('%s%s' % (x[0] or 'nolabel', '' if x[1] else '(rev)') for x in xs)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--subset-dir', type=Path, required=True)
    p.add_argument('--work', type=Path, required=True, help='protein work dir (genes.tsv.gz, proteins.faa, mmseqs/)')
    p.add_argument('--mmseqs', default='/global/scratch/users/kh36969/longread_callout/env/s4/bin/mmseqs')
    p.add_argument('--threads', type=int, default=8)
    a = p.parse_args()
    pe = a.subset_dir / 'protein_entry'
    rv = pe / 'core_review'
    out = pe / 's6'
    out.mkdir(exist_ok=True)

    genes = pc.read_tsv(a.work / 'genes.tsv.gz')
    gene = {g['gene_uid']: g for g in genes}
    assign, _, _ = pc.assignments(Namespace(work=a.work, mmseqs=a.mmseqs, threads=a.threads))
    labels = {g['gene_uid']: pc.label(g, assign['homology_group']) for g in genes}
    by_contig = defaultdict(list)
    for g in genes:
        by_contig[(g['sample'], g['contig'])].append(g)
    for gs in by_contig.values():
        gs.sort(key=lambda g: int(g['index']))

    def mobile(g):
        return 'keyword' if hint(g['product']) in MOBILE_HINTS else ''

    cores = pc.read_tsv(rv / 'review_cores.tsv')
    keys = {(c['structure'], c['core_key']) for c in cores}
    occs = {o['occurrence_id']: o for o in pc.read_tsv(pe / 'protein_occurrences.tsv.gz')
            if o['in_target_length'] == 'True' and (o['structure'], o['core_key']) in keys}
    review = {r['occurrence_id']: r for r in pc.read_tsv(rv / 'review_occurrences.tsv.gz')}
    for o in occs.values():
        windows(o, by_contig[(o['sample'], o['contig'])], labels, mobile)

    by_core = defaultdict(list)
    for oid, o in sorted(occs.items()):
        by_core[(o['structure'], o['core_key'])].append(oid)
    pair_rows, group_rows, evidence_rows, core_rows, site = [], [], [], [], {}
    for key in sorted(by_core):
        ids = by_core[key]
        verdicts = {}
        for i, x in enumerate(ids):
            for y in ids[i + 1:]:
                verdicts[frozenset((x, y))] = v = pair_verdict(occs[x], occs[y])
                pair_rows.append(dict(structure=key[0], core_key=key[1], occurrence_a=x, occurrence_b=y,
                                      sample_a=occs[x]['sample'], sample_b=occs[y]['sample'], **v))
        res, groups, amb = group(ids, {x: occs[x]['sample'] for x in ids}, verdicts)
        for r in pair_rows[len(pair_rows) - len(verdicts):]:
            r['ambiguous_multiple_correspondence'] = frozenset((r['occurrence_a'], r['occurrence_b'])) in amb
        gid = ['SG|%s|%s|%d' % (key[0][0], key[1], i) for i in range(len(groups))]
        for i, members in enumerate(groups):
            group_rows.append(dict(site_group=gid[i], structure=key[0], core_key=key[1], members=len(members),
                assemblies=len({occs[x]['sample'] for x in members}), occurrences=';'.join(members),
                locations=';'.join('%s:%s:%s-%s(%s)' % (occs[x]['sample'], occs[x]['contig'], int(occs[x]['start0']) + 1,
                                                       occs[x]['end0'], occs[x]['strand']) for x in members)))
            for j, x in enumerate(members):
                for y in members[j + 1:]:
                    v = verdicts[frozenset((x, y))]
                    evidence_rows.append(dict(site_group=gid[i], occurrence_a=x, occurrence_b=y, up_support=v['up_support'],
                                              down_support=v['down_support'], repeated_label=v['repeated_label']))
        for x in ids:
            st, why = res[x]
            site[x] = dict(status=st, site_group=gid[why] if st == 'site_group' else '',
                           reason=why if st == 'unresolved' else '')
        sts = Counter(site[x]['status'] for x in ids)
        supported = defaultdict(int)
        for x in ids:
            if site[x]['status'] == 'single_observation_site':
                supported[occs[x]['sample']] += 1
        for members in groups:
            for x in members:
                supported[occs[x]['sample']] += 1
        c = next(c for c in cores if (c['structure'], c['core_key']) == key)
        core_rows.append(dict(structure=key[0], core_key=key[1], products=c['products'], occurrences=len(ids),
            assemblies=len({occs[x]['sample'] for x in ids}), copies_per_assembly=c['copies_per_assembly'],
            copy_pattern=c['copy_pattern'], multi_observation_site_groups=len(groups),
            occurrences_in_site_groups=sts['site_group'], single_observation_sites=sts['single_observation_site'],
            unresolved=sts['unresolved'],
            unresolved_reasons=','.join('%s=%d' % kv for kv in sorted(Counter(site[x]['reason'] for x in ids
                                                                            if site[x]['status'] == 'unresolved').items())),
            assemblies_with_2plus_supported_positions=sum(n >= 2 for n in supported.values()),
            pair_verdicts=','.join('%s=%d' % kv for kv in sorted(Counter(v['verdict'] for v in verdicts.values()).items())),
            product_hint=c['product_hint']))

    occ_rows = []
    for oid, o in sorted(occs.items()):
        r = review.get(oid, {})
        occ_rows.append(dict(occurrence_id=oid, structure=o['structure'], core_key=o['core_key'], sample=o['sample'],
            contig=o['contig'], start0=o['start0'], end0=o['end0'], strand=o['strand'], coding_span=o['span'],
            s6_status=site[oid]['status'], site_group=site[oid]['site_group'], unresolved_reason=site[oid]['reason'],
            boundary_up=r.get('up_status', 'not_reviewed_no_different_partner'),
            boundary_down=r.get('down_status', 'not_reviewed_no_different_partner'),
            estimated_span=r.get('estimated_span', ''),
            near_up=fmt(o['near']['up']), near_down=fmt(o['near']['down']),
            near_mobile_keyword=sum(bool(mobile(gene[x[2]])) for s in ('up', 'down') for x in o['near'][s]),
            host_up=fmt(o['host']['up']), host_down=fmt(o['host']['down']),
            host_up_genes=','.join(x[2] for x in o['host']['up']), host_down_genes=','.join(x[2] for x in o['host']['down']),
            skipped_mobile_up=sum(bool(x[3]) for x in o['scan']['up']), skipped_mobile_down=sum(bool(x[3]) for x in o['scan']['down'])))
    pc.write_tsv(out / 's6_occurrences.tsv.gz', occ_rows, list(occ_rows[0]))
    pc.write_tsv(out / 's6_pairs.tsv.gz', pair_rows, list(pair_rows[0]))
    pc.write_tsv(out / 's6_site_groups.tsv', group_rows, ['site_group', 'structure', 'core_key', 'members', 'assemblies',
                                                          'occurrences', 'locations'])
    pc.write_tsv(out / 's6_group_evidence.tsv', evidence_rows, ['site_group', 'occurrence_a', 'occurrence_b', 'up_support',
                                                                'down_support', 'repeated_label'])
    pc.write_tsv(out / 's6_cores.tsv', core_rows, list(core_rows[0]))

    # stable-end list and boundary review list (occurrences with >= 1 different-context partner)
    partners = defaultdict(list)
    for r in pc.read_tsv(rv / 'review_pairs.tsv.gz'):
        for me, other, pre in ((r['occurrence_a'], r['occurrence_b'], 'a'), (r['occurrence_b'], r['occurrence_a'], 'b')):
            partners[me].append((other, r[pre + '_up'] != '', r[pre + '_down'] != ''))
    stable, unknown = [], []
    for oid, r in sorted(review.items()):
        ps = partners[oid]
        valid = [p for p in ps if p[1] and p[2]]

        def tag(x):
            return site[x]['site_group'] or ('single:' + x if site[x]['status'] == 'single_observation_site' else 'unresolved')
        est = r['estimated_span']
        row = dict(occurrence_id=oid, structure=r['structure'], core_key=r['core_key'], sample=r['sample'], contig=r['contig'],
            start0=r['start0'], end0=r['end0'], strand=r['strand'], coding_span=r['coding_span'],
            ext_up=r['ext_up'], up_status=r['up_status'], ext_down=r['ext_down'], down_status=r['down_status'],
            shared_estimate_span=est, strictly_lt2500='' if est == '' else int(est) < pc.TARGET_MAX,
            different_partners=len(ps), valid_partners_up=sum(p[1] for p in ps), valid_partners_down=sum(p[2] for p in ps),
            valid_partners_both=len(valid), valid_partner_assemblies=len({occs[p[0]]['sample'] for p in valid}),
            valid_partner_site_groups=len({tag(p[0]) for p in valid if tag(p[0]) != 'unresolved'}),
            valid_partners_unresolved_site=sum(tag(p[0]) == 'unresolved' for p in valid),
            single_partner=len(valid) == 1, s6_status=site[oid]['status'], site_group=site[oid]['site_group'],
            label='candidate_conserved_end' if r['boundary_resolved'] == 'True' else 'boundary_unknown')
        (stable if r['boundary_resolved'] == 'True' else unknown).append(row)
    cols = list((stable or unknown)[0])
    pc.write_tsv(out / 'stable_end_candidates.tsv', stable, cols)
    pc.write_tsv(out / 'boundary_review_unknown.tsv', unknown, cols)

    summary = dict(version=VERSION, subset=a.subset_dir.name, evidence_source='assembly_only', scope='target_lt2500',
        parameters=dict(N_host=N, K_near=pc.K, min_support=MIN_SUPPORT, max_scan=MAX_SCAN, mobile_skip='product keyword only',
                        mobile_hints=sorted(MOBILE_HINTS)),
        site_groups_are_independent_insertion_events=False, production_ready=False, structures={})
    for st in ('single_gene', 'two_gene'):
        cs = [c for c in core_rows if c['structure'] == st]
        os_ = [r for r in occ_rows if r['structure'] == st]
        ps = [r for r in pair_rows if r['structure'] == st]
        gs = [g for g in group_rows if g['structure'] == st]
        se, uk = [r for r in stable if r['structure'] == st], [r for r in unknown if r['structure'] == st]
        summary['structures'][st] = dict(
            cores=len(cs), occurrences=len(os_), assemblies=len({r['sample'] for r in os_}),
            occurrence_status=dict(Counter(r['s6_status'] for r in os_)),
            unresolved_reasons=dict(Counter(r['unresolved_reason'] for r in os_ if r['s6_status'] == 'unresolved')),
            multi_observation_site_groups=len(gs), site_group_size=dict(Counter(str(g['members']) for g in gs)),
            single_observation_sites=sum(r['s6_status'] == 'single_observation_site' for r in os_),
            cores_with_site_group=sum(c['multi_observation_site_groups'] > 0 for c in cs),
            cores_with_assembly_2plus_supported_positions=sum(c['assemblies_with_2plus_supported_positions'] > 0 for c in cs),
            pair_verdicts=dict(Counter(r['verdict'] for r in ps)),
            pairs_ambiguous_multiple_correspondence=sum(r['ambiguous_multiple_correspondence'] for r in ps),
            pairs_repeated_flank_same_assembly=sum(r['repeated_flank'] for r in ps),
            stable_end_candidates=len(se), stable_end_strictly_lt2500=sum(r['strictly_lt2500'] is True for r in se),
            stable_end_single_partner=sum(r['single_partner'] for r in se),
            stable_end_in_site_group=sum(r['s6_status'] == 'site_group' for r in se),
            boundary_unknown=len(uk))
    cases = pick(core_rows)
    summary['cases'] = [c['core_key'] for c in cases]
    summary['case_rule'] = pick.__doc__
    (out / 's6_summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    (out / 'cases_s6.svg').write_text(draw(cases, occs, site, occ_rows, gene, a.subset_dir.name))
    print(json.dumps(summary, indent=2))


def pick(core_rows):
    """Fixed rule: per structure, the core with the most multi-observation site groups (ties: occurrences in groups,
    then key); plus the single-gene core with the most assemblies holding >= 2 supported positions, if different."""
    out = []
    for st in ('single_gene', 'two_gene'):
        xs = [c for c in core_rows if c['structure'] == st and c['multi_observation_site_groups'] > 0]
        if xs:
            out.append(min(xs, key=lambda c: (-c['multi_observation_site_groups'], -c['occurrences_in_site_groups'], c['core_key'])))
    xs = [c for c in core_rows if c['structure'] == 'single_gene' and c['assemblies_with_2plus_supported_positions'] > 0
          and all(c is not o for o in out)]
    if xs:
        out.append(min(xs, key=lambda c: (-c['assemblies_with_2plus_supported_positions'], -c['occurrences'], c['core_key'])))
    return out


def draw(cases, occs, site, occ_rows, gene, name):
    """Per case: up to 3 site groups x 2 members, up to 2 single-observation sites, up to 2 unresolved occurrences."""
    W_, y, box, gap = 1500, 76, 70, 4
    parts = []
    for c in cases:
        key = (c['structure'], c['core_key'])
        ids = sorted(x for x, o in occs.items() if (o['structure'], o['core_key']) == key)
        rows = []
        for g in sorted({site[x]['site_group'] for x in ids if site[x]['site_group']})[:3]:
            rows += [(x, 'site group ' + g.rsplit('|', 1)[1]) for x in ids if site[x]['site_group'] == g][:2]
        rows += [(x, 'single-observation site') for x in ids if site[x]['status'] == 'single_observation_site'][:2]
        rows += [(x, 'unresolved: ' + site[x]['reason']) for x in ids if site[x]['status'] == 'unresolved'][:2]
        parts.append('<text x="14" y="%d" font-size="11" font-weight="bold">%s %s: %s | %s occ / %s assemblies, copies %s | '
                     'site groups %s (%s occ), single-observation %s, unresolved %s (%s)</text>' % (
                         y, esc(c['structure']), esc(c['core_key']), esc(c['products'][:45]), c['occurrences'], c['assemblies'],
                         c['copies_per_assembly'], c['multi_observation_site_groups'], c['occurrences_in_site_groups'],
                         c['single_observation_sites'], c['unresolved'], esc(c['unresolved_reasons'])))
        y += 10
        for x, tag in rows:
            o = occs[x]
            parts.append('<text x="14" y="%d" font-size="10">%s %s:%s-%s (%s)</text>' % (
                y + 14, o['sample'], o['contig'], int(o['start0']) + 1, o['end0'], o['strand']))
            parts.append('<text x="14" y="%d" font-size="9" fill="#555">%s</text>' % (y + 25, esc(tag)))
            up, down = o['scan']['up'][:7][::-1], o['scan']['down'][:7]
            xpos = 300 + (7 - len(up)) * (box + gap)
            for item in up + ['CORE'] + down:
                core = item == 'CORE'
                lab, same, uid, mob = ('core', True, '', '') if core else item
                fill = '#edc06f' if core else colour(lab or 'nolabel')
                pts = ('%d,%d %d,%d %d,%d %d,%d %d,%d' % (xpos, y + 6, xpos + box - 9, y + 6, xpos + box, y + 15, xpos + box - 9, y + 24, xpos, y + 24)
                       if same else '%d,%d %d,%d %d,%d %d,%d %d,%d' % (xpos + 9, y + 6, xpos + box, y + 6, xpos + box, y + 24, xpos + 9, y + 24, xpos, y + 15))
                stroke = 'stroke="#c0392b" stroke-width="2" stroke-dasharray="3,2"' if mob else 'stroke="#666" stroke-width="0.5"'
                tip = 'core' if core else '%s | %s%s' % (lab or 'nolabel', gene[uid]['product'], ' | skipped (mobile keyword)' if mob else '')
                parts.append('<polygon points="%s" fill="%s" %s><title>%s</title></polygon>' % (pts, fill, stroke, esc(tip)))
                parts.append('<text x="%d" y="%d" font-size="7">%s</text>' % (xpos + 10, y + 18, esc('core' if core else gene[uid]['product'][:12])))
                xpos += box + gap
            y += 32
        y += 14
    head = ['<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d">' % (W_, y + 10),
            '<rect width="100%" height="100%" fill="white"/>',
            '<text x="14" y="24" font-size="16">%s: S6 v0.1 site groups (target scope, assembly only)</text>' % esc(name),
            '<text x="14" y="42" font-size="10">Neighbour genes in core orientation, coloured by homology group (same colour = '
            'homologous); red dashed = mobile-keyword gene skipped for the host window. Site groups are repeated observations '
            'of one position, not independent insertion events.</text>',
            '<text x="14" y="56" font-size="10">Same position needs >= 2 ordered, orientation-consistent unique host matches on '
            'each side; groups are cliques of such pairs with one occurrence per assembly. Hover a gene for its product.</text>']
    return '\n'.join(head + parts + ['</svg>']) + '\n'


if __name__ == '__main__':
    main()
