"""Candidate review at scale (after scale_entry compare): S6 site groups with physical-overlap flags. Evidence stays
assembly_only; production_ready=false; site groups are not independent insertion events.

windows   per batch, from the genes tables: for every target-scope occurrence (single- and two-gene) the near window
          (K=3 genes as they are) and the S6 host window (N=5 genes per side after skipping mobile-keyword genes,
          at most 15 scanned), in core orientation with order, orientation and gene ids (site_groups.windows,
          unchanged). The first-pass K=3 label sets are not used as S6 evidence.
s6        per compare shard. Reviewed cores: every target-scope core with at least one pair that is not
          similar_context (different_context_candidate, insufficient_context_labels, one_side_shared_unresolved or
          edge_or_missing_context), so occurrences the first pass could not label still get host review. Per core,
          over all its occurrences: site_groups.pair_verdict on every pair and site_groups.group (clique-only, one
          occurrence per assembly, ambiguity kept); s6_review.unit_relation between position units and the largest
          pairwise-different unit sets (s6_review.max_clique). Every same-assembly pair gets the context class
          (s6_review.context_class), the copy gap on the contig and the number of host genes shared physically.
          Two separate dimensions per pair:
            local structure  U-core-D structural pair: one assembly, host context different on both sides, no host
                             gene shared physically;
            S6 grouping      supported_positions (both copies in a site group or a single-observation site) or
                             grouping_unresolved. It answers which position the copy is in other assemblies and whether
                             it can be counted once; it never vetoes the local structure. Grouping-unresolved copies are
                             not used for site or event counts.
          The strict rule v0 (structure and supported positions) is kept as ucd_pair_candidate_strict_v0.
          Cross-assembly single-copy pairs: two assemblies holding one copy each of the core, host context different
          on both sides (different_position_candidate), with the same S6 grouping field.
boundary  per shard, on the structural pairs (both S6 grouping tiers) and the cross-assembly single-copy pairs: core
          protein identity (blastp) and core DNA identity / coverage, and the coding span +/- W bp of both copies
          aligned (core_review.blast / chain / extensions, unchanged), one all-vs-all per core: a side is measured only
          if the core end is covered in both copies, unknown otherwise (never 0). Per occurrence and pair scope,
          core_review.side_boundary over all its partners in that scope; estimated span only when both sides are
          resolved, then classed < 2500 or >= 2500 bp; unknown boundaries are kept as unknown.
summarize per structure: reviewed cores, occurrence status, unit relations, same-assembly context classes with the
          physical flags, U-core-D pair candidates, the frozen regression set.
"""
import argparse
import gzip
import json
from collections import Counter, defaultdict
from pathlib import Path

import core_review as cr
import protein_context as pc
import s6_review as rv
import scale_entry as se
import site_groups as sg
from review_protein import hint

VERSION = 'scale_review_v0.1'
WIN_COLS = ['occurrence_id', 'near', 'host', 'skipped_mobile_up', 'skipped_mobile_down', 'near_mobile_keyword']
SUPPORTED = {'site_group', 'single_observation_site'}


def mobile(g):
    return 'keyword' if hint(g['product']) in sg.MOBILE_HINTS else ''


def cmd_windows(a):
    cat = se.load_catalog(a.work)
    hg = {h: r['homology_group'] for h, r in cat.items()}
    del cat
    out = a.work / 'review' / 'windows' / ('b%d' % a.batch)
    out.mkdir(parents=True, exist_ok=False)
    fh = [gzip.open(out / ('shard_%03d.tsv.gz' % k), 'wt') for k in range(se.SHARDS)]
    for f in fh:
        f.write('\t'.join(WIN_COLS) + '\n')
    n = 0
    want = defaultdict(list)
    for k in range(se.SHARDS):
        for o in pc.read_tsv(a.work / 'occ' / ('b%d' % a.batch) / ('shard_%03d.tsv.gz' % k)):
            if o['in_target_length'] == 'True':
                want[o['sample']].append((k, o['occurrence_id'], o['genes'], o['strand'], o['contig']))
    for smp, os_ in sorted(want.items()):
        genes = pc.read_tsv(a.work / 'genes' / (smp + '.tsv.gz'))
        labels = {g['gene_uid']: pc.label(g, hg) for g in genes}
        by_contig = defaultdict(list)
        for g in genes:
            by_contig[g['contig']].append(g)
        for gs in by_contig.values():
            gs.sort(key=lambda g: int(g['index']))
        for k, oid, gids, strand, contig in os_:
            o = dict(genes=gids, strand=strand)
            sg.windows(o, by_contig[contig], labels, mobile)
            fh[k].write('\t'.join([oid, json.dumps(o['near']), json.dumps(o['host']),
                                   str(sum(bool(x[3]) for x in o['scan']['up'])), str(sum(bool(x[3]) for x in o['scan']['down'])),
                                   str(sum(bool(mobile(by_contig[contig][int(x[2].rsplit('|', 1)[1])]))
                                           for s in ('up', 'down') for x in o['near'][s]))]) + '\n')
            n += 1
    for f in fh:
        f.close()
    (out / 'windows.json').write_text(json.dumps(dict(batch=a.batch, occurrences=n, assemblies=len(want)), indent=1) + '\n')


def side_list(js):
    return {s: [tuple(x) for x in v] for s, v in json.loads(js).items()}


def review_core(ys):
    """ys: occurrence dicts of one core with host windows. Returns occurrence rows, groups, same-assembly pair rows,
    unit-relation rows and a core summary."""
    ids = [o['occurrence_id'] for o in ys]
    occ = {o['occurrence_id']: o for o in ys}
    verdicts = {}
    for i, x in enumerate(ids):
        for y in ids[i + 1:]:
            verdicts[frozenset((x, y))] = sg.pair_verdict(occ[x], occ[y])
    res, groups, amb = sg.group(ids, {x: occ[x]['sample'] for x in ids}, verdicts)
    unit = {}
    for gi, members in enumerate(groups):
        for x in members:
            unit[x] = 'G%d' % gi
    for x, (st, _) in res.items():
        if st == 'single_observation_site':
            unit[x] = 'S:' + x
    members = defaultdict(list)
    for x, u in unit.items():
        members[u].append(x)
    units = sorted(members)
    adj_d, adj_c = defaultdict(set), defaultdict(set)
    rels, rel_rows = Counter(), []
    for i, u in enumerate(units):
        for v in units[i + 1:]:
            cross = [dict(verdicts[frozenset((x, y))], sample_a=occ[x]['sample'], sample_b=occ[y]['sample'])
                     for x in members[u] for y in members[v]]
            rel, ctx = rv.unit_relation(cross)
            rels[(rel, ctx)] += 1
            if rel == 'different':
                adj_d[u].add(v); adj_d[v].add(u)
                if ctx == 'context_different_both_sides':
                    adj_c[u].add(v); adj_c[v].add(u)
            if u.startswith('G') and v.startswith('G'):
                rel_rows.append(dict(unit_a=u, unit_b=v, relation=rel, context=ctx, cross_pairs=len(cross)))
    sa_rows, cross_rows = [], []
    hostg = lambda o: {x[2] for s in ('up', 'down') for x in o['host'][s]}
    status = lambda x: res[x][0] if res[x][0] != 'unresolved' else 'unresolved:' + res[x][1]
    copies = Counter(o['sample'] for o in ys)
    for k, v in sorted(verdicts.items(), key=lambda kv: sorted(kv[0])):
        x, y = sorted(k)
        X, Y = occ[x], occ[y]
        both = res[x][0] in SUPPORTED and res[y][0] in SUPPORTED
        pos = dict(s6_status_a=status(x), s6_status_b=status(y), both_supported_positions=both,
                   s6_grouping='supported_positions' if both else 'grouping_unresolved')
        if v['verdict'] == 'different_position_same_assembly':
            gap = (max(int(X['start0']), int(Y['start0'])) - min(int(X['end0']), int(Y['end0']))) if X['contig'] == Y['contig'] else ''
            phys = len(hostg(X) & hostg(Y))
            ctx = rv.context_class(v['up_side'], v['down_side'])
            # local structure of the pair only; S6 grouping is a separate field and never vetoes it
            structural = ctx == 'context_different_both_sides' and phys == 0
            sa_rows.append(dict(occurrence_a=x, occurrence_b=y, sample=X['sample'], context_class=ctx,
                                up_support=v['up_support'], down_support=v['down_support'], copy_gap_bp=gap,
                                host_genes_shared_physically=phys, ucd_structural_pair=structural, **pos,
                                ucd_pair_candidate_strict_v0=structural and both))
        elif v['verdict'] == 'different_position_candidate' and copies[X['sample']] == copies[Y['sample']] == 1:
            # one copy in each assembly: host context different on both sides (no same-assembly pair can carry it)
            cross_rows.append(dict(occurrence_a=x, occurrence_b=y, sample_a=X['sample'], sample_b=Y['sample'],
                                   up_support=v['up_support'], down_support=v['down_support'], **pos))
    occ_rows = []
    for x in ids:
        st, why = res[x]
        o = occ[x]
        occ_rows.append(dict(occurrence_id=x, sample=o['sample'], contig=o['contig'], start0=o['start0'], end0=o['end0'],
                             strand=o['strand'], s6_status=st, unit=unit.get(x, ''),
                             unresolved_reason=why if st == 'unresolved' else '',
                             host_up=sg.fmt(o['host']['up']), host_down=sg.fmt(o['host']['down']),
                             near_mobile_keyword=o['near_mobile_keyword'],
                             skipped_mobile=int(o['skipped_mobile_up']) + int(o['skipped_mobile_down'])))
    sts = Counter(r['s6_status'] for r in occ_rows)
    cls = Counter(r['context_class'] for r in sa_rows)
    summary = dict(occurrences=len(ids), assemblies=len({o['sample'] for o in ys}), site_groups=len(groups),
                   single_observation_sites=sts['single_observation_site'], unresolved=sts['unresolved'],
                   unresolved_reasons=','.join('%s=%d' % kv for kv in sorted(Counter(r['unresolved_reason'] for r in occ_rows
                                                                                     if r['unresolved_reason']).items())),
                   pair_verdicts=','.join('%s=%d' % kv for kv in sorted(Counter(v['verdict'] for v in verdicts.values()).items())),
                   ambiguous_links=len(amb), position_units=len(units),
                   unit_pairs_different=sum(n for (r, _), n in rels.items() if r == 'different'),
                   unit_pairs_unknown=rels[('unknown', '')], unit_pairs_conflicting=rels[('conflicting', '')],
                   unit_pairs_ambiguous=rels[('ambiguous', '')],
                   max_pairwise_different_units=rv.max_clique(units, adj_d) if units else 0,
                   max_pairwise_context_different_units=rv.max_clique(units, adj_c) if units else 0,
                   same_assembly_pairs=len(sa_rows),
                   **{'sa_' + c: cls[c] for c in ('context_different_both_sides', 'context_repeated_both_sides',
                                                  'context_repeated_one_side', 'context_unknown')},
                   sa_physically_shared_host_genes=sum(r['host_genes_shared_physically'] > 0 for r in sa_rows),
                   ucd_structural_pairs=sum(r['ucd_structural_pair'] for r in sa_rows),
                   ucd_structural_pairs_s6_supported=sum(r['ucd_pair_candidate_strict_v0'] for r in sa_rows),
                   ucd_assemblies=len({r['sample'] for r in sa_rows if r['ucd_structural_pair']}),
                   cross_different_position_pairs=sum(v['verdict'] == 'different_position_candidate' for v in verdicts.values()),
                   cross_single_copy_pairs=len(cross_rows),
                   cross_single_copy_pairs_s6_supported=sum(r['both_supported_positions'] for r in cross_rows))
    return occ_rows, groups, sa_rows, cross_rows, rel_rows, summary


def cmd_s6(a):
    batches = [int(x) for x in a.batches.split(',')]
    tag = 'b' + '-'.join(map(str, batches))
    cores = [c for c in pc.read_tsv(a.work / 'compare' / tag / ('shard_%03d' % a.shard) / 'cores.tsv.gz')
             if c['scope'] == 'target_lt2500' and int(c['similar_context']) < int(c['pairs_total'])]
    keys = {(c['structure'], c['core_key']): c for c in cores}
    win = {}
    for b in batches:
        for r in pc.read_tsv(a.work / 'review' / 'windows' / ('b%d' % b) / ('shard_%03d.tsv.gz' % a.shard)):
            win[r['occurrence_id']] = r
    groups_of = defaultdict(list)
    for o in se.read_shard(a.work, batches, a.shard):
        if o['in_target_length'] == 'True' and (o['structure'], o['core_key']) in keys:
            w = win[o['occurrence_id']]
            o.update(near=side_list(w['near']), host=side_list(w['host']), skipped_mobile_up=w['skipped_mobile_up'],
                     skipped_mobile_down=w['skipped_mobile_down'], near_mobile_keyword=w['near_mobile_keyword'])
            del o['up'], o['down']
            groups_of[(o['structure'], o['core_key'])].append(o)
    out = a.work / 'review' / 's6' / tag / ('shard_%03d' % a.shard)
    out.mkdir(parents=True, exist_ok=True)
    O, G, P, X, R, C = [], [], [], [], [], []
    for key in sorted(groups_of):
        ys = sorted(groups_of[key], key=lambda o: o['occurrence_id'])
        occ_rows, groups, sa_rows, cross_rows, rel_rows, summary = review_core(ys)
        base = dict(structure=key[0], core_key=key[1])
        O += [dict(base, **r) for r in occ_rows]
        for gi, members in enumerate(groups):
            G.append(dict(base, unit='G%d' % gi, members=len(members), occurrences=';'.join(members)))
        P += [dict(base, **r) for r in sa_rows]
        X += [dict(base, **r) for r in cross_rows]
        R += [dict(base, **r) for r in rel_rows]
        c = keys[key]
        C.append(dict(base, products=c['products'], product_hint=hint(c['products'].split(' | ')[0]),
                      k3_different_pairs=c['different_context_candidate'],
                      k3_different_pairs_same_assembly=c['different_context_candidate_same_assembly'],
                      k3_insufficient_pairs=c['insufficient_context_labels'], k3_edge_pairs=c['edge_or_missing_context'],
                      max_copies_in_one_assembly=c['max_copies_in_one_assembly'], **summary))
    for name, rows in (('occurrences', O), ('site_groups', G), ('same_assembly_pairs', P),
                       ('cross_assembly_single_copy_pairs', X), ('group_relations', R), ('cores', C)):
        pc.write_tsv(out / (name + '.tsv.gz'), rows, list(rows[0]) if rows else ['structure'])
    (out / 'run.json').write_text(json.dumps(dict(shard=a.shard, batches=batches, reviewed_cores=len(C),
                                                  occurrences=len(O)), indent=1) + '\n')


def structural_pairs(d):
    """Pairs that enter the structural review, both S6 grouping tiers: same-assembly U-core-D structural pairs and
    cross-assembly single-copy pairs."""
    same = [dict(r, pair_scope='same_assembly', sample_a=r['sample'], sample_b=r['sample'])
            for r in pc.read_tsv(d / 'same_assembly_pairs.tsv.gz') if r['ucd_structural_pair'] == 'True']
    cross = [dict(r, pair_scope='cross_assembly_single_copy')
             for r in pc.read_tsv(d / 'cross_assembly_single_copy_pairs.tsv.gz')]
    return same + cross


def span_class(est):
    if est is None:
        return 'unknown'
    return 'lt2500' if est < pc.TARGET_MAX else 'ge2500'


def cmd_boundary(a):
    batches = [int(x) for x in a.batches.split(',')]
    tag = 'b' + '-'.join(map(str, batches))
    pairs = structural_pairs(a.work / 'review' / 's6' / tag / ('shard_%03d' % a.shard))
    need = {x for r in pairs for x in (r['occurrence_a'], r['occurrence_b'])}
    occ = {o['occurrence_id']: o for o in se.read_shard(a.work, batches, a.shard) if o['occurrence_id'] in need}
    fasta = {m['sample']: m['fasta'] for m in pc.read_tsv(a.manifest)}
    by_sample = defaultdict(list)
    for x in need:
        by_sample[occ[x]['sample']].append(x)
    for smp in sorted(by_sample):
        seqs = pc.read_fna(fasta[smp])
        genes = {g['gene_uid']: g for g in pc.read_tsv(a.work / 'genes' / (smp + '.tsv.gz'))}
        prot = se.read_faa_gz(a.work / 'genes' / (smp + '.faa.gz'))
        for x in by_sample[smp]:
            o = occ[x]
            g = genes[o['genes'].split(';')[0]]
            o['ext'], o['c0'], o['c1'] = cr.window(seqs[o['contig']], int(o['start0']), int(o['end0']), o['strand'] == '+',
                                                   cr.W, g['circular_contig'] == 'True')
            o['core_dna'] = o['ext'][o['c0']:o['c1']]
            o['prots'] = [prot[genes[u]['protein_sha256']] for u in cr.core_genes(o)]
        del seqs
    rows_of = defaultdict(list)
    for r in pairs:
        rows_of[(r['structure'], r['core_key'])].append(r)
    P, partners = [], defaultdict(list)
    for key in sorted(rows_of):
        # one all-vs-all per core; blast runs with -subject, so each pair is scored on its own
        ids = sorted({x for r in rows_of[key] for x in (r['occurrence_a'], r['occurrence_b'])})
        num = {x: str(i) for i, x in enumerate(ids)}
        npos = 1 if key[0] == 'single_gene' else 2
        php = [cr.blast('blastp', {num[x] + '_%d' % i: occ[x]['prots'][i] for x in ids}, None) for i in range(npos)]
        dcore = cr.blast('blastn', {num[x]: occ[x]['core_dna'] for x in ids}, 'blastn')
        dext = cr.blast('blastn', {num[x]: occ[x]['ext'] for x in ids}, 'blastn')
        for r in rows_of[key]:
            A, B = occ[r['occurrence_a']], occ[r['occurrence_b']]
            ia, ib = num[A['occurrence_id']], num[B['occurrence_id']]
            pb = [cr.best(php[i][(ia + '_%d' % i, ib + '_%d' % i)]) for i in range(npos)]
            db = cr.best(dcore[(ia, ib)])
            blk = cr.chain([(h['qs'], h['qe'], h['ss'], h['se']) for h in dext[(ia, ib)]], A['c0'], A['c1'], B['c0'], B['c1'])
            ext = cr.extensions(blk, A['c0'], A['c1'], len(A['ext']), B['c0'], B['c1'], len(B['ext'])) if blk else \
                dict(up_cover='no_core_alignment', down_cover='no_core_alignment', a_up=None, b_up=None, a_down=None,
                     b_down=None, up_reaches_end=False, down_reaches_end=False)
            P.append(dict(structure=key[0], core_key=key[1], pair_scope=r['pair_scope'], s6_grouping=r['s6_grouping'],
                sample_a=r['sample_a'], sample_b=r['sample_b'], occurrence_a=A['occurrence_id'],
                occurrence_b=B['occurrence_id'], identical_protein=A['prots'] == B['prots'],
                protein_identity_min=min(x[0] for x in pb) if all(pb) else '',
                protein_coverage_min=round(min(min(x[1], x[2]) for x in pb), 4) if all(pb) else '',
                identical_core_dna=A['core_dna'].upper() == B['core_dna'].upper(),
                core_dna_identity=db[0] if db else '', core_dna_coverage_min=round(min(db[1], db[2]), 4) if db else '',
                up_cover=ext['up_cover'], a_up=ext['a_up'], b_up=ext['b_up'],
                down_cover=ext['down_cover'], a_down=ext['a_down'], b_down=ext['b_down']))
            for o, pre in ((A, 'a'), (B, 'b')):
                partners[(o['occurrence_id'], r['pair_scope'])].append(
                    (ext[pre + '_up'], ext['up_reaches_end'], ext[pre + '_down'], ext['down_reaches_end'],
                     r['s6_grouping'] == 'supported_positions'))
    O = []
    for x, scope in sorted(partners):
        o, ps = occ[x], partners[(x, scope)]
        up = cr.side_boundary([p[0] for p in ps], any(p[1] for p in ps))
        down = cr.side_boundary([p[2] for p in ps], any(p[3] for p in ps))
        span = int(o['span'])
        res = up[1] == down[1] == 'resolved'
        est = span + up[0] + down[0] if res else None
        O.append(dict(structure=o['structure'], core_key=o['core_key'], occurrence_id=x, sample=o['sample'],
                      pair_scope=scope, coding_span=span, partners=len(ps),
                      partners_s6_supported=sum(p[4] for p in ps), up_status=up[1],
                      ext_up='' if up[0] is None else up[0], up_partners_end_not_covered=up[2], down_status=down[1],
                      ext_down='' if down[0] is None else down[0], down_partners_end_not_covered=down[2],
                      boundary_resolved=res, estimated_span='' if est is None else est, span_class=span_class(est)))
    out = a.work / 'review' / 'boundary' / tag / ('shard_%03d' % a.shard)
    out.mkdir(parents=True, exist_ok=True)
    pc.write_tsv(out / 'pairs.tsv.gz', P, list(P[0]) if P else ['structure'])
    pc.write_tsv(out / 'occurrences.tsv.gz', O, list(O[0]) if O else ['structure'])
    (out / 'run.json').write_text(json.dumps(dict(shard=a.shard, pairs=len(P), occurrences=len(O)), indent=1) + '\n')


def boundary_stats(bp, bo):
    num = lambda xs: [float(x) for x in xs if x != '']
    pairs = lambda rs: dict(
        pairs=len(rs), identical_core_dna=sum(r['identical_core_dna'] == 'True' for r in rs),
        core_dna_identity_min=min(num(r['core_dna_identity'] for r in rs), default=''),
        protein_identity_min=min(num(r['protein_identity_min'] for r in rs), default=''),
        core_dna_coverage_lt_0_9=sum(r['core_dna_coverage_min'] != '' and float(r['core_dna_coverage_min']) < 0.9 for r in rs),
        no_core_dna_hit=sum(r['core_dna_identity'] == '' for r in rs),
        pair_side_cover=dict(Counter('%s/%s' % (r['up_cover'], r['down_cover']) for r in rs).most_common(6)))
    return dict(
        pairs=pairs(bp), pairs_by_s6_grouping={g: pairs([r for r in bp if r['s6_grouping'] == g])
                                               for g in ('supported_positions', 'grouping_unresolved')},
        occurrences=len(bo), both_sides_resolved=sum(r['boundary_resolved'] == 'True' for r in bo),
        span_class=dict(Counter(r['span_class'] for r in bo)),
        occurrences_with_s6_supported_partner=sum(r['partners_s6_supported'] != '0' for r in bo),
        span_class_with_s6_supported_partner=dict(Counter(r['span_class'] for r in bo if r['partners_s6_supported'] != '0')),
        side_status=dict(Counter(r['up_status'] for r in bo) + Counter(r['down_status'] for r in bo)),
        single_partner=sum(r['partners'] == '1' for r in bo))


def cmd_summarize(a):
    batches = [int(x) for x in a.batches.split(',')]
    tag = 'b' + '-'.join(map(str, batches))
    root = a.work / 'review' / 's6' / tag
    C, P, X, O = [], [], [], []
    missing = []
    for k in range(se.SHARDS):
        d = root / ('shard_%03d' % k)
        if not (d / 'run.json').exists():
            missing.append(k); continue
        C += pc.read_tsv(d / 'cores.tsv.gz')
        P += pc.read_tsv(d / 'same_assembly_pairs.tsv.gz')
        X += [dict(structure=r['structure'], core_key=r['core_key'], s6_grouping=r['s6_grouping'])
              for r in pc.read_tsv(d / 'cross_assembly_single_copy_pairs.tsv.gz')]
        O += [dict(structure=r['structure'], s6_status=r['s6_status'], unresolved_reason=r['unresolved_reason'])
              for r in pc.read_tsv(d / 'occurrences.tsv.gz')]
    out = dict(version=VERSION, batches=batches, evidence_source='assembly_only', production_ready=False,
               site_groups_are_independent_insertion_events=False,
               counts_are_pairs_not_sites_or_events=True,
               s6_grouping_unresolved_not_used_for_site_or_event_counts=True, shards_missing=missing, structures={})
    for st in ('single_gene', 'two_gene'):
        cs = [c for c in C if c['structure'] == st]
        ps = [p for p in P if p['structure'] == st]
        xs = [x for x in X if x['structure'] == st]
        os_ = [o for o in O if o['structure'] == st]
        ucd = [p for p in ps if p['ucd_structural_pair'] == 'True']
        hints = {c['core_key']: c['product_hint'] for c in cs}
        k3 = {c['core_key']: int(c['k3_different_pairs']) for c in cs}
        def cores(rows):
            ks = {r['core_key'] for r in rows}
            sup = {r['core_key'] for r in rows if r['s6_grouping'] == 'supported_positions'}
            return dict(cores=len(ks), cores_with_s6_supported_pair=len(sup), cores_only_grouping_unresolved=len(ks - sup),
                        cores_by_hint=dict(Counter(hints[k] for k in ks)),
                        cores_where_k3_had_no_different_pair=sum(k3[k] == 0 for k in ks))
        out['structures'][st] = dict(
            reviewed_cores=len(cs), reviewed_cores_k3_different=sum(int(c['k3_different_pairs']) > 0 for c in cs),
            reviewed_cores_k3_without_different=sum(int(c['k3_different_pairs']) == 0 for c in cs),
            occurrences=len(os_), occurrence_status=dict(Counter(o['s6_status'] for o in os_)),
            unresolved_reasons=dict(Counter(o['unresolved_reason'] for o in os_ if o['unresolved_reason'])),
            site_groups=sum(int(c['site_groups']) for c in cs),
            single_observation_sites=sum(int(c['single_observation_sites']) for c in cs),
            unit_pairs=dict(different=sum(int(c['unit_pairs_different']) for c in cs),
                            unknown=sum(int(c['unit_pairs_unknown']) for c in cs),
                            conflicting=sum(int(c['unit_pairs_conflicting']) for c in cs),
                            ambiguous=sum(int(c['unit_pairs_ambiguous']) for c in cs)),
            cores_with_2plus_pairwise_context_different_units=sum(int(c['max_pairwise_context_different_units']) >= 2 for c in cs),
            same_assembly_pairs=len(ps),
            same_assembly_context=dict(Counter(p['context_class'] for p in ps)),
            same_assembly_physically_shared_host_genes=dict(Counter(p['context_class'] for p in ps
                                                                    if int(p['host_genes_shared_physically']) > 0)),
            ucd_structural_pairs=len(ucd),
            ucd_structural_pairs_by_s6_grouping=dict(Counter(p['s6_grouping'] for p in ucd)),
            ucd_structural_pairs_s6_status=dict(Counter(s for p in ucd for s in (p['s6_status_a'], p['s6_status_b']))),
            ucd_structural=cores(ucd),
            cross_different_position_pairs=sum(int(c['cross_different_position_pairs']) for c in cs),
            cross_single_copy_pairs=len(xs),
            cross_single_copy_pairs_by_s6_grouping=dict(Counter(x['s6_grouping'] for x in xs)),
            cross_single_copy=cores(xs))
    broot = a.work / 'review' / 'boundary' / tag
    BP, BO, bmissing = [], [], []
    for k in range(se.SHARDS):
        d = broot / ('shard_%03d' % k)
        if not (d / 'run.json').exists():
            bmissing.append(k); continue
        BP += pc.read_tsv(d / 'pairs.tsv.gz'); BO += pc.read_tsv(d / 'occurrences.tsv.gz')
    out['boundary'] = dict(shards_missing=bmissing)
    for st in ('single_gene', 'two_gene'):
        out['boundary'][st] = {sc: boundary_stats([r for r in BP if r['structure'] == st and r['pair_scope'] == sc],
                                                  [r for r in BO if r['structure'] == st and r['pair_scope'] == sc])
                               for sc in ('same_assembly', 'cross_assembly_single_copy')}
    a.repo_out.mkdir(parents=True, exist_ok=True)
    if a.regression:
        cands = pc.read_tsv(a.regression)
        loc = {}
        for k in range(se.SHARDS):
            for o in se.read_shard(a.work, batches, k):
                if o['structure'] == 'single_gene':
                    loc[(o['sample'], o['contig'], str(int(o['start0']) + 1), o['end0'])] = o['occurrence_id']
        sa = {frozenset((r['occurrence_a'], r['occurrence_b'])): r for r in P}
        bpair = {frozenset((r['occurrence_a'], r['occurrence_b'])): r for r in BP if r['pair_scope'] == 'same_assembly'}
        rows = []
        for c in cands:
            ids = []
            for x in ('loc_a', 'loc_b'):
                contig, rest = c[x].rsplit(':', 1)
                st_, en = rest.split('(')[0].split('-')
                ids.append(loc.get((c['sample'], contig, st_, en)))
            k = frozenset(ids)
            p_, b_ = sa.get(k, {}), bpair.get(k, {})
            rows.append(dict(candidate=c['candidate'], batch=c['batch'], occurrence_a=ids[0], occurrence_b=ids[1],
                             same_assembly_pair_in_s6=bool(p_), context_class=p_.get('context_class', ''),
                             host_genes_shared_physically=p_.get('host_genes_shared_physically', ''),
                             ucd_structural_pair=p_.get('ucd_structural_pair', ''),
                             s6_grouping=p_.get('s6_grouping', ''), s6_status_a=p_.get('s6_status_a', ''),
                             s6_status_b=p_.get('s6_status_b', ''),
                             ucd_pair_candidate_strict_v0=p_.get('ucd_pair_candidate_strict_v0', ''),
                             core_dna_identity=b_.get('core_dna_identity', ''), core_dna_coverage_min=b_.get('core_dna_coverage_min', ''),
                             protein_identity_min=b_.get('protein_identity_min', ''),
                             side_cover='%s/%s' % (b_['up_cover'], b_['down_cover']) if b_ else ''))
        pc.write_tsv(a.repo_out / ('s6_regression_%s.tsv' % tag), rows, list(rows[0]))
        state = lambda r: (r['s6_grouping'] if r['ucd_structural_pair'] == 'True' else
                           'not_structural' if r['same_assembly_pair_in_s6'] else 'not_a_same_assembly_pair')
        out['regression'] = dict(Counter('%s:%s' % (r['batch'], state(r)) for r in rows))
    pc.write_tsv(a.repo_out / ('s6_cores_%s.tsv.gz' % tag), C, list(C[0]))
    pc.write_tsv(a.repo_out / ('s6_ucd_structural_pairs_%s.tsv.gz' % tag),
                 [p for p in P if p['ucd_structural_pair'] == 'True'], list(P[0]))
    for sc, name in (('same_assembly', 'same_assembly'), ('cross_assembly_single_copy', 'cross_single_copy')):
        bp = [r for r in BP if r['pair_scope'] == sc]
        bo = [r for r in BO if r['pair_scope'] == sc]
        if bo:
            pc.write_tsv(a.repo_out / ('boundary_occurrences_%s_%s.tsv.gz' % (name, tag)), bo, list(bo[0]))
        if bp and sc == 'same_assembly':
            pc.write_tsv(a.repo_out / ('boundary_pairs_%s_%s.tsv.gz' % (name, tag)), bp, list(bp[0]))
    (a.repo_out / ('s6_summary_%s.json' % tag)).write_text(json.dumps(out, indent=1) + '\n')
    print(json.dumps(out, indent=1))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest='mode', required=True)
    s = sub.add_parser('windows')
    s.add_argument('--work', type=Path, required=True)
    s.add_argument('--batch', type=int, required=True)
    for name in ('s6', 'boundary', 'summarize'):
        s = sub.add_parser(name)
        s.add_argument('--work', type=Path, required=True)
        s.add_argument('--batches', required=True)
        if name in ('s6', 'boundary'):
            s.add_argument('--shard', type=int, required=True)
        if name == 'boundary':
            s.add_argument('--manifest', type=Path, required=True)
        if name == 'summarize':
            s.add_argument('--repo-out', type=Path, required=True)
            s.add_argument('--regression', type=Path)
    a = p.parse_args()
    dict(windows=cmd_windows, s6=cmd_s6, boundary=cmd_boundary, summarize=cmd_summarize)[a.mode](a)


if __name__ == '__main__':
    main()
