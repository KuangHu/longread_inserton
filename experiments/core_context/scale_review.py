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
          U-core-D pair candidates (the frozen-candidate rule, now over every core): one assembly, both copies
          supported positions, host context different on both sides, no host gene shared physically.
boundary  per shard, on the U-core-D pair candidates only: core protein identity (sha256 or blastp) and core DNA identity /
          coverage, and the coding span +/- W bp of both copies aligned (core_review.blast / chain / extensions,
          unchanged): a side is measured only if the core end is covered in both copies, unknown otherwise (never 0).
          Per occurrence, core_review.side_boundary over all its U-core-D partners; estimated span and the strict
          < 2500 bp re-check only when both sides are resolved. Partners are the U-core-D pair candidates of that
          assembly, not every different-context partner, so "resolved" means resolved among those partners.
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
    sa_rows = []
    hostg = lambda o: {x[2] for s in ('up', 'down') for x in o['host'][s]}
    for k, v in verdicts.items():
        if v['verdict'] != 'different_position_same_assembly':
            continue
        x, y = sorted(k)
        X, Y = occ[x], occ[y]
        gap = (max(int(X['start0']), int(Y['start0'])) - min(int(X['end0']), int(Y['end0']))) if X['contig'] == Y['contig'] else ''
        phys = len(hostg(X) & hostg(Y))
        ctx = rv.context_class(v['up_side'], v['down_side'])
        both = res[x][0] in SUPPORTED and res[y][0] in SUPPORTED
        sa_rows.append(dict(occurrence_a=x, occurrence_b=y, sample=X['sample'], context_class=ctx,
                            up_support=v['up_support'], down_support=v['down_support'], copy_gap_bp=gap,
                            host_genes_shared_physically=phys, both_supported_positions=both,
                            ucd_pair_candidate=both and ctx == 'context_different_both_sides' and phys == 0))
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
                   ucd_pair_candidates=sum(r['ucd_pair_candidate'] for r in sa_rows),
                   ucd_assemblies=len({r['sample'] for r in sa_rows if r['ucd_pair_candidate']}))
    return occ_rows, groups, sa_rows, rel_rows, summary


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
    O, G, P, R, C = [], [], [], [], []
    for key in sorted(groups_of):
        ys = sorted(groups_of[key], key=lambda o: o['occurrence_id'])
        occ_rows, groups, sa_rows, rel_rows, summary = review_core(ys)
        base = dict(structure=key[0], core_key=key[1])
        O += [dict(base, **r) for r in occ_rows]
        for gi, members in enumerate(groups):
            G.append(dict(base, unit='G%d' % gi, members=len(members), occurrences=';'.join(members)))
        P += [dict(base, **r) for r in sa_rows]
        R += [dict(base, **r) for r in rel_rows]
        c = keys[key]
        C.append(dict(base, products=c['products'], product_hint=hint(c['products'].split(' | ')[0]),
                      k3_different_pairs=c['different_context_candidate'],
                      k3_different_pairs_same_assembly=c['different_context_candidate_same_assembly'],
                      k3_insufficient_pairs=c['insufficient_context_labels'], k3_edge_pairs=c['edge_or_missing_context'],
                      max_copies_in_one_assembly=c['max_copies_in_one_assembly'], **summary))
    for name, rows in (('occurrences', O), ('site_groups', G), ('same_assembly_pairs', P), ('group_relations', R),
                       ('cores', C)):
        pc.write_tsv(out / (name + '.tsv.gz'), rows, list(rows[0]) if rows else ['structure'])
    (out / 'run.json').write_text(json.dumps(dict(shard=a.shard, batches=batches, reviewed_cores=len(C),
                                                  occurrences=len(O)), indent=1) + '\n')


def cmd_boundary(a):
    batches = [int(x) for x in a.batches.split(',')]
    tag = 'b' + '-'.join(map(str, batches))
    d = a.work / 'review' / 's6' / tag / ('shard_%03d' % a.shard)
    pairs = [r for r in pc.read_tsv(d / 'same_assembly_pairs.tsv.gz') if r['ucd_pair_candidate'] == 'True']
    need = {x for r in pairs for x in (r['occurrence_a'], r['occurrence_b'])}
    occ = {o['occurrence_id']: o for o in se.read_shard(a.work, batches, a.shard) if o['occurrence_id'] in need}
    fasta = {m['sample']: m['fasta'] for m in pc.read_tsv(a.manifest)}
    groups = defaultdict(list)
    for r in pairs:
        groups[(r['structure'], r['core_key'], r['sample'])].append(r)
    P, partners = [], defaultdict(list)
    for smp in sorted({k[2] for k in groups}):
        seqs = pc.read_fna(fasta[smp])
        genes = {g['gene_uid']: g for g in pc.read_tsv(a.work / 'genes' / (smp + '.tsv.gz'))}
        prot = se.read_faa_gz(a.work / 'genes' / (smp + '.faa.gz'))
        for key in sorted(k for k in groups if k[2] == smp):
            ids = sorted({x for r in groups[key] for x in (r['occurrence_a'], r['occurrence_b'])})
            num = {x: str(i) for i, x in enumerate(ids)}
            for x in ids:
                o = occ[x]
                g = genes[o['genes'].split(';')[0]]
                o['ext'], o['c0'], o['c1'] = cr.window(seqs[o['contig']], int(o['start0']), int(o['end0']), o['strand'] == '+',
                                                       cr.W, g['circular_contig'] == 'True')
                o['core_dna'] = o['ext'][o['c0']:o['c1']]
                o['prots'] = [prot[genes[u]['protein_sha256']] for u in cr.core_genes(o)]
            npos = 1 if key[0] == 'single_gene' else 2
            php = [cr.blast('blastp', {num[x] + '_%d' % i: occ[x]['prots'][i] for x in ids}, None) for i in range(npos)]
            dcore = cr.blast('blastn', {num[x]: occ[x]['core_dna'] for x in ids}, 'blastn')
            dext = cr.blast('blastn', {num[x]: occ[x]['ext'] for x in ids}, 'blastn')
            for r in groups[key]:
                A, B = occ[r['occurrence_a']], occ[r['occurrence_b']]
                ia, ib = num[A['occurrence_id']], num[B['occurrence_id']]
                pb = [cr.best(php[i][(ia + '_%d' % i, ib + '_%d' % i)]) for i in range(npos)]
                db = cr.best(dcore[(ia, ib)])
                blk = cr.chain([(h['qs'], h['qe'], h['ss'], h['se']) for h in dext[(ia, ib)]], A['c0'], A['c1'], B['c0'], B['c1'])
                ext = cr.extensions(blk, A['c0'], A['c1'], len(A['ext']), B['c0'], B['c1'], len(B['ext'])) if blk else \
                    dict(up_cover='no_core_alignment', down_cover='no_core_alignment', a_up=None, b_up=None, a_down=None,
                         b_down=None, up_reaches_end=False, down_reaches_end=False)
                P.append(dict(structure=key[0], core_key=key[1], sample=smp, occurrence_a=A['occurrence_id'],
                    occurrence_b=B['occurrence_id'], identical_protein=A['prots'] == B['prots'],
                    protein_identity_min=min(x[0] for x in pb) if all(pb) else '',
                    protein_coverage_min=round(min(min(x[1], x[2]) for x in pb), 4) if all(pb) else '',
                    identical_core_dna=A['core_dna'].upper() == B['core_dna'].upper(),
                    core_dna_identity=db[0] if db else '', core_dna_coverage_min=round(min(db[1], db[2]), 4) if db else '',
                    up_cover=ext['up_cover'], a_up=ext['a_up'], b_up=ext['b_up'],
                    down_cover=ext['down_cover'], a_down=ext['a_down'], b_down=ext['b_down']))
                for o, pre in ((A, 'a'), (B, 'b')):
                    partners[o['occurrence_id']].append((ext[pre + '_up'], ext['up_reaches_end'], ext[pre + '_down'],
                                                         ext['down_reaches_end']))
    O = []
    for x in sorted(partners):
        o, ps = occ[x], partners[x]
        up = cr.side_boundary([p[0] for p in ps], any(p[1] for p in ps))
        down = cr.side_boundary([p[2] for p in ps], any(p[3] for p in ps))
        span = int(o['span'])
        res = up[1] == down[1] == 'resolved'
        est = span + up[0] + down[0] if res else None
        O.append(dict(structure=o['structure'], core_key=o['core_key'], occurrence_id=x, sample=o['sample'],
                      coding_span=span, ucd_partners=len(ps), up_status=up[1], ext_up='' if up[0] is None else up[0],
                      up_partners_end_not_covered=up[2], down_status=down[1], ext_down='' if down[0] is None else down[0],
                      down_partners_end_not_covered=down[2], boundary_resolved=res,
                      estimated_span='' if est is None else est, strictly_lt2500='' if est is None else est < pc.TARGET_MAX))
    out = a.work / 'review' / 'boundary' / tag / ('shard_%03d' % a.shard)
    out.mkdir(parents=True, exist_ok=True)
    pc.write_tsv(out / 'pairs.tsv.gz', P, list(P[0]) if P else ['structure'])
    pc.write_tsv(out / 'occurrences.tsv.gz', O, list(O[0]) if O else ['structure'])
    (out / 'run.json').write_text(json.dumps(dict(shard=a.shard, pairs=len(P), occurrences=len(O)), indent=1) + '\n')


def cmd_summarize(a):
    batches = [int(x) for x in a.batches.split(',')]
    tag = 'b' + '-'.join(map(str, batches))
    root = a.work / 'review' / 's6' / tag
    C, P, O = [], [], []
    missing = []
    for k in range(se.SHARDS):
        d = root / ('shard_%03d' % k)
        if not (d / 'run.json').exists():
            missing.append(k); continue
        C += pc.read_tsv(d / 'cores.tsv.gz')
        P += [r for r in pc.read_tsv(d / 'same_assembly_pairs.tsv.gz')]
        O += [dict(structure=r['structure'], s6_status=r['s6_status'], unresolved_reason=r['unresolved_reason'])
              for r in pc.read_tsv(d / 'occurrences.tsv.gz')]
    out = dict(version=VERSION, batches=batches, evidence_source='assembly_only', production_ready=False,
               site_groups_are_independent_insertion_events=False, shards_missing=missing, structures={})
    for st in ('single_gene', 'two_gene'):
        cs = [c for c in C if c['structure'] == st]
        ps = [p for p in P if p['structure'] == st]
        os_ = [o for o in O if o['structure'] == st]
        ucd = [p for p in ps if p['ucd_pair_candidate'] == 'True']
        ucd_cores = {p['core_key'] for p in ucd}
        hints = {c['core_key']: c['product_hint'] for c in cs}
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
            ucd_pair_candidates=len(ucd), ucd_cores=len(ucd_cores),
            ucd_cores_by_hint=dict(Counter(hints[k] for k in ucd_cores)),
            ucd_cores_where_k3_had_no_different_pair=sum(1 for c in cs if c['core_key'] in ucd_cores
                                                          and int(c['k3_different_pairs']) == 0))
    broot = a.work / 'review' / 'boundary' / tag
    BP, BO, bmissing = [], [], []
    for k in range(se.SHARDS):
        d = broot / ('shard_%03d' % k)
        if not (d / 'run.json').exists():
            bmissing.append(k); continue
        BP += pc.read_tsv(d / 'pairs.tsv.gz'); BO += pc.read_tsv(d / 'occurrences.tsv.gz')
    out['boundary'] = dict(shards_missing=bmissing)
    for st in ('single_gene', 'two_gene'):
        bp = [r for r in BP if r['structure'] == st]
        bo = [r for r in BO if r['structure'] == st]
        num = lambda xs: [float(x) for x in xs if x != '']
        out['boundary'][st] = dict(
            pairs=len(bp), identical_core_dna=sum(r['identical_core_dna'] == 'True' for r in bp),
            core_dna_identity_min=min(num(r['core_dna_identity'] for r in bp), default=''),
            core_dna_coverage_lt_0_9=sum(r['core_dna_coverage_min'] != '' and float(r['core_dna_coverage_min']) < 0.9 for r in bp),
            no_core_dna_hit=sum(r['core_dna_identity'] == '' for r in bp),
            pair_side_cover=dict(Counter('%s/%s' % (r['up_cover'], r['down_cover']) for r in bp).most_common(6)),
            occurrences=len(bo), both_sides_resolved=sum(r['boundary_resolved'] == 'True' for r in bo),
            resolved_strictly_lt2500=sum(r['strictly_lt2500'] == 'True' for r in bo),
            side_status=dict(Counter(r['up_status'] for r in bo) + Counter(r['down_status'] for r in bo)),
            single_partner=sum(r['ucd_partners'] == '1' for r in bo))
    if a.regression:
        cands = pc.read_tsv(a.regression)
        loc = {}
        for k in range(se.SHARDS):
            for o in se.read_shard(a.work, batches, k):
                if o['structure'] == 'single_gene':
                    loc[(o['sample'], o['contig'], str(int(o['start0']) + 1), o['end0'])] = o['occurrence_id']
        ucd = {frozenset((r['occurrence_a'], r['occurrence_b'])): r for r in P}
        bpair = {frozenset((r['occurrence_a'], r['occurrence_b'])): r for r in BP}
        rows = []
        for c in cands:
            ids = []
            for x in ('loc_a', 'loc_b'):
                contig, rest = c[x].rsplit(':', 1)
                st_, en = rest.split('(')[0].split('-')
                ids.append(loc.get((c['sample'], contig, st_, en)))
            k = frozenset(ids)
            p_, b_ = ucd.get(k, {}), bpair.get(k, {})
            rows.append(dict(candidate=c['candidate'], batch=c['batch'], occurrence_a=ids[0], occurrence_b=ids[1],
                             same_assembly_pair_in_s6=bool(p_), context_class=p_.get('context_class', ''),
                             host_genes_shared_physically=p_.get('host_genes_shared_physically', ''),
                             ucd_pair_candidate=p_.get('ucd_pair_candidate', ''),
                             core_dna_identity=b_.get('core_dna_identity', ''), core_dna_coverage_min=b_.get('core_dna_coverage_min', ''),
                             protein_identity_min=b_.get('protein_identity_min', ''),
                             side_cover='%s/%s' % (b_['up_cover'], b_['down_cover']) if b_ else ''))
        pc.write_tsv(a.repo_out / ('s6_regression_%s.tsv' % tag), rows, list(rows[0]))
        out['regression'] = dict(Counter('%s:%s' % (r['batch'], r['ucd_pair_candidate'] or 'not_a_same_assembly_pair')
                                         for r in rows))
    a.repo_out.mkdir(parents=True, exist_ok=True)
    pc.write_tsv(a.repo_out / ('s6_cores_%s.tsv.gz' % tag), C, list(C[0]))
    if BO:
        pc.write_tsv(a.repo_out / ('boundary_occurrences_%s.tsv.gz' % tag), BO, list(BO[0]))
        pc.write_tsv(a.repo_out / ('boundary_pairs_%s.tsv.gz' % tag), BP, list(BP[0]))
    pc.write_tsv(a.repo_out / ('s6_ucd_pair_candidates_%s.tsv.gz' % tag),
                 [p for p in P if p['ucd_pair_candidate'] == 'True'], list(P[0]))
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
        else:
            s.add_argument('--repo-out', type=Path, required=True)
            s.add_argument('--regression', type=Path)
    a = p.parse_args()
    dict(windows=cmd_windows, s6=cmd_s6, boundary=cmd_boundary, summarize=cmd_summarize)[a.mode](a)


if __name__ == '__main__':
    main()
