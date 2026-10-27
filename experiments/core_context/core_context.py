"""Gap-independent repeated-block discovery from a verified PanGraph 1.4 JSON.

This is an assembly discovery prototype, not a mobile-element or orthology classifier.
One graph node is one occurrence, including secondary copies within a sample.
Graph block boundaries are provisional; no GFF/CDS or read support is inferred.
Requires Python >=3.10 and Biopython. No FASTQ, gap calls, or empty site is an input.
"""
import argparse
import csv
import hashlib
import itertools
import json
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

from Bio import Align, SeqIO
from Bio.Seq import Seq

VERSION = 'core_context_graph_v0.1'


def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1048576), b''):
            h.update(b)
    return h.hexdigest()


def sha(sequence):
    return hashlib.sha256(sequence.encode()).hexdigest()


def write_json(path, data):
    Path(path).write_text(json.dumps(data, indent=2) + '\n', encoding='utf-8')


def write_tsv(path, rows, columns):
    with open(path, 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, columns, delimiter='\t', lineterminator='\n', extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)


def read_fasta(path):
    out = {}
    for r in SeqIO.parse(path, 'fasta'):
        if r.id in out:
            raise ValueError('Duplicate FASTA id: ' + r.id)
        out[r.id] = str(r.seq).upper()
    return out


def identity(a, b):
    """Global affine identity, denominator includes gaps. Not an orthology test."""
    if not a or not b:
        return None
    aligner = Align.PairwiseAligner(mode='global', match_score=2, mismatch_score=-3,
                                   open_gap_score=-5, extend_gap_score=-2)
    c = aligner.align(a, b)[0].coordinates
    matches = columns = 0
    for j in range(c.shape[1] - 1):
        s, t = map(int, c[:, j]); e, u = map(int, c[:, j + 1])
        columns += max(e - s, u - t)
        if e > s and u > t:
            matches += sum(x == y and x in 'ACGT' for x, y in zip(a[s:e], b[t:u]))
    return matches / columns


def oriented_segments(seq, start, end, strand, flank):
    if not 0 <= start < end <= len(seq):
        raise ValueError('Invalid linear interval')
    up, core, down = seq[max(0, start-flank):start], seq[start:end], seq[end:end+flank]
    if strand == '-':
        up, core, down = (str(Seq(x).reverse_complement()) for x in (down, core, up))
    return up, core, down


def relation(core_id, up_id, down_id, anchors_a, anchors_b, complete,
             min_core=0.9, similar_flank=0.9, low_flank=0.7):
    # Low flank identity alone is deliberately insufficient. Distinct adjacent graph
    # anchors are a second diagnostic; even then orthology remains untested.
    if not complete:
        return 'edge_or_missing_context'
    if core_id is None or core_id < min_core:
        return 'core_similarity_below_threshold'
    if up_id >= similar_flank and down_id >= similar_flank:
        return 'similar_flanks_candidate_same_context'
    known = all(x is not None for x in anchors_a + anchors_b)
    if known and anchors_a[0] != anchors_b[0] and anchors_a[1] != anchors_b[1] and max(up_id, down_id) < low_flank:
        return 'different_adjacency_low_flank_similarity_candidate'
    return 'context_unresolved'


def build(manifest, out, pangraph, threads):
    """Build a fresh graph; sample ids and each contig retain an explicit mapping."""
    out.mkdir(parents=True, exist_ok=False)
    with open(manifest, encoding='utf-8-sig') as f:
        samples = list(csv.DictReader(f, delimiter='\t'))
    if not samples or len({r['sample'] for r in samples}) != len(samples):
        raise ValueError('Manifest needs unique sample ids and fasta paths')
    mapping, inputs = [], []
    with open(out / 'assemblies.fna', 'w') as h:
        for row in samples:
            path = Path(row['fasta'])
            if not path.is_absolute():
                path = Path(manifest).resolve().parent / path
            inputs.append(dict(sample=row['sample'], fasta=str(path.resolve()), sha256=digest(path)))
            for contig, seq in read_fasta(path).items():
                name = 'ctg_' + sha(row['sample'] + '\0' + contig)[:24]
                mapping.append(dict(graph_name=name, sample=row['sample'], contig=contig,
                                    length=len(seq), sha256=sha(seq)))
                h.write('>' + name + '\n' + seq + '\n')
    if len({r['graph_name'] for r in mapping}) != len(mapping):
        raise ValueError('Graph record id collision')
    write_tsv(out / 'contigs.tsv', mapping, ['graph_name', 'sample', 'contig', 'length', 'sha256'])
    version = subprocess.run([pangraph, '--version'], capture_output=True, text=True, check=True).stdout.strip()
    command = [pangraph, 'build', '-j', str(threads), '--verify', '--no-progress-bar',
               str(out / 'assemblies.fna'), '-o', str(out / 'graph.json')]
    with open(out / 'build.stderr.txt', 'w') as err:
        subprocess.run(command, check=True, stderr=err)
    write_json(out / 'graph_provenance.json', dict(tool_version=version, command=command,
               inputs=inputs, graph_sha256=digest(out / 'graph.json'),
               joint_fasta_sha256=digest(out / 'assemblies.fna'),
               contigs_sha256=digest(out / 'contigs.tsv'), graph_verify_passed=True))


def extract(graph, genomes, mapping, flank, min_length, max_length):
    paths, nodes, blocks = graph['paths'], graph['nodes'], graph['blocks']
    grouped, excluded = defaultdict(list), []
    seen = set()
    for path in paths.values():
        name = path['name']
        if name not in mapping or name not in genomes:
            raise ValueError('Unmapped graph path: ' + name)
        seq, meta = genomes[name], mapping[name]
        if len(seq) != path['tot_len'] or len(seq) != int(meta['length']) or sha(seq) != meta['sha256']:
            raise ValueError('FASTA / graph mapping mismatch: ' + name)
        ids = path['nodes']
        for idx, nid in enumerate(ids):
            n = nodes[str(nid)]
            if str(n['path_id']) != str(path['id']) or nid in seen:
                raise ValueError('Node path mismatch or duplicate node')
            seen.add(nid)
            start, end = n['position']
            # Circular and wrapped contexts need a separate implementation; do not
            # silently slice across origins or label the missing side as different.
            if path['circular'] or not 0 <= start < end <= len(seq):
                excluded.append(dict(node_id=nid, reason='circular_or_wrapped_context_pending'))
                continue
            if not min_length <= end-start < max_length:
                continue
            up, core, down = oriented_segments(seq, start, end, n['strand'], flank)
            neighbors = []
            for j in (idx - 1, idx + 1):
                other = nodes[str(ids[j])] if 0 <= j < len(ids) else None
                neighbors.append(None if other is None else (str(other['block_id']), other['strand']))
            if n['strand'] == '-':
                neighbors = [None if x is None else (x[0], '+' if x[1] == '-' else '-') for x in reversed(neighbors)]
            fid = 'PG_' + str(n['block_id'])
            if str(n['block_id']) not in blocks:
                raise ValueError('Missing graph block')
            grouped[fid].append(dict(family_id=fid, occurrence_id='PGN_' + str(nid), sample=meta['sample'],
                contig=meta['contig'], start0=start, end0=end, strand=n['strand'], core_length=end-start,
                core_sha256=sha(core), upstream=up, core=core, downstream=down, anchors=neighbors,
                flank_complete=len(up) == flank and len(down) == flank,
                boundary_status='graph_block_boundary_pending', annotation_status='not_reviewed',
                evidence_source='assembly_only', discovery_source='pangraph_repeated_block', gap_required=False))
    if len(seen) != len(nodes):
        raise ValueError('Orphan graph nodes')
    return {f: sorted(xs, key=lambda x: x['occurrence_id']) for f, xs in grouped.items() if len(xs) >= 2}, excluded


def call(bundle, out, flank=500, min_length=50, max_length=2500, max_pairs=20000):
    out.mkdir(parents=True, exist_ok=False)
    provenance = json.loads((bundle / 'graph_provenance.json').read_text())
    for file, key in [('graph.json', 'graph_sha256'), ('assemblies.fna', 'joint_fasta_sha256'), ('contigs.tsv', 'contigs_sha256')]:
        if digest(bundle / file) != provenance[key]:
            raise ValueError('Bundle hash mismatch: ' + file)
    if provenance.get('graph_verify_passed') is not True:
        raise ValueError('Graph reconstruction verification required')
    graph = json.loads((bundle / 'graph.json').read_text())
    with open(bundle / 'contigs.tsv') as h:
        rows = list(csv.DictReader(h, delimiter='\t'))
    mapping = {r['graph_name']: r for r in rows}
    if len(mapping) != len(rows):
        raise ValueError('Duplicate mapping')
    groups, excluded = extract(graph, read_fasta(bundle / 'assemblies.fna'), mapping, flank, min_length, max_length)
    pairs, family_rows, occurrences = [], [], []
    for family, xs in sorted(groups.items()):
        total = len(xs)*(len(xs)-1)//2
        statuses = Counter()
        for a, b in itertools.islice(itertools.combinations(xs, 2), max_pairs):
            c, u, d = (identity(a[k], b[k]) for k in ('core', 'upstream', 'downstream'))
            status = relation(c, u, d, a['anchors'], b['anchors'], a['flank_complete'] and b['flank_complete'])
            statuses[status] += 1
            pairs.append(dict(family_id=family, occurrence_a=a['occurrence_id'], occurrence_b=b['occurrence_id'],
                core_identity=c, upstream_identity=u, downstream_identity=d, relation=status,
                orthology_status='not_tested', empty_site_required=False))
        family_rows.append(dict(family_id=family, occurrences=len(xs), samples=len({x['sample'] for x in xs}),
            pairs_total=total, pairs_checked=sum(statuses.values()), pair_review_complete=total <= max_pairs,
            different_context_candidate_pairs=statuses['different_adjacency_low_flank_similarity_candidate'],
            similar_context_pairs=statuses['similar_flanks_candidate_same_context']))
        occurrences.extend(xs)
    write_tsv(out / 'occurrences.tsv', occurrences, ['family_id','occurrence_id','sample','contig','start0','end0','strand',
        'core_length','core_sha256','flank_complete','boundary_status','annotation_status','evidence_source','discovery_source','gap_required'])
    write_tsv(out / 'context_pairs.tsv', pairs, ['family_id','occurrence_a','occurrence_b','core_identity','upstream_identity',
        'downstream_identity','relation','orthology_status','empty_site_required'])
    write_tsv(out / 'families.tsv', family_rows, ['family_id','occurrences','samples','pairs_total','pairs_checked',
        'pair_review_complete','different_context_candidate_pairs','similar_context_pairs'])
    write_json(out / 'occurrences.json', occurrences)
    write_json(out / 'excluded_nodes.json', excluded)
    with open(out / 'contexts.fna', 'w') as f:
        for x in occurrences:
            for key in ('upstream','core','downstream'):
                if x[key]:
                    f.write('>' + x['occurrence_id'] + '|' + key + '\n' + x[key] + '\n')
    summary = dict(version=VERSION, discovery_entry='all_repeated_pangraph_blocks', gap_inputs=[], reads_inputs=[],
        bounds=dict(min_length=min_length, max_length_exclusive=max_length, flank=flank),
        diagnostic_thresholds=dict(core_identity=0.9, similar_flank=0.9, low_flank=0.7),
        graph_sha256=provenance['graph_sha256'], bundle_provenance_sha256=digest(bundle/'graph_provenance.json'),
        families=len(groups), occurrences=len(occurrences), context_pairs=len(pairs),
        relations=dict(Counter(p['relation'] for p in pairs)), excluded_nodes=len(excluded),
        truncated_pair_reviews=sum(not r['pair_review_complete'] for r in family_rows),
        annotation_complete=False, site_groups_computed=False, production_ready=False,
        limitations=['Graph blocks can split or combine biological cores; 1-2 CDS not yet checked.',
          'Distinct graph adjacency plus low global flank identity is a diagnostic candidate, not proof of different orthologous sites.',
          'Circular graph paths are deferred. Linear graph paths do not prove biologically linear contigs.',
          'Divergent cores missed by graph construction cannot be recovered here; protein-based discovery remains required.',
          'No mobile activity, novelty, independent event, or long-read validation is inferred.'])
    write_json(out / 'summary.json', summary)
    return summary


def main():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest='mode', required=True)
    b = sub.add_parser('build'); b.add_argument('--manifest', type=Path, required=True)
    b.add_argument('--out', type=Path, required=True); b.add_argument('--pangraph', default='pangraph')
    b.add_argument('--threads', type=int, default=4)
    c = sub.add_parser('call'); c.add_argument('--bundle', type=Path, required=True); c.add_argument('--out', type=Path, required=True)
    c.add_argument('--flank', type=int, default=500); c.add_argument('--min-length', type=int, default=50)
    c.add_argument('--max-length', type=int, default=2500); c.add_argument('--max-pairs-per-family', type=int, default=20000)
    a = p.parse_args()
    if a.mode == 'build':
        if a.threads < 1: p.error('threads must be positive')
        build(a.manifest, a.out, a.pangraph, a.threads)
    else:
        if not 0 < a.min_length < a.max_length or a.flank < 1 or a.max_pairs_per_family < 1:
            p.error('invalid bounds')
        print(json.dumps(call(a.bundle, a.out, a.flank, a.min_length, a.max_length, a.max_pairs_per_family), indent=2))


if __name__ == '__main__':
    main()
