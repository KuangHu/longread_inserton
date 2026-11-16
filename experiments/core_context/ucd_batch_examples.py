"""Write a few real U-X-D examples of a merged ucd_batch run, with their sequence and alignment evidence (Markdown).

Selection rule (fixed, no hand-picking): for each pair scope x structure, among final_status=pass rows of the
candidate_insertion kind, take the pair with the longest insertion (min of both X lengths), then smallest pair_id.
Only reads the merged outputs (all rows, evidence_sequences, evidence_hsps); ucd_batch.py itself is left unchanged.

  python ucd_batch_examples.py --out scale_comamonadaceae_v0.1/ucd_batch_v0.1 \
      --work /global/scratch/users/kh36969/longread_callout/core_context/ucd_batch_v0.1
"""
import argparse
import gzip
from collections import defaultdict
from pathlib import Path

import ucd_compare as uc

J = 40  # bp shown either side of each junction


def _pick(best, r):
    k = (r['pair_scope'], r['structure'])
    key = (-min(int(r['x_len_a']), int(r['x_len_b'])), r['pair_id'])
    if k not in best or key < best[k][0]:
        best[k] = (key, r)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--out', required=True)
    p.add_argument('--work', required=True)
    p.add_argument('--kind', default='candidate_insertion')
    p.add_argument('--rows', default='', help='rows file (default <work>/ucd_all_rows_<kind>.tsv.gz)')
    p.add_argument('--evidence', default='', help='comma-separated dirs with evidence files (default <work>)')
    p.add_argument('--name', default='examples.md')
    a = p.parse_args()
    out, work = Path(a.out), Path(a.work)
    rows_p = Path(a.rows) if a.rows else work / ('ucd_all_rows_%s.tsv.gz' % a.kind)
    ev = [Path(x) for x in (a.evidence.split(',') if a.evidence else [str(work)])]
    best = {}
    with gzip.open(rows_p, 'rt') as f:
        hdr = f.readline().rstrip('\n').split('\t')
        rows = (dict(zip(hdr, line.rstrip('\n').split('\t'))) for line in f)
        rows = (r for r in rows if r['final_status'] == 'pass')
        for r in rows:
            _pick(best, r)
    pick = {r['pair_key']: r for _, r in best.values()}
    sids = {}
    for r in pick.values():
        for s in 'ab':
            lo, hi = r[s + '_x_interval'].split('-')
            sids['%s|%s|%s-%s' % (r[s + '_occurrence'], a.kind, lo, hi)] = None
    seqs = {}
    hsps = defaultdict(list)
    for d in ev:
        with gzip.open(d / 'evidence_sequences.fa.gz', 'rt') as f:
            for h, s in zip(f, f):
                sid, part = h[1:].split()[0].rsplit('|', 1)
                if sid in sids:
                    seqs[(sid, part)] = (h[1:].strip(), s.strip())
        with gzip.open(d / 'evidence_hsps.tsv.gz', 'rt') as f:
            cols = f.readline().rstrip('\n').split('\t')
            for line in f:
                if line.split('\t', 1)[0] in pick:
                    hsps[line.split('\t', 1)[0]].append(dict(zip(cols, line.rstrip('\n').split('\t'))))
    L = ['# Real U–X–D examples (`%s`, final_status = pass)\n' % a.kind,
         'One per pair scope × structure, chosen by a fixed rule: the longest insertion, then the smallest pair_id '
         '(`ucd_batch_examples.py`).',
         '- Sequences come from the Check-1 evidence: re-extracted from the verified source FASTA and identical to the '
         'screen\'s, base by base.',
         '- HSPs are those of Check 2. `chain` marks the collinear HSPs that make up the insertion match.',
         '- Coordinates are 0-based half-open on the contig forward strand. Sequences are shown in the insertion '
         'orientation. B is shown as read; when `b_flipped` is true, B\'s insertion is the reverse complement of A\'s.',
         '- This is a structural pattern only: it says nothing about mobility or an empty site.\n']
    for (_, r) in sorted(best.values(), key=lambda x: (x[1]['pair_scope'], x[1]['structure'])):
        k = r['pair_key']
        L.append('## %s / %s: `%s`\n' % (r['pair_scope'], r['structure'], r['pair_id']))
        L.append('core `%s`, core match `%s` (core DNA identity %s), s6_grouping label `%s` (not used)\n' % (
            r['core_key'], r['core_match_class'], r['core_dna_identity'], r['s6_grouping']))
        L.append('| side | occurrence | assembly | contig (length) | strand | U | X | D | FASTA sha256 |')
        L.append('|---|---|---|---|---|---|---|---|---|')
        for s in 'ab':
            L.append('| %s | `%s` | %s | %s (%s) | %s | %s | %s | %s | `%s…` |' % (
                s.upper(), r[s + '_occurrence'], r[s + '_assembly'], r[s + '_contig'], r[s + '_contig_len'],
                r[s + '_strand'], r[s + '_u_interval'], r[s + '_x_interval'], r[s + '_d_interval'],
                r[s + '_fasta_sha256'][:16]))
        L.append('')
        L.append('- **Insertion:** X/X′ identity %s%%, coverage %s / %s (collinear chain of %s HSP, %s strand); '
                 'lengths %s / %s bp.' % (
                     r['chk_insertion_identity'], r['chk_insertion_cov_a'], r['chk_insertion_cov_b'],
                     r['chk_insertion_chain_hsps'], r['chk_insertion_orientation'], r['x_len_a'], r['x_len_b']))
        L.append('- **Contexts** (coverage A / B by HSPs ≥ 80%): ' + '; '.join(
            '%s %s (%s / %s)' % (c, r['chk_%s_class' % c], r['chk_%s_cov_a' % c], r['chk_%s_cov_b' % c])
            for c in ('UU', 'DD', 'UU_inv', 'DD_inv', 'UD', 'UD_inv', 'DU', 'DU_inv')))
        L.append('- **Verdicts:** screen `%s`; source_check `%s`; ucd_check `%s`.\n' % (
            r['screen_verdict'], r['source_check'], r['ucd_check']))
        L.append('Junctions (%d bp of U | %d bp of X … %d bp of X | %d bp of D):\n' % (J, J, J, J))
        L.append('```')
        for s in 'ab':
            lo, hi = r[s + '_x_interval'].split('-')
            sid = '%s|%s|%s-%s' % (r[s + '_occurrence'], a.kind, lo, hi)
            u, x, d = (seqs[(sid, p)][1] for p in 'UXD')
            L.append('%s  %s | %s ... %s | %s' % (s.upper(), u[-J:].lower(), x[:J], x[-J:], d[:J].lower()))
        L.append('```\n')
        L.append('Check-2 HSPs (q = A, s = B as read):\n')
        L.append('| q | s | strand | identity | length | q range | s range | chain |')
        L.append('|---|---|---|---|---|---|---|---|')
        for h in sorted(hsps[k], key=lambda h: (h['q'] != 'X', h['q'], h['s'], -float(h['bitscore'])))[:12]:
            L.append('| %s | %s | %s | %s | %s | %s–%s | %s–%s | %s |' % (
                h['q'], h['s'], h['strand'], h['pident'], h['length'], h['qs'], h['qe'], h['ss'], h['se'],
                'yes' if h['in_insertion_chain'] == 'True' else ''))
        if not hsps[k]:
            L.append('| – | – | – | – | – | – | – | – |')
        L.append('\nNo HSP other than X/X′ means no reliable similarity between any pair of flanks.\n'
                 if all(h['q'] == 'X' and h['s'] == 'X' for h in hsps[k]) else '')
    (out / a.name).write_text('\n'.join(L) + '\n')
    print('\n'.join(L))


if __name__ == '__main__':
    main()
