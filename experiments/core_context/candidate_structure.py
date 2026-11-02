"""Full-insertion check of the candidate structures (not_in_assembly:core_absent) of a read-site search batch: the
inserted sequence from several molecules, its core coverage, its gene composition and its boundary evidence.
Three levels are reported apart:
  core hit         the read core against every panel sequence of the core: coverage of the panel sequence, identity,
                   and whether both ends of the panel sequence are covered (core_both_ends)
  structure        molecules joining U - insertion - D (strict P/A, candidate_evidence.py) and the boundary counts by
                   molecule (read_site_search insert_fields)
  complete target  the insertion between the two junctions is < MAX_INSERT bp and holds 1-2 complete genes
Inserted sequence: a consensus, never one read. The draft is the P allele window (the best record's insert +/- PAD bp
of assembly) from candidate_evidence; every record of a molecule the strict check called P (P or P_matched) is aligned
to it (minimap2 -a, read preset) and samtools consensus (simple mode, call fraction CALL_FRACTION, min depth
MIN_DEPTH, insertions shown, deletions removed) gives a new draft; ROUNDS rounds. Per-base depth of the last round is
reported (min over the insertion). The junctions are where PAD bp of assembly U ends and PAD bp of D begins in the
consensus (blastn, best HSP reaching the inner end); the insertion is the sequence between them (a target-site
duplication stays in U and D).
Genes: ORFs of the insertion on both strands (translation table 11, starts ATG/GTG/TTG, >= MIN_ORF_AA aa, longest
first, overlaps <= ORF_OVERLAP bp allowed); complete = start and stop inside the insertion. Each ORF is searched
(blastp) against the annotated proteins of the sample and every protein of the core's families (pool proteins);
best hit gives family, homology group and product; a complete gene = complete ORF whose best hit covers >=
GENE_COVERAGE of the subject and of the ORF. Gene count = complete genes; partial ORFs and ORFs without a hit are
listed, not counted. Positions where the consensus is ambiguous (N: depth < MIN_DEPTH or no base reaches
CALL_FRACTION) are counted (insert_ambiguous_bases). genes_status: resolved (none); resolved_ambiguous_bases_in_genes
(every one lies inside a complete gene: an N does not shift the frame, the gene is still full length); otherwise
unresolved_ambiguous_consensus, and then the gene count is not interpreted (target_1_2_genes unresolved): a frameshift
or stop of a low-depth consensus is never read as a broken gene. Assembled copies: the insertion is aligned to the sample's own assembly (minimap2 asm20, all
copies); copies covering >= GENE_COVERAGE of the insertion are listed with the annotated CDS inside them and whether
they are catalog occurrences, so the consensus genes can be compared with an annotated copy.
Not insertion candidates: not_in_assembly:junction positions are listed apart as structural anomalies.
Usage: python candidate_structure.py --out <batch repo dir> --work <batch scratch dir> [--threads N]
"""
import argparse
import gzip
import json
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

from Bio.Seq import Seq

import protein_context as pc
import read_site_search as rs
import s3_readcheck as rc3
import ucd_readcheck as u

PAD = 300  # bp of assembly kept on each side of the insertion in the consensus
ROUNDS = 2
CALL_FRACTION = 0.6
MIN_DEPTH = 3
MIN_ORF_AA = 50
ORF_OVERLAP = 60
GENE_COVERAGE = 0.9
FAMILIES = u.PROJ / 'core_context' / 'scale_comamonadaceae_v0.1' / 'families' / 'catalog.tsv.gz'
POOL_PROTEINS = u.PROJ / 'core_context' / 'scale_comamonadaceae_v0.1' / 'pool_proteins.faa.gz'
GENES = u.PROJ / 'core_context' / 'scale_comamonadaceae_v0.1' / 'genes'


def run(cmd, out=None):
    if out:
        with open(str(out), 'w') as h:
            subprocess.run(cmd, stdout=h, stderr=subprocess.DEVNULL, check=True)
    else:
        return subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, check=True,
                              universal_newlines=True).stdout


def fasta(path):
    return rc3.read_fasta(str(path))


def consensus(draft, reads_fa, preset, wd, threads):
    """ROUNDS of align + samtools consensus; returns (sequence, per-base depth of the last round)."""
    cur = draft
    for i in range(ROUNDS):
        ref = wd / ('draft%d.fa' % i)
        ref.write_text('>draft\n%s\n' % cur)
        bam = wd / ('round%d.bam' % i)
        subprocess.run('minimap2 -ax %s -t %d %s %s 2>/dev/null | samtools sort -o %s - 2>/dev/null && samtools index %s'
                       % (preset, threads, ref, reads_fa, bam, bam), shell=True, check=True)
        out = wd / ('cons%d.fa' % i)
        run(['samtools', 'consensus', '-a', '--mode', 'simple', '-c', str(CALL_FRACTION), '-d', str(MIN_DEPTH),
             '--show-ins', 'yes', '--show-del', 'no', '-f', 'fasta', str(bam)], out)
        cur = ''.join(fasta(out).values()).upper()
    ref = wd / 'final.fa'
    ref.write_text('>cons\n%s\n' % cur)
    bam = wd / 'final.bam'
    subprocess.run('minimap2 -ax %s -t %d %s %s 2>/dev/null | samtools sort -o %s - 2>/dev/null && samtools index %s'
                   % (preset, threads, ref, reads_fa, bam, bam), shell=True, check=True)
    depth = [0] * len(cur)
    for line in run(['samtools', 'depth', '-a', str(bam)]).splitlines():
        f = line.split('\t')
        depth[int(f[1]) - 1] = int(f[2])
    return cur, depth


def blastn(q, s, wd, name):
    qf, sf = wd / (name + '_q.fa'), wd / (name + '_s.fa')
    qf.write_text('>q\n%s\n' % q)
    sf.write_text('>s\n%s\n' % s)
    out = run(['blastn', '-query', str(qf), '-subject', str(sf), '-outfmt',
               '6 pident length qstart qend sstart send qlen slen bitscore'])
    return [dict(zip(('pident', 'length', 'qstart', 'qend', 'sstart', 'send', 'qlen', 'slen', 'bitscore'),
                     map(float, l.split('\t')))) for l in out.splitlines()]


def junctions(cons, useq, dseq, wd):
    """(insert start, insert end) in the consensus: end of U's best HSP reaching U's inner end, start of D's best HSP
    reaching D's inner end (plus strand, the consensus is in assembly orientation); None when not found."""
    ue = [h for h in blastn(useq, cons, wd, 'u') if h['qend'] >= len(useq) - 5 and h['sstart'] < h['send']]
    ds = [h for h in blastn(dseq, cons, wd, 'd') if h['qstart'] <= 6 and h['sstart'] < h['send']]
    if not ue or not ds:
        return None
    a = int(max(ue, key=lambda h: h['bitscore'])['send'])
    b = int(max(ds, key=lambda h: h['bitscore'])['sstart']) - 1
    return a, b


def orfs(seq):
    """ORFs >= MIN_ORF_AA on both strands, complete ones first, then longest first, without overlaps > ORF_OVERLAP: start, end
    (insertion coordinates, 0-based half-open), strand, protein, complete. In each stop-to-stop segment the ORF starts
    at the first start codon (complete when the segment ends with a stop inside the insertion); a segment that runs
    off either end of the insertion is partial (no start codon before it reaches the start, or no stop)."""
    out, n = [], len(seq)
    for strand, s in (('+', seq), ('-', str(Seq(seq).reverse_complement()))):
        for frame in range(3):
            codons = [s[i:i + 3] for i in range(frame, n - 2, 3)]
            i = 0
            while i < len(codons):
                j = i
                while j < len(codons) and codons[j] not in ('TAA', 'TAG', 'TGA'):
                    j += 1
                seg, has_stop = codons[i:j], j < len(codons)
                st = next((k for k, c in enumerate(seg) if c in ('ATG', 'GTG', 'TTG')), None)
                k0 = st if st is not None else (0 if i == 0 else None)  # no start: partial only at the left end
                if k0 is not None and len(seg) - k0 >= MIN_ORF_AA:
                    a, b = frame + 3 * (i + k0), frame + 3 * j + (3 if has_stop else 0)
                    prot = str(Seq(''.join(seg[k0:])).translate(table=11))
                    if st is not None:
                        prot = 'M' + prot[1:]
                    lo, hi = (a, b) if strand == '+' else (n - b, n - a)
                    out.append(dict(start=lo, end=hi, strand=strand, protein=prot, aa=len(prot),
                                    complete=st is not None and has_stop))
                i = j + 1
    keep = []
    for o in sorted(out, key=lambda o: (not o['complete'], -o['aa'], o['start'])):
        if all(min(o['end'], k['end']) - max(o['start'], k['start']) <= ORF_OVERLAP for k in keep):
            keep.append(o)
    return sorted(keep, key=lambda o: o['start'])


def protein_db(samples, core_fams, wd):
    """blastp subject FASTA: annotated proteins of the samples and every pool protein of the core families; with
    sha -> (family, homology group) and sha -> product."""
    want, prod = set(), {}
    for s in samples:
        for g in pc.read_tsv(GENES / ('%s.tsv.gz' % s)):
            if g['type'] == 'CDS' and g['protein_sha256']:
                want.add(g['protein_sha256'])
                prod.setdefault(g['protein_sha256'], g['product'])
    fam = {}
    with gzip.open(str(FAMILIES), 'rt') as h:
        next(h)
        for line in h:
            f = line.split('\t', 4)
            if f[1] in core_fams or f[0] in want:
                fam[f[0]] = (f[1], f[3])
                if f[1] in core_fams:
                    want.add(f[0])
    db = wd / 'proteins.faa'
    with gzip.open(str(POOL_PROTEINS), 'rt') as h, open(str(db), 'w') as o:
        keep = False
        for line in h:
            if line.startswith('>'):
                keep = line[1:].split()[0] in want
            if keep:
                o.write(line)
    return db, fam, prod


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--work', type=Path, default=rs.WORK)
    p.add_argument('--out', type=Path, default=rs.REPO / 'read_site_search_pilot')
    p.add_argument('--threads', type=int, default=16)
    a = p.parse_args()
    od = a.out / 'candidates'
    summ = list(pc.read_tsv(od / 'evidence_summary.tsv'))
    mols = defaultdict(list)
    for r in pc.read_tsv(od / 'evidence_molecules.tsv'):
        mols[r['candidate']].append(r)
    pos = list(pc.read_tsv(a.out / 'positions.tsv.gz'))
    members = {m['member']: m for m in pc.read_tsv(a.work / 'panel' / 'members.tsv')}
    panel = rs.panel_members(a.work / 'panel' / 'panel.fa', set(members))
    catalog = list(pc.read_tsv(a.work / 'panel' / 'catalog.tsv'))
    tested = [c for c in summ if c['status'] == 'tested']
    core_fams = {f for c in tested for f in rs.core_families(c['core_key'])}
    wd = a.work / 'structure'
    wd.mkdir(parents=True, exist_ok=True)
    db, fam, prod = protein_db(sorted({c['sample'] for c in tested}), core_fams, wd)
    rows, orf_rows, copy_rows, fa = [], [], [], []
    for c in tested:
        cid, s = c['candidate'], c['sample']
        cw = wd / cid.replace(':', '_')
        cw.mkdir(exist_ok=True)
        runs = u.run_info()[s]
        preset = rc3.PRESET[runs[0]['read_type']]
        # records of the molecules the strict check called P
        names = {x.split(':')[0] for r in mols[cid] if r['strict_call'] in ('P', 'P_matched')
                 for x in r['record_calls'].split(';') if x.split(':', 1)[1].startswith('P')}
        rfa = cw / 'reads.fa'
        with open(str(rfa), 'w') as o:
            for r in runs:
                with gzip.open(r['local_file'], 'rt') as h:
                    keep = False
                    for i, line in enumerate(h):
                        if i % 4 == 0:
                            keep = line[1:].split(None, 1)[0] in names
                            if keep:
                                o.write('>' + line[1:].split(None, 1)[0] + '\n')
                        elif i % 4 == 1 and keep:
                            o.write(line)
        al = fasta(a.work / 'candidates' / s / 'alleles.fa')
        pal = al['%s|C|P' % cid]
        L = int(c['insert_len'])
        F = rc3.FLANK
        draft = pal[F - PAD:F + L + PAD]
        cons, depth = consensus(draft, rfa, preset, cw, a.threads)
        contig, lp, rp = c['contig'], int(c['left_pos']), int(c['right_pos'])
        seq = fasta(rs.manifest()[s]['fasta'])[contig]
        j = junctions(cons, seq[lp - PAD:lp], seq[rp:rp + PAD], cw)
        row = dict(candidate=cid, sample=s, core_key=c['core_key'], P_molecules=sum(
            r['strict_call'] in ('P', 'P_matched') for r in mols[cid]), records_used=len(names),
            A_matched=c['A_matched'], P_matched=c['P_matched'], boundary=c.get('boundary', ''),
            both_agree_molecules=c.get('both_agree_molecules', ''), boundary_molecules=c.get('boundary_molecules', ''),
            search_insert_len=L)
        if j is None:
            rows.append(dict(row, status='junction_not_found_in_consensus'))
            continue
        ins = cons[j[0]:j[1]]
        fa.append((cid, ins))
        dep = depth[j[0]:j[1]] or [0]
        # core hit level: every panel sequence of this core against the insertion
        best = None
        for name, m in members.items():
            if m['core_key'] != c['core_key'] or not name.startswith('core|'):
                continue
            for h in blastn(panel[name], ins, cw, 'core'):
                k = (h['bitscore'], name)
                if best is None or k > best[0]:
                    best = (k, name, h)
        core = {}
        if best:
            h = best[2]
            core = dict(core_member=best[1], core_member_len=int(h['qlen']),
                        core_coverage=round((h['qend'] - h['qstart'] + 1) / h['qlen'], 3), core_identity=h['pident'],
                        core_both_ends=h['qstart'] <= 10 and h['qend'] >= h['qlen'] - 10,
                        core_in_insert='%d-%d' % (min(h['sstart'], h['send']), max(h['sstart'], h['send'])))
        # genes
        os_ = orfs(ins)
        qf = cw / 'orfs.faa'
        qf.write_text(''.join('>orf%d\n%s\n' % (i, o['protein'].rstrip('*')) for i, o in enumerate(os_)))
        hits = {}
        if os_:
            out = run(['blastp', '-query', str(qf), '-subject', str(db), '-outfmt',
                       '6 qseqid sseqid pident length qstart qend sstart send qlen slen evalue bitscore',
                       '-evalue', '1e-5', '-max_target_seqs', '5'])
            for l in out.splitlines():
                f = l.split('\t')
                if f[0] not in hits or float(f[11]) > float(hits[f[0]][11]):
                    hits[f[0]] = f
        genes, gene_iv = 0, []
        for i, o in enumerate(os_):
            h = hits.get('orf%d' % i)
            fm = fam.get(h[1], ('', '')) if h else ('', '')
            scov = (int(h[7]) - int(h[6]) + 1) / int(h[9]) if h else 0
            qcov = (int(h[5]) - int(h[4]) + 1) / int(h[8]) if h else 0
            full = bool(o['complete'] and h and scov >= GENE_COVERAGE and qcov >= GENE_COVERAGE)
            genes += full
            if full:
                gene_iv.append((o['start'], o['end']))
            orf_rows.append(dict(candidate=cid, orf='orf%d' % i, start=o['start'], end=o['end'], strand=o['strand'],
                                 aa=o['aa'], complete_orf=o['complete'], best_hit=h[1][:16] if h else '',
                                 identity=h[2] if h else '', subject_coverage=round(scov, 3), orf_coverage=round(qcov, 3),
                                 family=fm[0], homology_group=fm[1], core_family=fm[0] in core_fams,
                                 product=prod.get(h[1], '') if h else '', complete_gene=full))
        # assembled copies of the insertion in the sample's own assembly
        qi = cw / 'insert.fa'
        qi.write_text('>ins\n%s\n' % ins)
        out = run(['minimap2', '-cx', 'asm20', '--secondary=yes', '-N', '50', '-p', '0.5', rs.manifest()[s]['fasta'],
                   str(qi)])
        gl = defaultdict(list)
        for g in pc.read_tsv(GENES / ('%s.tsv.gz' % s)):
            if g['type'] == 'CDS':
                gl[g['contig']].append(g)
        copies = []
        for l in out.splitlines():
            f = l.split('\t')
            cov = (int(f[3]) - int(f[2])) / int(f[1])
            if cov < GENE_COVERAGE:
                continue
            t, ts, te = f[5], int(f[7]), int(f[8])
            inside = [g for g in gl[t] if int(g['start0']) >= ts - 50 and int(g['end0']) <= te + 50]
            cat = [x['occurrence_id'] for x in catalog if x['sample'] == s and x['contig'] == t
                   and int(x['start0']) < te and int(x['end0']) > ts]
            copies.append(dict(candidate=cid, contig=t, start=ts, end=te, strand=f[4], insert_coverage=round(cov, 3),
                               identity=round(int(f[9]) / int(f[10]), 3), annotated_cds=len(inside),
                               complete_cds=sum(g['complete_cds'] == 'True' and g['partial'] == 'False'
                                                and g['pseudo'] == 'False' for g in inside),
                               products=';'.join(g['product'] for g in inside), catalog=';'.join(cat)))
        copy_rows += copies
        amb_at = [k for k, ch in enumerate(ins) if ch not in 'ACGT']
        amb = len(amb_at) if amb_at and not all(any(x <= k < y for x, y in gene_iv) for k in amb_at) else 0
        rows.append(dict(row, status='checked', insert_len=len(ins), insert_ambiguous_bases=len(amb_at),
                         ambiguous_positions=','.join(map(str, amb_at)),
                         genes_status='resolved' if not amb_at else 'resolved_ambiguous_bases_in_genes' if not amb
                         else 'unresolved_ambiguous_consensus', insert_size='le_%d' % rs.MAX_INSERT
                         if len(ins) <= rs.MAX_INSERT else 'gt_%d' % rs.MAX_INSERT,
                         insert_min_depth=min(dep), insert_median_depth=sorted(dep)[len(dep) // 2],
                         **core, orfs=len(os_), complete_genes=genes,
                         target_1_2_genes=('unresolved' if amb else 1 <= genes <= 2 and len(ins) <= rs.MAX_INSERT),
                         assembled_copies=len(copies),
                         assembled_copy_best=max(copies, key=lambda x: x['identity'])['contig'] + ':%d-%d' % (
                             max(copies, key=lambda x: x['identity'])['start'],
                             max(copies, key=lambda x: x['identity'])['end']) if copies else ''))
    cols = ['candidate', 'sample', 'core_key', 'status', 'P_molecules', 'records_used', 'P_matched', 'A_matched',
            'boundary', 'boundary_molecules', 'both_agree_molecules', 'search_insert_len', 'insert_len', 'insert_size',
            'insert_min_depth', 'insert_median_depth', 'insert_ambiguous_bases', 'ambiguous_positions', 'genes_status', 'core_member', 'core_member_len', 'core_coverage',
            'core_identity', 'core_both_ends', 'core_in_insert', 'orfs', 'complete_genes', 'target_1_2_genes',
            'assembled_copies', 'assembled_copy_best']
    pc.write_tsv(od / 'structure_summary.tsv', rows, cols)
    if orf_rows:
        pc.write_tsv(od / 'structure_orfs.tsv', orf_rows, list(orf_rows[0]))
    if copy_rows:
        pc.write_tsv(od / 'structure_assembled_copies.tsv', copy_rows, list(copy_rows[0]))
    (od / 'inserts_consensus.fa').write_text(''.join('>%s\n%s\n' % x for x in fa))
    junc = [r for r in pos if int(r['flank']) == rs.MAIN_FLANK and r['cls'] == 'not_in_assembly:junction'
            and int(r['molecules']) >= rs.MAIN_COUNT]
    pc.write_tsv(od / 'structural_anomalies.tsv', junc, list(pos[0]) if pos else [])
    (od / 'structure.json').write_text(json.dumps(dict(
        pad=PAD, rounds=ROUNDS, call_fraction=CALL_FRACTION, min_depth=MIN_DEPTH, min_orf_aa=MIN_ORF_AA,
        gene_coverage=GENE_COVERAGE, candidates=len(tested), checked=sum(r['status'] == 'checked' for r in rows),
        structural_anomalies_not_treated_as_insertions=len(junc)), indent=1) + '\n')
    for r in rows:
        print({k: r.get(k, '') for k in cols})


if __name__ == '__main__':
    main()
