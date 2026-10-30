#!/bin/bash
#SBATCH --job-name=rc01_inputs
#SBATCH --account=pc_rubinlab
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --array=1-4
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=8:00:00
#SBATCH --output=/global/scratch/users/kh36969/longread_callout/logs/%x_%A_%a.out
# Input check for one sample: reads and GCA/GCF packages were downloaded and verified on lrc-xfer by 00_prepare.sh.
# Nothing is downloaded here; a missing or unverified input stops the task with input_not_ready (exit 3).
source /global/home/users/kh36969/middle_insertion_callout/longread_inserton/experiments/read_callout/common.sh
FASTQ=$(python $CODE/scripts/reads_ready.py check --runs $RUN | cut -f3 -d' ') || exit 3
for acc in $GCA $GCF; do
    [ -d $ASM_DIR/$acc ] || { echo "input_not_ready assembly $acc missing in $ASM_DIR: run 00_prepare.sh on lrc-xfer" >&2; exit 3; }
done
mkdir -p $OUT
python $HERE/molecule_reads.py diagnose $FASTQ --limit 20000 > $OUT/read_unit.json

# Assemblies: GCA (submitted) and GCF (RefSeq, PGAP annotation used for the callout).
FNA=$(ls $ASM_DIR/$GCF/ncbi_dataset/data/$GCF/*_genomic.fna)
GBFF=$ASM_DIR/$GCF/ncbi_dataset/data/$GCF/genomic.gbff

python - "$ASM_DIR" "$GCA" "$GCF" "$SAMPLE" "$FNA" "$GBFF" <<'EOF'
import glob, hashlib, json, sys
from pathlib import Path
from Bio import SeqIO
asm, gca, gcf, sample, fna, gbff = sys.argv[1:]
def seqs(acc):
    path = glob.glob(f"{asm}/{acc}/ncbi_dataset/data/{acc}/*_genomic.fna")[0]
    return {r.id: (len(r.seq), hashlib.sha256(str(r.seq).upper().encode()).hexdigest()) for r in SeqIO.parse(path, "fasta")}
a, b = seqs(gca), seqs(gcf)
same = sorted(v[1] for v in a.values()) == sorted(v[1] for v in b.values())
cds = sum(1 for r in SeqIO.parse(gbff, "genbank") for f in r.features if f.type == "CDS")
check = dict(sample=sample, gca=gca, gcf=gcf, gca_records=a, gcf_records=b,
             identical_sequence_sets=same, gcf_cds_features=cds,
             fasta=fna, fasta_sha256=hashlib.sha256(Path(fna).read_bytes()).hexdigest(),
             gbff=gbff, gbff_sha256=hashlib.sha256(Path(gbff).read_bytes()).hexdigest())
Path(asm, "assembly_check.json").write_text(json.dumps(check, indent=2) + "\n")
print(json.dumps({k: v for k, v in check.items() if "records" not in k}, indent=2))
if not same:
    sys.exit("GCA and GCF sequence sets differ; stop before callout")
manifest = dict(schema_version=1, samples=[dict(id=sample, fasta=fna, annotations=[gbff], annotation_format="genbank")])
Path(asm, "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
EOF
echo "[$(date -Is)] done"
