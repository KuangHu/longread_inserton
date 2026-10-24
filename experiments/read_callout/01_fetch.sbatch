#!/bin/bash
#SBATCH --job-name=rc01_fetch
#SBATCH --account=pc_rubinlab
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --array=1-4
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=8:00:00
#SBATCH --output=/global/scratch/users/kh36969/longread_callout/logs/%x_%A_%a.out
# Download and verify reads (ENA) and the GCA/GCF assemblies (NCBI) for one sample.
source /global/home/users/kh36969/middle_insertion_callout/longread_inserton/experiments/read_callout/common.sh
mkdir -p $ASM_DIR $OUT

# Reads: archive size/MD5, spot and base counts against ENA metadata.
if [ -f $READ_DIR/$RUN/verify.json ]; then
    echo "reads already verified: $READ_DIR/$RUN/verify.json"
else
    python $CODE/scripts/fetch_long_reads.py fetch --out $READ_DIR --runs $RUN
fi
FASTQ=$(ls $READ_DIR/$RUN/*.fastq.gz)
python $HERE/molecule_reads.py diagnose $FASTQ --limit 20000 > $OUT/read_unit.json

# Assemblies: GCA (submitted) and GCF (RefSeq, PGAP annotation used for the callout).
$DATASETS --version > $ASM_DIR/datasets_version.txt
for acc in $GCA $GCF; do
    if [ ! -d $ASM_DIR/$acc ]; then
        $DATASETS download genome accession $acc --include genome,gbff --filename $ASM_DIR/$acc.zip --no-progressbar
        unzip -q -o $ASM_DIR/$acc.zip -d $ASM_DIR/$acc
    fi
done
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
