##This is script to generate cellrange counts for RNA and ADT assay
#Creating a feature_ref.csv file is must
#Total_SEqC panel for ADTs used to generate feature_ref.csv is added as well

#!/bin/bash -l
#SBATCH --job-name=MDS_CRmulti
#SBATCH --output=MDS_CRmulti_%A_%a.out
#SBATCH --error=MDS_CRmulti_%A_%a.err
#SBATCH --time=7-00:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=500G
#SBATCH --partition=defq
#SBATCH --array=1-18%4

module purge || true
module load cellranger/9.0.1
set -eox pipefail

# INPUTS 
GEX_DIR="/net/beegfs/users/P089417/MDS_Data/MDS_GEx/data/FR34044381_10x-GEX-Library-MDS-submission-form-1/B22W323LT3"
ADT_DIR="/net/beegfs/users/P089417/MDS_Data/MDS_prot/data/FR34044382_10x-Protein-Library-MDS-submission-form-2/B22W323LT3"
REF="/net/beegfs/users/P089417/Refs/refdata-gex-GRCh38-2024-A"
FEATURE_REF="/net/beegfs/users/P089417/MDS_CRfiles/feature_ref.csv"   # id,name,read,pattern,sequence,feature_type
OUT_BASE="/net/beegfs/users/P089417/MDS_OUTs"
CHEM="SC5P-R2"

SAMPLES=(MDS001 MDS005 MDS006 MDS010 MDS016 MDS023 MDS029 MDS038 MDS059 MDS065 MDS154 MDS155 MDS167 MDS169 MDS180 MDS189 MDS201 MDS212)
S="${SAMPLES[$((SLURM_ARRAY_TASK_ID-1))]}"

echo "[INFO] $(date) Host=$(hostname) User=$(whoami)"
echo "[INFO] Cell Ranger: $(which cellranger)"; cellranger --version

test -d "$GEX_DIR"; test -d "$ADT_DIR"; test -f "$FEATURE_REF"

# Confirm feature_ref header (Cell Ranger 9 expects exact order)
FR_HDR="$(head -n1 "$FEATURE_REF" | tr -d '\r')"
if [[ "$FR_HDR" != "id,name,read,pattern,sequence,feature_type" ]]; then
  echo "[FATAL] feature_ref.csv header '$FR_HDR' != 'id,name,read,pattern,sequence,feature_type'" >&2
  exit 10
fi

#Per-sample dirs & links
ROOT="${OUT_BASE}/${S}"
GEX_OUT="${ROOT}/GEX"
ADT_OUT="${ROOT}/ADT"
mkdir -p "$GEX_OUT" "$ADT_OUT"

( cd "$GEX_OUT" && ln -sf ${GEX_DIR}/${S}-*_R*.fastq.gz . || true )
( cd "$ADT_OUT" && ln -sf ${ADT_DIR}/${S}-*_R*.fastq.gz . || true )

ls -lh "$GEX_OUT" | sed -n '1,20p' || true
ls -lh "$ADT_OUT" | sed -n '1,20p' || true

# Require at least one R1 file in each
n_gex=$(ls -1 "$GEX_OUT"/*_R1_*.fastq.gz 2>/dev/null | wc -l || true)
n_adt=$(ls -1 "$ADT_OUT"/*_R1_*.fastq.gz 2>/dev/null | wc -l || true)
if [[ "$n_gex" -eq 0 ]]; then echo "[FATAL] No GEX FASTQs for ${S}" >&2; exit 3; fi
if [[ "$n_adt" -eq 0 ]]; then echo "[FATAL] No ADT FASTQs for ${S}" >&2; exit 3; fi

# AUTO-DETECT fastq_id from filenames (prefix before `_S`)
# Example: MDS001-09-203_S17_L001_R1_001.fastq.gz -> fastq_id = MDS001-09-203
gex_first=$(ls "$GEX_OUT"/*_R1_*.fastq.gz | head -n1)
adt_first=$(ls "$ADT_OUT"/*_R1_*.fastq.gz | head -n1)
GEX_ID=$(basename "$gex_first"); GEX_ID="${GEX_ID%%_S*}"
ADT_ID=$(basename "$adt_first"); ADT_ID="${ADT_ID%%_S*}"

echo "[INFO] Detected fastq_id (GEX): $GEX_ID"
echo "[INFO] Detected fastq_id (ADT): $ADT_ID"

# --- Write v9-compliant multi_config.csv ---
CFG="${ROOT}/multi_config.csv"
cat > "$CFG" <<EOF
[gene-expression]
reference,${REF}
chemistry,${CHEM}
create-bam,true

[feature]
reference,${FEATURE_REF}

[libraries]
fastq_id,fastqs,feature_types
${GEX_ID},${GEX_OUT},Gene Expression
${ADT_ID},${ADT_OUT},Antibody Capture
EOF

echo "[INFO] multi_config.csv:"
head -n 20 "$CFG"

# Skip if finished already
if [[ -d "${ROOT}/${S}_multi/outs" ]]; then
  echo "[INFO] Outputs exist for ${S}; skipping."
  exit 0
fi

# --- Run ---
cd "$ROOT"
cellranger multi \
  --id="${S}_multi" \
  --csv="$CFG" \
  --localcores="${SLURM_CPUS_PER_TASK}" \
  --localmem="$((SLURM_MEM_PER_NODE/1024))" \
  | tee "${S}_multi.stdout"

# Quick QC print
MS="${ROOT}/${S}_multi/outs/metrics_summary.csv"
if [[ -f "$MS" ]]; then
  grep -i 'Antibody: Median UMI counts per cell' "$MS" || true
  grep -i 'Feature Reads Mapped to Features' "$MS" || true
fi

echo "[INFO] $(date) Done ${S}"
