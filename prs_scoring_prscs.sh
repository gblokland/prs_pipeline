#!/usr/bin/env bash
set -euo pipefail

############################################
# Arguments
############################################

if [[ $# -lt 2 ]]; then
  echo "Usage: $0 <PHENOTYPE> <GWAS_SUMSTATS>"
  exit 1
fi

PHENO="$1"
GWAS="$2"

############################################
# Base paths
############################################

BASE_PATH="/root/persistent"

# PRS-CS installation
PRSCS_DIR="${BASE_PATH}/opt/PRScs"
PRSCS_PY="${PRSCS_DIR}/PRScs.py"

# LD reference (match GWAS ancestry!)
REF_DIR="${PRSCS_DIR}/Reference/1KG/ldblk_1kg_eur"

# Target genotype (PLINK prefix, no extension)
TARGET_PREFIX="${BASE_PATH}/data/imp/dmsa_imp_eur_1000GP_P3_auto_xchr/\
cobg_dir_genome_wide/t2d_dmsa_eur_gb-qc1.hg19.ch.fl.bgn.reid.randomid.filtered_geno_hwe_mind"

# Output
OUTDIR="${BASE_PATH}/data/prs/PRScs/${PHENO}"
mkdir -p "$OUTDIR"

############################################
# PRS-CS parameters
############################################

A=1
B=0.5
PHI=1e-4          # fixed; tune outside pipeline if needed
N_ITER=1000
N_BURNIN=500
THIN=5
SEED=12345
CHROM_LIST=$(seq -s, 1 22)

############################################
# GWAS sample size (required by PRS-CS)
############################################
# Adjust or externalise later

case "$PHENO" in
  CRP)          N_GWAS=575531 ;;
  lymphocyte)  N_GWAS=173480 ;;
  granulocyte) N_GWAS=173480 ;;
  *)
    echo "ERROR: GWAS sample size not defined for phenotype '$PHENO'"
    exit 1
    ;;
esac

############################################
# Checks
############################################

[[ ! -f "$GWAS" ]] && { echo "ERROR: GWAS file not found"; exit 1; }
[[ ! -f "$PRSCS_PY" ]] && { echo "ERROR: PRScs.py not found"; exit 1; }
[[ ! -d "$REF_DIR" ]] && { echo "ERROR: LD reference not found"; exit 1; }

############################################
# Run PRS-CS
############################################

echo "=========================================="
echo "PRS-CS | Phenotype: $PHENO"
echo "GWAS: $GWAS"
echo "N_GWAS: $N_GWAS"
echo "PHI: $PHI"
echo "=========================================="

python "$PRSCS_PY" \
  --ref_dir "$REF_DIR" \
  --bim_prefix "$TARGET_PREFIX" \
  --sst_file "$GWAS" \
  --n_gwas "$N_GWAS" \
  --a "$A" \
  --b "$B" \
  --phi "$PHI" \
  --n_iter "$N_ITER" \
  --n_burnin "$N_BURNIN" \
  --thin "$THIN" \
  --chrom "$CHROM_LIST" \
  --out_dir "$OUTDIR" \
  --seed "$SEED"

############################################
# Merge chromosome-specific effect sizes
############################################

MERGED="${OUTDIR}/${PHENO}_pst_eff_a${A}_b${B}_phi${PHI}_merged.txt"

cat "${OUTDIR}/${PHENO}_pst_eff_a${A}_b${B}_phi${PHI}_chr"*.txt > "$MERGED"

echo "Merged effect sizes: $MERGED"

############################################
# Calculate PRS using PLINK2
############################################

plink2 \
  --bfile "$TARGET_PREFIX" \
  --score "$MERGED" 2 4 6 \
  --out "${OUTDIR}/${PHENO}_PRSCS_scores"

echo "PRS-CS completed for $PHENO"
