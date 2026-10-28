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
# Paths & settings
############################################

BASE_PATH="/root/persistent"

OUTDIR="${BASE_PATH}/data/prs/SBayesRC/${PHENO}"
mkdir -p "$OUTDIR"

OUT_PREFIX="${OUTDIR}/sbayesrc_${PHENO}"

LD_DIR="${BASE_PATH}/opt/gctb_refs/LD_Reference/ukbEUR_Imputed"
ANNOT="${BASE_PATH}/opt/gctb_refs/annot_baseline2.2.txt"

THREADS=4
export OMP_NUM_THREADS=$THREADS

############################################
# Checks
############################################

[[ ! -f "$GWAS" ]] && { echo "ERROR: GWAS file not found: $GWAS"; exit 1; }
[[ ! -d "$LD_DIR" ]] && { echo "ERROR: LD reference not found: $LD_DIR"; exit 1; }
[[ ! -f "$ANNOT" ]] && { echo "ERROR: Annotation file not found: $ANNOT"; exit 1; }

############################################
# Run SBayesRC
############################################

echo "=========================================="
echo "SBayesRC | Phenotype: $PHENO"
echo "GWAS: $GWAS"
echo "=========================================="

# 1. Tidy
Rscript -e "
library(SBayesRC)
SBayesRC::tidy(
  mafile = '$GWAS',
  LDdir  = '$LD_DIR',
  output = '${OUT_PREFIX}_tidy.ma',
  log2file = TRUE
)
"

# 2. Impute
Rscript -e "
library(SBayesRC)
SBayesRC::impute(
  mafile = '${OUT_PREFIX}_tidy.ma',
  LDdir  = '$LD_DIR',
  output = '${OUT_PREFIX}_imp.ma',
  log2file = TRUE
)
"

# 3. Main SBayesRC
Rscript -e "
library(SBayesRC)
SBayesRC::sbayesrc(
  mafile   = '${OUT_PREFIX}_imp.ma',
  LDdir    = '$LD_DIR',
  outPrefix= '${OUT_PREFIX}_sbrc',
  annot    = '$ANNOT',
  log2file = TRUE
)
"

echo "SBayesRC completed for $PHENO"
