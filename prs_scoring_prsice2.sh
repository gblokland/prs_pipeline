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
# Base paths (adjust once, not per phenotype)
############################################

BASE_PATH="/root/persistent"

# Target genotype (PRS target cohort)
TARGET_PREFIX="${BASE_PATH}/data/UKB500k_241121/UKB500k_chr#_241121_Qced"

# Covariates & phenotype files
PHENO_FILE="${BASE_PATH}/data/prs/UKB_pheno.txt"
COV_FILE="${BASE_PATH}/data/prs/UKB_cov.txt"

# PRSice installation
OPT_DIR="${BASE_PATH}/opt"
PRSICE_BIN="${OPT_DIR}/PRSice_linux"
PRSICE_R="${OPT_DIR}/PRSice.R"
R_LIB="/usr/lib/R/library"

# Output
OUTDIR="${BASE_PATH}/data/prs/PRSice2/${PHENO}"
mkdir -p "$OUTDIR"
cd "$OUTDIR"

############################################
# Checks
############################################

[[ ! -f "$GWAS" ]] && { echo "ERROR: GWAS not found: $GWAS"; exit 1; }
[[ ! -f "$PRSICE_BIN" ]] && { echo "ERROR: PRSice binary not found"; exit 1; }
[[ ! -f "$PRSICE_R" ]] && { echo "ERROR: PRSice.R not found"; exit 1; }

############################################
# Trait type inference (simple heuristic)
############################################
# You can override this later with a flag

BINARY_TARGET="F"
STAT="BETA"
BETA_FLAG="--beta"

if [[ "$PHENO" =~ ^(SZ|MDD|BD|ADHD|ASD|CAD|T2D|ANX|AD).* ]]; then
  BINARY_TARGET="T"
  STAT="OR"
  BETA_FLAG="--or"
fi

############################################
# Common PRSice options
############################################

COMMON_OPTS=(
  --dir "$R_LIB"
  --prsice "$PRSICE_BIN"
  --base "$GWAS"
  --target "$TARGET_PREFIX"
  --binary-target "$BINARY_TARGET"
  --stat "$STAT"
  $BETA_FLAG
  --clump-kb 500kb
  --clump-r2 0.2
  --clump-p 1
  --base-maf MAF:0.01
  --missing MEAN_IMPUTE
  --pheno "$PHENO_FILE"
  --cov "$COV_FILE"
  --bar-levels 5e-8,5e-7,5e-6,5e-5,5e-4,5e-3,5e-2,0.1,0.2,0.3,0.4,0.5
  --quantile 20
  --quant-break 1,5,10,15,20
  --score sum
  --all-score
  --out "${OUTDIR}/${PHENO}"
)

############################################
# Run PRSice (1st pass: create .valid)
############################################

echo "=========================================="
echo "PRSice2 | Phenotype: $PHENO"
echo "GWAS: $GWAS"
echo "Binary trait: $BINARY_TARGET"
echo "=========================================="

Rscript "$PRSICE_R" "${COMMON_OPTS[@]}"

############################################
# Run PRSice (2nd pass: extract valid SNPs)
############################################

VALID_FILE="${OUTDIR}/${PHENO}.valid"

if [[ -f "$VALID_FILE" ]]; then
  echo "Found .valid file — rerunning with --extract"
  Rscript "$PRSICE_R" \
    "${COMMON_OPTS[@]}" \
    --extract "$VALID_FILE"
else
  echo "WARNING: .valid file not produced — skipping extract step"
fi

echo "PRSice2 completed for $PHENO"
