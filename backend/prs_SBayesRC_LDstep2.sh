#!/usr/bin/env bash
# ------------------------------------------------------------
# SBayesRC LD step 2 – generate LD matrices per block
# Pipeline-compatible version
# Author: Gabriëlla Blokland, Maastricht University
# ------------------------------------------------------------

set -euo pipefail

############################################
# Arguments
############################################

if [[ $# -lt 2 ]]; then
  echo "Usage:"
  echo "  prs_SBayesRC_LDstep2.sh <PHENOTYPE> <BLOCK_START> [BLOCK_END]"
  echo
  echo "Example:"
  echo "  prs_SBayesRC_LDstep2.sh CRP 1 591"
  exit 1
fi

PHENO="$1"
BLOCK_START="$2"
BLOCK_END="${3:-$BLOCK_START}"

############################################
# Fixed paths (pipeline-level config)
############################################

BASE_PATH="/root/persistent"

OUTDIR="${BASE_PATH}/data/prs/SBayesRC/sbayesrc_ld_${PHENO}"

############################################
# Sanity checks
############################################

if [[ ! -d "${OUTDIR}" ]]; then
  echo "ERROR: LD directory not found:"
  echo "  ${OUTDIR}"
  echo "Run prepare_sbayesrc_ld.sh first."
  exit 1
fi

############################################
# Run LDstep2
############################################

echo "=========================================="
echo "SBayesRC LDstep2"
echo "Phenotype   : ${PHENO}"
echo "Block range : ${BLOCK_START}–${BLOCK_END}"
echo "Output dir  : ${OUTDIR}"
echo "=========================================="

for (( IDX=${BLOCK_START}; IDX<=${BLOCK_END}; IDX++ )); do

  BIN_FILE="${OUTDIR}/b${IDX}.ldm.full.bin"
  LOG_FILE="${OUTDIR}/LDstep2_block${IDX}.log"

  if [[ -f "${BIN_FILE}" ]]; then
    echo "✔ Block ${IDX} already completed — skipping"
    continue
  fi

  echo "▶ Running LDstep2 for block ${IDX}"

  Rscript -e "
    suppressPackageStartupMessages(library(SBayesRC))
    SBayesRC::LDstep2(
      outDir='${OUTDIR}',
      blockIndex=${IDX},
      log2file=TRUE
    )
  " > "${LOG_FILE}" 2>&1

done

echo "=========================================="
echo "LDstep2 completed for ${PHENO}"
echo "=========================================="
