#!/usr/bin/env bash
# ------------------------------------------------------------
# SBayesRC LD step 3 – eigen decomposition per LD block
# Pipeline-compatible version
# Author: Gabriëlla Blokland, Maastricht University
# ------------------------------------------------------------

set -euo pipefail

############################################
# Arguments
############################################

if [[ $# -lt 2 ]]; then
  echo "Usage:"
  echo "  prs_SBayesRC_LDstep3.sh <PHENOTYPE> <BLOCK_START> [BLOCK_END]"
  echo
  echo "Example:"
  echo "  prs_SBayesRC_LDstep3.sh CRP 1 591"
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
# Environment
############################################

# Control threading for eigen decomposition
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-4}"

############################################
# Sanity checks
############################################

if [[ ! -d "${OUTDIR}" ]]; then
  echo "ERROR: LD directory not found:"
  echo "  ${OUTDIR}"
  echo "Run prepare_sbayesrc_ld.sh and LDstep2 first."
  exit 1
fi

############################################
# Run LDstep3
############################################

echo "=========================================="
echo "SBayesRC LDstep3"
echo "Phenotype   : ${PHENO}"
echo "Block range : ${BLOCK_START}–${BLOCK_END}"
echo "Threads     : ${OMP_NUM_THREADS}"
echo "Output dir  : ${OUTDIR}"
echo "=========================================="

for (( IDX=${BLOCK_START}; IDX<=${BLOCK_END}; IDX++ )); do

  EIGEN_FILE="${OUTDIR}/block${IDX}.eigen.bin"
  LOG_FILE="${OUTDIR}/LDstep3_block${IDX}.log"

  if [[ -f "${EIGEN_FILE}" ]]; then
    echo "✔ Block ${IDX} eigen decomposition exists — skipping"
    continue
  fi

  echo "▶ Running LDstep3 for block ${IDX}"

  Rscript -e "
    suppressPackageStartupMessages(library(SBayesRC))
