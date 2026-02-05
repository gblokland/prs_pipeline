#!/usr/bin/env bash
# ------------------------------------------------------------
# Wrapper script to run PRS pipeline for multiple phenotypes
# Author: Gabriëlla Blokland
# ------------------------------------------------------------

set -euo pipefail

############################################
# CONFIGURATION - edit only this part to specify your jobs and phenotypes
############################################

PIPELINE_SCRIPT="./run_prs_pipeline.sh"

# Which PRS methods to run:
#   ldpred2 | prscs | prsice2 | sbayesrc | all
METHOD="all"

# Base directory for GWAS summary statistics
SUMSTATS_BASE="/root/persistent/sumstats/Biomarkers"

# Phenotypes to process
PHENOTYPES=(
  CRP
  lymphocyte
  granulocyte
)

############################################
# Sanity checks
############################################

if [[ ! -x "${PIPELINE_SCRIPT}" ]]; then
  echo "ERROR: ${PIPELINE_SCRIPT} not found or not executable"
  exit 1
fi

############################################
# Run pipeline
############################################

echo "=========================================="
echo "PRS WRAPPER START"
echo "Methods   : ${METHOD}"
echo "Phenotypes: ${PHENOTYPES[*]}"
echo "=========================================="

for PHENO in "${PHENOTYPES[@]}"; do

  SUMSTATS="${SUMSTATS_BASE}/${PHENO}/${PHENO}_sumstats.QC.gz"

  if [[ ! -f "${SUMSTATS}" ]]; then
    echo "WARNING: Sumstats not found for ${PHENO}"
    echo "Expected: ${SUMSTATS}"
    echo "Skipping ${PHENO}"
    continue
  fi

  echo
  echo "------------------------------------------"
  echo "Running PRS for phenotype: ${PHENO}"
  echo "------------------------------------------"

  bash "${PIPELINE_SCRIPT}" \
    --phenotype "${PHENO}" \
    --sumstats "${SUMSTATS}" \
    --method "${METHOD}"

done

echo
echo "=========================================="
echo "PRS WRAPPER FINISHED"
echo "=========================================="
