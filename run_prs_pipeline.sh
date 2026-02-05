#!/usr/bin/env bash
set -euo pipefail

############################################
# Usage & argument parsing
############################################

usage() {
  cat <<EOF
Usage:
  ./run_prs_pipeline.sh \
    -p <phenotype> \
    -g <gwas_sumstats_path> \
    -m <method>

Required arguments:
  -p   Phenotype name (e.g. MDD, SCZ)
  -g   Path to GWAS summary statistics
  -m   Method to run:
         ldpred2
         prscs
         prsice2
         sbayesrc
         all

Example:
  ./run_prs_pipeline.sh -p MDD -g gwas/MDD.sumstats.gz -m all
EOF
  exit 1
}

while getopts ":p:g:m:h" opt; do
  case $opt in
    p) PHENO="$OPTARG" ;;
    g) GWAS="$OPTARG" ;;
    m) METHOD="$OPTARG" ;;
    h) usage ;;
    \?) echo "Invalid option: -$OPTARG" >&2; usage ;;
    :) echo "Option -$OPTARG requires an argument." >&2; usage ;;
  esac
done

############################################
# Input checks
############################################

[[ -z "${PHENO:-}" ]] && usage
[[ -z "${GWAS:-}" ]] && usage
[[ -z "${METHOD:-}" ]] && usage

if [[ ! -f "$GWAS" ]]; then
  echo "ERROR: GWAS file not found: $GWAS"
  exit 1
fi

############################################
# Script locations
############################################

LDPRED2_SCRIPT="prs_scoring_ldpred2.R"
PRSCS_SCRIPT="prs_scoring_prscs.sh"
PRSICE2_SCRIPT="prs_scoring_prsice2.sh"
SBAYESRC_R_SCRIPT="prs_scoring_sbayesrc.R"
SBAYESRC_SH_SCRIPT="prs_scoring_sbayesrc.sh"

############################################
# Helper function
############################################

run_step() {
  echo
  echo "=========================================="
  echo "Running $1 for phenotype $PHENO"
  echo "=========================================="
}

############################################
# Method dispatch
############################################

case "$METHOD" in
  ldpred2)
    run_step "LDpred2"
    Rscript "$LDPRED2_SCRIPT" "$PHENO" "$GWAS"
    ;;

  prscs)
    run_step "PRS-CS"
    bash "$PRSCS_SCRIPT" "$PHENO" "$GWAS"
    ;;

  prsice2)
    run_step "PRSice-2"
    bash "$PRSICE2_SCRIPT" "$PHENO" "$GWAS"
    ;;

  sbayesrc)
    run_step "SBayesRC"
    Rscript "$SBAYESRC_R_SCRIPT" "$PHENO" "$GWAS"
    bash "$SBAYESRC_SH_SCRIPT" "$PHENO" "$GWAS"
    ;;

  all)
    run_step "LDpred2"
    Rscript "$LDPRED2_SCRIPT" "$PHENO" "$GWAS"

    run_step "PRS-CS"
    bash "$PRSCS_SCRIPT" "$PHENO" "$GWAS"

    run_step "PRSice-2"
    bash "$PRSICE2_SCRIPT" "$PHENO" "$GWAS"

    run_step "SBayesRC"
    Rscript "$SBAYESRC_R_SCRIPT" "$PHENO" "$GWAS"
    bash "$SBAYESRC_SH_SCRIPT" "$PHENO" "$GWAS"
    ;;

  *)
    echo "ERROR: Unknown method '$METHOD'"
    usage
    ;;
esac

echo
echo "PRS pipeline finished successfully for phenotype: $PHENO"

