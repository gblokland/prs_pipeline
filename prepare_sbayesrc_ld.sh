#!/usr/bin/env bash
set -euo pipefail

############################################
# Purpose:
#   Prepare LD reference for SBayesRC
#   (LDstep1–LDstep4)
#
# Run ONCE per ancestry + SNP panel
############################################

############################################
# User-configurable paths
############################################

BASE_PATH="/root/persistent"

# Genotype data used as LD reference (PLINK prefix)
GENO_PREFIX="${BASE_PATH}/data/imp/dmsa_imp_eur_1000GP_P3_auto_xchr/cobg_dir_genome_wide/t2d_dmsa_eur_gb-qc1.hg19.ch.fl.bgn.reid.randomid.filtered_geno_hwe_mind_maf01"

# A representative GWAS summary file in COJO format
# (only used to determine SNP list & alleles)
GWAS_EXAMPLE="${BASE_PATH}/sumstats/Biomarkers/CRP/CRP_sumstats_COJO.tsv"

# Output directory for LD
OUTDIR="${BASE_PATH}/data/prs/SBayesRC/ld_ukbEUR_imputed"

THREADS=4
export OMP_NUM_THREADS=$THREADS

############################################
# Optional parameters
############################################

GENO_CHR=""        # e.g. "1-22,X" if split by chr
BLOCK_REF=""       # leave empty to use default GRCh37 blocks
MAX_JOBS=25        # parallel LDstep2 jobs

############################################
# Sanity checks
############################################

[[ ! -f "$GWAS_EXAMPLE" ]] && { echo "ERROR: GWAS example not found"; exit 1; }
[[ ! -f "${GENO_PREFIX}.bed" ]] && { echo "ERROR: PLINK bed not found"; exit 1; }

mkdir -p "$OUTDIR"

############################################
# Step 1: LDstep1 (block definition)
############################################

if [[ ! -f "${OUTDIR}/ldm.info" ]]; then
  echo "Running LDstep1..."
  Rscript -e "
  library(SBayesRC)
  SBayesRC::LDstep1(
    mafile    = '$GWAS_EXAMPLE',
    genoPrefix= '$GENO_PREFIX',
    outDir    = '$OUTDIR',
    genoCHR   = '$GENO_CHR',
    blockRef  = '$BLOCK_REF',
    log2file  = TRUE
  )
  "
else
  echo "LDstep1 already completed — skipping"
fi

############################################
# Step 2: LDstep2 (LD matrices per block)
############################################

NUM_BLOCKS=$(ls "$OUTDIR/snplist"/*.snplist | wc -l)

echo "Detected $NUM_BLOCKS LD blocks"

job_pids=()

check_jobs() {
  while [[ ${#job_pids[@]} -ge $MAX_JOBS ]]; do
    wait -n
    job_pids=($(jobs -rp))
  done
}

for idx in $(seq 1 "$NUM_BLOCKS"); do
  if [[ ! -f "${OUTDIR}/b${idx}.ldm.full.bin" ]]; then
    echo "Submitting LDstep2 block $idx"
    check_jobs
    Rscript -e "
      library(SBayesRC)
      SBayesRC::LDstep2(
        outDir     = '$OUTDIR',
        blockIndex = $idx,
        log2file   = TRUE
      )
    " &> "${OUTDIR}/LDstep2_${idx}.log" &
    job_pids+=($!)
  fi
done

wait
echo "LDstep2 completed"

############################################
# Step 3: LDstep3 (eigen decomposition)
############################################

for idx in $(seq 1 "$NUM_BLOCKS"); do
  if [[ ! -f "${OUTDIR}/block${idx}.eigen.bin" ]]; then
    echo "Running LDstep3 block $idx"
    Rscript -e "
      library(SBayesRC)
      SBayesRC::LDstep3(
        outDir     = '$OUTDIR',
        blockIndex = $idx,
        log2file   = TRUE
      )
    " &> "${OUTDIR}/LDstep3_${idx}.log"
  fi
done

############################################
# Step 4: LDstep4 (merge LD info)
############################################

if [[ ! -f "${OUTDIR}/snp.info" ]]; then
  echo "Running LDstep4..."
  Rscript -e "
  library(SBayesRC)
  SBayesRC::LDstep4(
    outDir   = '$OUTDIR',
    log2file = TRUE
  )
  "
else
  echo "LDstep4 already completed — skipping"
fi

############################################
# Done
############################################

echo "=========================================="
echo "SBayesRC LD preparation completed"
echo "LD directory: $OUTDIR"
echo "=========================================="
