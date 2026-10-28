#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(SBayesRC)
  library(data.table)
})

############################################################
# Arguments
############################################################

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop(
    "Usage:\n",
    "  Rscript prs_scoring_sbayesrc.R <PHENOTYPE> <GWAS_SUMSTATS_COJO> <OUTDIR>\n\n",
    "Example:\n",
    "  Rscript prs_scoring_sbayesrc.R CRP CRP_sumstats_COJO.tsv /root/persistent/data/prs/SBayesRC"
  )
}

PHENO        <- args[1]
GWAS_COJO   <- args[2]
OUTDIR      <- args[3]

############################################################
# Fixed paths (pipeline-level configuration)
############################################################

BASE_PATH <- "/root/persistent"

LD_FOLDER <- file.path(
  BASE_PATH,
  "opt/gctb_refs/LD_Reference/ukbEUR_Imputed"
)

ANNOT_FILE <- file.path(
  BASE_PATH,
  "opt/gctb_refs/annot_baseline2.2.txt"
)

GENO_PREFIX <- file.path(
  BASE_PATH,
  "data/imp/dmsa_imp_eur_1000GP_P3_auto_xchr",
  "cobg_dir_genome_wide",
  "t2d_dmsa_eur_gb-qc1.hg19.ch.fl.bgn.reid.randomid.filtered_geno_hwe_mind"
)

############################################################
# Output structure
############################################################

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

OUT_PREFIX <- file.path(
  OUTDIR,
  paste0("sbayesrc_", PHENO)
)

############################################################
# Step 1: Tidy GWAS summary statistics
############################################################

TIDY_MA <- paste0(OUT_PREFIX, "_tidy.ma")

if (!file.exists(TIDY_MA)) {
  message("▶ Tidying GWAS summary statistics for ", PHENO)
  SBayesRC::tidy(
    mafile  = GWAS_COJO,
    LDdir   = LD_FOLDER,
    output  = TIDY_MA,
    log2file = TRUE
  )
} else {
  message("✔ Tidy GWAS file exists, skipping")
}

############################################################
# Step 2: (Optional) Imputation
# Uncomment if your GWAS does NOT fully overlap LD panel
############################################################

# IMP_MA <- paste0(OUT_PREFIX, "_imp.ma")
# if (!file.exists(IMP_MA)) {
#   message("▶ Imputing missing SNPs")
#   SBayesRC::impute(
#     mafile  = TIDY_MA,
#     LDdir   = LD_FOLDER,
#     output  = IMP_MA,
#     log2file = TRUE
#   )
# }
# MA_INPUT <- IMP_MA

MA_INPUT <- TIDY_MA

############################################################
# Step 3: Run SBayesRC
############################################################

SBR_OUT_PREFIX <- paste0(OUT_PREFIX, "_model")

if (!file.exists(paste0(SBR_OUT_PREFIX, ".txt"))) {
  message("▶ Running SBayesRC for ", PHENO)
  SBayesRC::sbayesrc(
    mafile   = MA_INPUT,
    LDdir    = LD_FOLDER,
    outPrefix = SBR_OUT_PREFIX,
    annot    = ANNOT_FILE,
    log2file = TRUE
  )
} else {
  message("✔ SBayesRC output exists, skipping")
}

############################################################
# Step 4: PRS calculation
############################################################

PRS_OUT_PREFIX <- paste0(OUT_PREFIX, "_prs")

if (!file.exists(paste0(PRS_OUT_PREFIX, ".score.txt"))) {
  message("▶ Calculating PRS for ", PHENO)
  SBayesRC::prs(
    weight     = paste0(SBR_OUT_PREFIX, ".txt"),
    genoPrefix = GENO_PREFIX,
    genoCHR    = "",
    outPrefix  = PRS_OUT_PREFIX
  )
} else {
  message("✔ PRS score exists, skipping")
}

############################################################
# Done
############################################################

message("✅ SBayesRC pipeline completed for phenotype: ", PHENO)
