#!/usr/bin/env Rscript
set.seed(12345)

library(bigsnpr)
library(data.table)
library(magrittr)

############################################
# Command line arguments
############################################

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript prs_scoring_ldpred2.R <PHENOTYPE> <GWAS_SUMSTATS>")
}

PHENO <- args[1]
GWAS <- args[2]

############################################
# Base paths
############################################

BASE_PATH <- "/root/persistent"

# LDpred2 file paths
SUM_STATS_FILE <- paste0(BASE_PATH, "/sumstats/Biomarkers/", PHENO, "/", PHENO, "_sumstats_QC.gz")
GENO_PREFIX <- paste0(BASE_PATH, "/data/imp/dmsa_imp_eur_1000GP_P3_auto_xchr/cobg_dir_genome_wide/t2d_dmsa_eur_gb-qc1.hg19.ch.fl.bgn.reid.randomid")

# Output directory for PRS scores
OUTDIR <- paste0(BASE_PATH, "/data/prs/LDpred2/", PHENO)
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# LD reference
LD_REF_DIR <- paste0(BASE_PATH, "/opt/gctb_refs/LD_Reference/ukbEUR_Imputed")

# PCA & Covariates
PHENO_FILE <- paste0(BASE_PATH, "/data/prs/UKB_pheno.txt")
COV_FILE <- paste0(BASE_PATH, "/data/prs/UKB_cov.txt")

############################################
# Read input data
############################################

# Read phenotype, covariate, and PC data
phenotype <- fread(PHENO_FILE)
covariate <- fread(COV_FILE)
pcs <- fread("EUR.eigenvec")  # Assuming this file exists and is correctly formatted

# Merge phenotype, covariate, and PCs
colnames(pcs) <- c("FID", "IID", paste0("PC", 1:6))
pheno_data <- merge(phenotype, covariate) %>%
  merge(pcs, by = c("FID", "IID"))

# Read summary statistics file and format it for LDpred
sumstats <- bigreadr::fread2(SUM_STATS_FILE)
names(sumstats) <- c("chr", "pos", "rsid", "a1", "a0", "n_eff", "beta_se", "p", "OR", "INFO", "MAF")
sumstats$beta <- log(sumstats$OR)
sumstats <- sumstats[sumstats$rsid %in% info$rsid, ]

# Read the HapMap3 SNP information
info <- readRDS(runonce::download_file("https://ndownloader.figshare.com/files/25503788", fname = "map_hm3_ldpred2.rds"))

############################################
# Calculate the LD matrix (using SNPs from the HapMap3 reference)
############################################

NCORES <- nb_cores()  # Get the maximum number of cores available
tmp <- tempfile(tmpdir = "tmp-data")
on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE)

snp_readBed("EUR.QC.bed")
obj.bigSNP <- snp_attach("EUR.QC.rds")
map <- obj.bigSNP$map[-3]
names(map) <- c("chr", "rsid", "pos", "a1", "a0")

# Match SNPs in the genotype and summary statistics
info_snp <- snp_match(sumstats, map)

genotype <- obj.bigSNP$genotypes
CHR <- map$chr
POS <- map$pos

POS2 <- snp_asGeneticPos(CHR, POS, dir = ".")

# Calculate LD for each chromosome
ld <- NULL
corr <- NULL

for (chr in 1:22) {
    ind.chr <- which(info_snp$chr == chr)
    ind.chr2 <- info_snp$`_NUM_ID_`[ind.chr]
    corr0 <- snp_cor(
        genotype,
        ind.col = ind.chr2,
        ncores = NCORES,
        infos.pos = POS2[ind.chr2],
        size = 3 / 1000
    )
    if (chr == 1) {
        ld <- Matrix::colSums(corr0^2)
        corr <- as_SFBM(corr0, tmp)
    } else {
        ld <- c(ld, Matrix::colSums(corr0^2))
        corr$add_columns(corr0, nrow(corr))
    }
}

# Store fam order
fam.order <- as.data.table(obj.bigSNP$fam)
setnames(fam.order, c("family.ID", "sample.ID"), c("FID", "IID"))

############################################
# LD Score Regression to estimate heritability
############################################

df_beta <- info_snp[, c("beta", "beta_se", "n_eff", "_NUM_ID_")]
ldsc <- snp_ldsc(ld, length(ld), chi2 = (df_beta$beta / df_beta$beta_se)^2,
                 sample_size = df_beta$n_eff, blocks = NULL)

h2_est <- ldsc[["h2"]]

############################################
# Calculate PRS with LDpred2
############################################

# LDpred2: Infinitesimal model
beta_inf <- snp_ldpred2_inf(corr, df_beta, h2 = h2_est)

# Calculate PRS for all samples
genotype <- obj.bigSNP$genotypes
ind.test <- 1:nrow(genotype)
pred_inf <- big_prodVec(genotype, beta_inf, ind.row = ind.test, ind.col = info_snp$`_NUM_ID_`)

# Sum the PRS across chromosomes
for (chr in 1:22) {
    obj.bigSNP <- snp_attach(paste0("EUR_chr", chr, ".rds"))
    genotype <- obj.bigSNP$genotypes
    ind.test <- 1:nrow(genotype)
    chr.idx <- which(info_snp$chr == chr)
    ind.chr <- info_snp$`_NUM_ID_`[chr.idx]
    tmp <- big_prodVec(genotype, beta_inf, ind.row = ind.test, ind.col = ind.chr)
    pred_inf <- pred_inf + tmp
}

############################################
# Final performance evaluation
############################################

# Reformat phenotype data
y <- pheno_data[fam.order, on = c("FID", "IID")]

# Linear regression for quantitative traits
reg.formula <- paste("PC", 1:6, sep = "", collapse = "+") %>%
  paste0("Height~PRS+Sex+", .) %>%
  as.formula

reg.dat <- y
reg.dat$PRS <- pred_inf

inf.model <- lm(reg.formula, data = reg.dat) %>%
  summary

result <- data.table(
  infinitesimal = inf.model$r.squared - null.r2,
  null = null.r2
)

# Save results
write.table(result, file = paste0(OUTDIR, "/LDpred2_results_", PHENO, ".txt"), sep = "\t", row.names = FALSE)

cat("LDpred2 PRS calculation completed for", PHENO, "\n")
