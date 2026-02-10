# prs_pipeline

Usage: ./run_prs_pipeline.sh
-p
-g <gwas_sumstats_path>
-m

Required arguments: -p Phenotype name (e.g. MDD, SCZ) -g Path to GWAS summary statistics -m Method to run: ldpred2 prscs prsice2 sbayesrc all

Example: ./run_prs_pipeline.sh -p MDD -g gwas/MDD.sumstats.gz -m all
