# prs_pipeline

Edit the prs_pipeline_wrapper_script.sh to specify your training GWAS sumstats file and other specifics.

Run the script like this:

```
./prs_pipeline_wrapper_script.sh
```

When you run this wrapper script it will call the run_prs_pipeline.sh script like this:

Usage: ./run_prs_pipeline.sh
-p <phenotype_name>
-g <gwas_sumstats_path>
-m <prs scoring method(s)>

Required arguments: -p Phenotype name (e.g. MDD, SCZ) -g Path to GWAS summary statistics -m Method to run: ldpred2 prscs prsice2 sbayesrc all

Example:
```
./run_prs_pipeline.sh -p MDD -g gwas/MDD.sumstats.gz -m all
```
