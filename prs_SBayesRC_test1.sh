#!/bin/bash

MAX_JOBS=40
job_pids=()  # Array to store PIDs of background jobs

# Function to reap completed jobs
reap_jobs() {
    for pid in "${job_pids[@]}"; do
        if ! kill -0 $pid 2>/dev/null; then
            # Remove finished PID from the array
            wait $pid 2>/dev/null
            echo "Reaped job with PID $pid."
            job_pids=("${job_pids[@]/$pid}")
        fi
    done
}

# Trap SIGCHLD to reap child processes as they exit
trap 'reap_jobs' SIGCHLD

# Function to check if job slots are available
check_jobs() {
    echo "Checking running jobs (${#job_pids[@]} currently running)..."
    while [ ${#job_pids[@]} -ge $MAX_JOBS ]; do
        echo "Maximum number of jobs (${MAX_JOBS}) reached. Waiting for a job to finish..."
        reap_jobs
        sleep 1
    done
    echo "Slots available for new jobs. Proceeding..."
}

base_path="/root/persistent/"

for pheno in CRP granulocyte lymphocyte; do
    ma_file="${base_path}sumstats/Biomarkers/${pheno}/${pheno}_sumstats_COJO.tsv"
    genoPrefix="${base_path}data/imp/dmsa_imp_eur_1000GP_P3_auto_xchr/cobg_dir_genome_wide/t2d_dmsa_eur_gb-qc1.hg19.ch.fl.bgn.reid.randomid.filtered_geno_hwe_mind"
    outDir="${base_path}data/prs/SBayesRC/sbayesrc_ld_${pheno}"
    
    # Ensure output directory exists
    mkdir -p "$outDir"

    threads=4  # Number of CPU cores

    for idx in {1..591}; do
        if [ ! -f "$outDir/b$idx.ldm.full.bin" ]; then
            echo "Preparing to submit job $idx..."
            check_jobs
            
            # Submit the job in the background
            echo "Submitting job $idx"
            Rscript -e "SBayesRC::LDstep2(outDir='$outDir', blockIndex=$idx, log2file=TRUE)" &> "$outDir/LDstep2_${idx}.log" &
            
            # Add the job's PID to the array
            job_pids+=($!)
            echo "Job $idx submitted with PID ${job_pids[-1]}."
        fi
    done

    # Wait for all jobs to finish
    echo "All jobs submitted for $pheno. Waiting for remaining jobs to finish..."
    while [ ${#job_pids[@]} -gt 0 ]; do
        reap_jobs
        sleep 1
    done

    echo "All LDstep2 jobs completed for $pheno!"
done
