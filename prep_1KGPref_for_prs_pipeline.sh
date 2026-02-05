#!/bin/bash

##################################################################################################
mkdir -p ~/ref_panels/1KGPref
cd ~/ref_panels/1KGPref
#wget -r -np -R "index.html*, supporting" -e robots=off --cut-dirs=5 ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/
#cd ftp.1000genomes.ebi.ac.uk/

#Download ALL.chr${i}.phase3_v5.shapeit2_mvncall_integrated.noSingleton.genotypes.vcf.gz 
if [ ! -f ~/ref_panels/1KGPref/integrated_call_samples_v3.20130502.ALL.panel ]; then
    mkdir -p ~/ref_panels/1KGPref
    cd ~/ref_panels/1KGPref/
    for i in {1..22} X Y MT; do
        wget -c "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
        wget -c "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi"
    done
    wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel
    wget ftp://yunlianon:anon@rc-ns-ftp.its.unc.edu/ALL.phase3_v5.shapeit2_mvncall_integrated.noSingleton.tgz
    tar xzf ALL.phase3_v5.shapeit2_mvncall_integrated.noSingleton.tgz
    mv all/* ./
    rmdir all
fi

#Adding allele frequencies if they are missing from the sumstats
grep 'CEU' ~/ref_panels/1KGPref/integrated_call_samples_v3.20130502.ALL.panel | awk '{print $1}' > ~/ref_panels/1KGPref/CEU_samples.txt
grep 'EUR' ~/ref_panels/1KGPref/integrated_call_samples_v3.20130502.ALL.panel | awk '{print $1}' > ~/ref_panels/1KGPref/EUR_samples.txt
grep 'EAS' ~/ref_panels/1KGPref/integrated_call_samples_v3.20130502.ALL.panel | awk '{print $1}' > ~/ref_panels/1KGPref/EAS_samples.txt
grep 'SAS' ~/ref_panels/1KGPref/integrated_call_samples_v3.20130502.ALL.panel | awk '{print $1}' > ~/ref_panels/1KGPref/SAS_samples.txt
grep 'AFR' ~/ref_panels/1KGPref/integrated_call_samples_v3.20130502.ALL.panel | awk '{print $1}' > ~/ref_panels/1KGPref/AFR_samples.txt
grep 'AMR' ~/ref_panels/1KGPref/integrated_call_samples_v3.20130502.ALL.panel | awk '{print $1}' > ~/ref_panels/1KGPref/AMR_samples.txt

#Merge ALL.chr${i}.phase3_v5.shapeit2_mvncall_integrated.noSingleton.genotypes.vcf.gz into merged_whole_genome.vcf.gz
cd ~/ref_panels/1KGPref/
~/code/prs_pipeline/vcf_concat_chr.sh

#Generate EUR_frequencies.frq
~/opt/vcftools --gzvcf ~/ref_panels/1KGPref/merged_whole_genome.vcf.gz --freq \
--keep ~/ref_panels/1KGPref/EUR_samples.txt --out ~/ref_panels/1KGPref/EUR_frequencies
~/opt/vcftools --gzvcf ~/ref_panels/1KGPref/merged_whole_genome.vcf.gz --freq \
--keep ~/ref_panels/1KGPref/EAS_samples.txt --out ~/ref_panels/1KGPref/EAS_frequencies
~/opt/vcftools --gzvcf ~/ref_panels/1KGPref/merged_whole_genome.vcf.gz --freq \
--keep ~/ref_panels/1KGPref/SAS_samples.txt --out ~/ref_panels/1KGPref/SAS_frequencies
~/opt/vcftools --gzvcf ~/ref_panels/1KGPref/merged_whole_genome.vcf.gz --freq \
--keep ~/ref_panels/1KGPref/AFR_samples.txt --out ~/ref_panels/1KGPref/AFR_frequencies
~/opt/vcftools --gzvcf ~/ref_panels/1KGPref/merged_whole_genome.vcf.gz --freq \
--keep ~/ref_panels/1KGPref/AMR_samples.txt --out ~/ref_panels/1KGPref/AMR_frequencies

cat ~/ref_panels/1KGPref/EUR_frequencies.frq | sed 's/{ALLELE:FREQ}/{ALLELE:FREQ}.1\t{ALLELE:FREQ}.2/g' > ~/ref_panels/1KGPref/EUR_frequencies.frq.edit
python ~/opt/prs_pipeline/split_allele_freq.py ~/ref_panels/1KGPref/EUR_frequencies.frq.edit ~/ref_panels/1KGPref/EUR_frequencies.frq.tsv

#Extract SNP rs ids and CHROM and POS from vcf file and write to tsv file
python ~/opt/prs_pipeline/extract_vcf_info.py ~/ref_panels/1KGPref/merged_whole_genome.vcf.gz ~/ref_panels/1KGPref/merged_whole_genome_snps.tsv

##Merge SNP info with frq file
#python ~/opt/prs_pipeline/merge_and_sort_files.py ~/ref_panels/1KGPref/merged_whole_genome_snps.tsv \
#~/ref_panels/1KGPref/EUR_frequencies.frq.tsv \
#~/ref_panels/1KGPref/EUR_frequencies.tsv

cd ~/ref_panels/
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/snp147Common.txt.gz
gunzip snp147Common.txt.gz
###cut -f2-5 snp147Common.txt | awk '{gsub("chr", "", $2); print $4"\t$1"\t"$2"\t"($3+1)}' > rsid_chr_bp_map.tsv
cut -f2-5 snp147Common.txt | awk 'BEGIN {print "SNP\tCHR\tBP\tBP_END"} {gsub("chr", "", $1); print $4"\t"$1"\t"$2"\t"$3}' > rsid_chr_bp_map.tsv

##################################################################################################
#Run format conversion
cd ~/ref_panels
mkdir -p ensembl_grch37 && cd ensembl_grch37

for chr in {1..22} X Y MT; do
wget ftp://ftp.ensembl.org/pub/grch37/release-105/variation/vcf/homo_sapiens/homo_sapiens-chr${chr}.vcf.gz
wget ftp://ftp.ensembl.org/pub/grch37/release-105/variation/vcf/homo_sapiens/homo_sapiens-chr${chr}.vcf.gz.csi
tabix -p vcf homo_sapiens-chr${chr}.vcf.gz #make *.tbi file
done

wget -r -nH --cut-dirs=5 -A "homo_sapiens-chr*.vcf.gz,homo_sapiens-chr*.vcf.gz.tbi" ftp://ftp.ensembl.org/pub/grch37/current/variation/vcf/homo_sapiens/

cd homo_sapiens

bcftools concat -Oz -o all_chromosomes_merged.vcf.gz homo_sapiens-chr{1..22}.vcf.gz homo_sapiens-chrX.vcf.gz homo_sapiens-chrY.vcf.gz
tabix -p vcf all_chromosomes_merged.vcf.gz

wget https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.25.gz -O all_chromosomes.vcf.gz
tabix -p vcf all_chromosomes.vcf.gz
