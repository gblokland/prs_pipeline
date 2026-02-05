#!/bin/bash

#Harmonizing formatting of GWAS summary statistics

#install a set op commonly used python libraries from the bash command line
cd ~/opt/
echo "
numpy
pandas
matplotlib
scikit-learn
seaborn
requests
flask
django
jupyterlab
" > python_requirements_shortlist.txt

pip install -r python_requirements_shortlist.txt

#Resources for semi-automatic QC and GWAS sumstats formatting:
#https://github.com/BioPsyk/cleansumstats 
#https://academic.oup.com/bioinformatics/article/28/3/444/189687
#https://github.com/precimed/python_convert
cd ~/opt/
if [ ! -d python_convert ]; then
	git clone https://github.com/precimed/python_convert.git
fi
cd python_convert/
chmod +x *.py
if [ ! -f 2558411_ref.bim ]; then
	wget https://precimed.s3-eu-west-1.amazonaws.com/python_convert/2558411_ref.bim
	wget https://precimed.s3-eu-west-1.amazonaws.com/python_convert/9279485_ref.bim
	wget https://precimed.s3-eu-west-1.amazonaws.com/python_convert/b149_RsMergeArch.bcp.gz
	wget https://precimed.s3-eu-west-1.amazonaws.com/python_convert/b149_SNPChrPosOnRef_105.bcp.gz
	wget https://precimed.s3-eu-west-1.amazonaws.com/python_convert/b149_SNPHistory.bcp.gz
	wget https://precimed.s3-eu-west-1.amazonaws.com/python_convert/hg18ToHg19.over.chain.gz
	wget https://precimed.s3-eu-west-1.amazonaws.com/python_convert/ref_1kG_phase3_EUR.tar.gz
	tar xzf ref_1kG_phase3_EUR.tar.gz
fi

#CRP - https://doi.org/10.1038/s41467-022-29650-5 
mkdir -p /root/persistent/sumstats/Biomarkers/CRP && cd /root/persistent/sumstats/Biomarkers/CRP/
wget -r -np -R "index.html*" -e robots=off --cut-dirs=8 https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90029001-GCST90030000/GCST90029070/harmonised/

#Lymphocyte count GWAS: PMID 27863252 x adaptive immune response - https://doi.org/10.1016/j.cell.2016.10.042
mkdir -p /root/persistent/sumstats/Biomarkers/lymphocyte && cd /root/persistent/sumstats/Biomarkers/lymphocyte/
wget -r -np -R "index.html*" -e robots=off --cut-dirs=8 https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST004001-GCST005000/GCST004632/harmonised/

#Granulocyte count GWAS: PMID 27863252 x innate immune response - https://doi.org/10.1016/j.cell.2016.10.042 
mkdir -p /root/persistent/sumstats/Biomarkers/granulocyte && cd /root/persistent/sumstats/Biomarkers/granulocyte/
wget -r -np -R "index.html*" -e robots=off --cut-dirs=8 https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST004001-GCST005000/GCST004614/harmonised/

#Granulocytes are a type of white blood cell that has small granules inside them. These granules contain proteins. 
#The specific types of granulocytes are neutrophils, eosinophils, and basophils.



cd /root/persistent/sumstats/Biomarkers/CRP/
python ~/opt/python_convert/sumstats.py csv \
--sumstats 35459240-GCST90029070-EFO_0004458-Build37.f.tsv.gz --out CRP_sumstats.tsv \
--force --auto --head 5 --chr chromosome --bp base_pair_location --snp variant_id --beta beta --se standard_error --pval p_value --n-val 575531
#UKB n=427367; CHARGE n = 148164; total n = 575531
###python $(python_convert)/sumstats.py mat --sumstats PGC_SCZ_2014.csv --out PGC_SCZ_2014.mat --ref 2558411_ref.bim --force


cd /root/persistent/sumstats/Biomarkers/lymphocyte/
python ~/opt/python_convert/sumstats.py csv \
--sumstats 27863252-GCST004632-EFO_0007993-Build37.f.tsv.gz --out lymphocyte_sumstats.tsv \
--force --auto --head 5 --chr chromosome --bp base_pair_location --snp variant_id --beta beta --se standard_error --pval p_value --n-val 173480


cd /root/persistent/sumstats/Biomarkers/granulocyte/
python ~/opt/python_convert/sumstats.py csv \
--sumstats 27863252-GCST004614-EFO_0007987-Build37.f.tsv.gz --out granulocyte_sumstats.tsv \
--force --auto --head 5 --chr chromosome --bp base_pair_location --snp variant_id --beta beta --se standard_error --pval p_value --n-val 173480


#For PRScs: SNP	A1	A2	BETA	SE
for phenotype in CRP lymphocyte granulocyte; do
	awk '{print $1, $5, $6, $8, $9}' /root/persistent/sumstats/Biomarkers/${phenotype}/${phenotype}_sumstats.tsv > /root/persistent/sumstats/Biomarkers/${phenotype}/${phenotype}_sumstats_PRScs.tsv
done


#For SBayesRC: 'SNP', 'A1', 'A2', 'Freq', 'b', 'se', 'p', 'N'

#Adding allele frequencies if they are missing from the sumstats
grep 'CEU' ~/persistent/ref_panels/1KGPref/integrated_call_samples_v3.20130502.ALL.panel | awk '{print $1}' > ~/persistent/ref_panels/1KGPref/CEU_samples.txt
grep 'EUR' ~/persistent/ref_panels/1KGPref/integrated_call_samples_v3.20130502.ALL.panel | awk '{print $1}' > ~/persistent/ref_panels/1KGPref/EUR_samples.txt
grep 'EAS' ~/persistent/ref_panels/1KGPref/integrated_call_samples_v3.20130502.ALL.panel | awk '{print $1}' > ~/persistent/ref_panels/1KGPref/EAS_samples.txt
grep 'SAS' ~/persistent/ref_panels/1KGPref/integrated_call_samples_v3.20130502.ALL.panel | awk '{print $1}' > ~/persistent/ref_panels/1KGPref/SAS_samples.txt
grep 'AFR' ~/persistent/ref_panels/1KGPref/integrated_call_samples_v3.20130502.ALL.panel | awk '{print $1}' > ~/persistent/ref_panels/1KGPref/AFR_samples.txt
grep 'AMR' ~/persistent/ref_panels/1KGPref/integrated_call_samples_v3.20130502.ALL.panel | awk '{print $1}' > ~/persistent/ref_panels/1KGPref/AMR_samples.txt

vcftools --gzvcf ~/persistent/ref_panels/1KGPref/merged_whole_genome.vcf.gz --freq \
--keep ~/persistent/ref_panels/1KGPref/EUR_samples.txt --out ~/persistent/ref_panels/1KGPref/EUR_frequencies
vcftools --gzvcf ~/persistent/ref_panels/1KGPref/merged_whole_genome.vcf.gz --freq \
--keep ~/persistent/ref_panels/1KGPref/EAS_samples.txt --out ~/persistent/ref_panels/1KGPref/EAS_frequencies
vcftools --gzvcf ~/persistent/ref_panels/1KGPref/merged_whole_genome.vcf.gz --freq \
--keep ~/persistent/ref_panels/1KGPref/SAS_samples.txt --out ~/persistent/ref_panels/1KGPref/SAS_frequencies
vcftools --gzvcf ~/persistent/ref_panels/1KGPref/merged_whole_genome.vcf.gz --freq \
--keep ~/persistent/ref_panels/1KGPref/AFR_samples.txt --out ~/persistent/ref_panels/1KGPref/AFR_frequencies
vcftools --gzvcf ~/persistent/ref_panels/1KGPref/merged_whole_genome.vcf.gz --freq \
--keep ~/persistent/ref_panels/1KGPref/AMR_samples.txt --out ~/persistent/ref_panels/1KGPref/AMR_frequencies

cat ~/persistent/ref_panels/1KGPref/EUR_frequencies.frq | sed 's/{ALLELE:FREQ}/{ALLELE:FREQ}.1\t{ALLELE:FREQ}.2/g' > ~/persistent/ref_panels/1KGPref/EUR_frequencies.frq.edit
python ~/ozpt/split_allele_freq.py ~/persistent/ref_panels/1KGPref/EUR_frequencies.frq.edit ~/persistent/ref_panels/1KGPref/EUR_frequencies.frq.tsv

#Extract SNP rs ids and CHROM and POS from vcf file and write to tsv file
python ~/opt/extract_vcf_info.py ~/persistent/ref_panels/1KGPref/merged_whole_genome.vcf.gz ~/persistent/ref_panels/1KGPref/merged_whole_genome_snps.tsv

##Merge SNP info with frq file
#python ~/opt/merge_and_sort_files.py ~/persistent/ref_panels/1KGPref/merged_whole_genome_snps.tsv \
#~/persistent/ref_panels/1KGPref/EUR_frequencies.frq.tsv \
#~/persistent/ref_panels/1KGPref/EUR_frequencies.tsv


for phenotype in CRP lymphocyte granulocyte; do
	python ~/opt/convert_sumstats_to_cojo.py /root/persistent/sumstats/Biomarkers/${phenotype}/${phenotype}_sumstats.tsv \
	~/persistent/ref_panels/1KGPref/EUR_frequencies.frq.tsv \
	~/persistent/sumstats/Biomarkers/${phenotype}/${phenotype}_sumstats_COJO.tsv
done

for phenotype in CRP lymphocyte granulocyte; do
	python ~/opt/OR_beta_conversion.py ~/persistent/sumstats/Biomarkers/${phenotype}/${phenotype}_sumstats_OR.tsv \
	~/persistent/sumstats/Biomarkers/${phenotype}/${phenotype}_sumstats_BETA.tsv

done




cd ~/persistent/ref_panels/1KGPref

wget -r -np -R "index.html*, supporting" -e robots=off --cut-dirs=5 ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/

cd ftp.1000genomes.ebi.ac.uk/

#  add_maf_to_sumstats(gwas_file, maf_file, output_file) #The following scripts don't work yet

gwas_file="/root/persistent/sumstats/Biomarkers/CRP/CRP_sumstats.tsv"  # Path to your GWAS summary statistics file
maf_file='path_to_maf_file.txt'
output_file="/root/persistent/sumstats/Biomarkers/CRP/CRP_sumstats_wMAF.tsv"  # Path to the output file

python ~/opt/add_maf_to_gwas_ensembl.py gwas_file output_file --snp_col SNP --population "1000GENOMES:phase_3:EUR"

python ~/opt/add_maf_to_gwas_ensembl.py /root/persistent/sumstats/Biomarkers/CRP/CRP_sumstats.tsv \
/root/persistent/sumstats/Biomarkers/CRP/CRP_sumstats_wMAF.tsv --snp_col SNP --population "1000GENOMES:phase_3:EUR"

python ~/opt/add_maf_to_gwas.py /root/persistent/sumstats/Biomarkers/CRP/CRP_sumstats.tsv \
/root/persistent/sumstats/Biomarkers/CRP/CRP_sumstats_wMAF.tsv \
--source vcf --vcf_file ~/persistent/ref_panels/1KGPref/merged_whole_genome.vcf.gz --snp_col SNP


python ~/opt/add_maf_to_sumstats_final.py /root/persistent/sumstats/Biomarkers/CRP/CRP_sumstats.tsv \
/root/persistent/sumstats/Biomarkers/CRP/CRP_sumstats_wMAF.tsv \
--source vcf --vcf_file ~/persistent/ref_panels/1KGPref/merged_whole_genome.vcf.gz --snp_col SNP
