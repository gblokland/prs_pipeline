#!/bin/bash
#Merge vcf files of chr 1-22 into one vcf
concat -o merged_whole_genome.vcf.gz -O z ALL.chr{1..22}.phase3_v5.shapeit2_mvncall_integrated.noSingleton.genotypes.vcf.gz
bcftools index merged_whole_genome.vcf.gz
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AF\n' merged_whole_genome.vcf.gz | awk '{ if ($6 < 0.5) print $0; else print $1, $2, $3, $4, $5, 1-$6 }' > merged_whole_genome_maf_output.txt #The python script to merge with sumstats is supposed to do this but it doesn't work

