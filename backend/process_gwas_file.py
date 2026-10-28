import pysam
import csv
import argparse
import os

# Cache opened tabix files per chromosome
tabix_cache = {}

def load_chrom_vcf(vcf_dir, chrom):
    chrom = chrom.replace("chr", "")  # Normalize: remove 'chr' prefix if present
    filename = os.path.join(vcf_dir, f"homo_sapiens_chr{chrom}.vcf.gz")
    
    if filename not in tabix_cache:
        if not os.path.exists(filename):
            raise FileNotFoundError(f"Missing VCF file for chromosome {chrom}: {filename}")
        tabix_cache[filename] = pysam.TabixFile(filename)
    return tabix_cache[filename]

def get_rsid(vcf_dir, chrom, pos, ref=None, alt=None):
    try:
        tabix_file = load_chrom_vcf(vcf_dir, chrom)
        if chrom.startswith("chr"):
            chrom = chrom[3:]

        records = tabix_file.fetch(chrom, int(pos) - 1, int(pos))
        for record in records:
            fields = record.strip().split('\t')
            vcf_pos = fields[1]
            rsid = fields[2]
            vcf_ref = fields[3]
            vcf_alt = fields[4]

            if vcf_pos == pos:
                if ref and alt:
                    if vcf_ref == ref and alt in vcf_alt.split(','):
                        return rsid
                else:
                    return rsid
    except Exception:
        return "NA"
    return "NA"

def process_gwas_file(input_file, output_file, vcf_dir):
    with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.reader(infile, delimiter='\t')
        writer = csv.writer(outfile, delimiter='\t')

        header = next(reader)
        writer.writerow(['rsID'] + header[1:])  # Drop original SNP ID

        for row in reader:
            snp = row[0]
            other_data = row[1:]

            try:
                chrom, pos, allele_info = snp.split(":")
                ref, alt = allele_info.split("_")
                rsid = get_rsid(vcf_dir, chrom, pos, ref, alt)
                writer.writerow([rsid] + other_data)
            except Exception as e:
                writer.writerow([f"Error: {e}"] + other_data)

def main():
    parser = argparse.ArgumentParser(description="Convert CHR:POS:REF_ALT to rsID using chromosome-split dbSNP VCFs.")
    parser.add_argument("input_file", help="Input GWAS file with SNPs in CHR:POS:REF_ALT format")
    parser.add_argument("output_file", help="Output file with rsIDs")
    parser.add_argument("vcf_dir", help="Directory containing chromosome-specific VCF files (e.g., homo_sapiens_chr1.vcf.gz)")

    args = parser.parse_args()
    process_gwas_file(args.input_file, args.output_file, args.vcf_dir)

if __name__ == "__main__":
    main()
