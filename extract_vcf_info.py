import pysam
import argparse
import re

def extract_vcf_info(vcf_file):
    # Open the VCF file (gzipped or not) using pysam
    vcf = pysam.VariantFile(vcf_file)
    
    # List to store tuples of (SNP, CHROM, POS)
    vcf_info = []

    # Iterate through each record in the VCF file
    for record in vcf.fetch():
        # Extract SNP (rsID), CHROM, and POS
        snp = record.id
        chrom = record.chrom
        pos = record.pos
        
        # Check if the SNP (rsID) is not missing ('.' represents missing in VCF)
        if snp != '.':
            # Use regex to extract the part of rsID before the first colon, if it exists
            match = re.match(r'(rs\d+)', snp)
            if match:
                snp = match.group(1)  # Extract the actual rsID (e.g., rs12345)
            vcf_info.append((snp, chrom, pos))
    
    return vcf_info

def save_to_file(vcf_info, output_file):
    # Open the output file for writing
    with open(output_file, 'w') as f:
        # Write header
        f.write("SNP\tCHROM\tPOS\n")
        
        # Write each SNP, CHROM, and POS in tab-separated format
        for snp, chrom, pos in vcf_info:
            f.write(f"{snp}\t{chrom}\t{pos}\n")

def main():
    # Set up command line argument parsing
    parser = argparse.ArgumentParser(description="Extract SNPs (rsIDs), chromosomes, and positions from a VCF file, splitting SNP after the actual rsID.")
    parser.add_argument("vcf_file", help="Path to the input VCF or VCF.gz file")
    parser.add_argument("output_file", help="Path to the output text file")
    
    # Parse the command line arguments
    args = parser.parse_args()
    
    # Extract the VCF information
    vcf_info = extract_vcf_info(args.vcf_file)

    # Save the output to a file
    save_to_file(vcf_info, args.output_file)

if __name__ == "__main__":
    main()
