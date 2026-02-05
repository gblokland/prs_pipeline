#!/usr/bin/env python3

import argparse
import pandas as pd

# === APOE region info based on GRCh37/hg19 ===
# rs429358: chr19:45411941
# rs7412:   chr19:45412079
# Region: 1 Mb before rs429358 to 1 Mb after rs7412
# => chr19:44411941 - 46412079

APOE_CHR = '19'
APOE_START = 44411941
APOE_END = 46412079

def remove_apoe_region(input_file, output_file, snp_col='SNP', chr_col='CHR', bp_col='BP', sep='\t'):
    # Load GWAS summary stats
    df = pd.read_csv(input_file, sep=sep)

    # Ensure chromosome column is string for consistent comparison
    df[chr_col] = df[chr_col].astype(str)

    # Filter out SNPs in the APOE region
    mask = ~((df[chr_col] == APOE_CHR) &
             (df[bp_col] >= APOE_START) &
             (df[bp_col] <= APOE_END))

    filtered_df = df[mask]

    # Save filtered file
    filtered_df.to_csv(output_file, sep=sep, index=False)
    print(f"Filtered summary statistics saved to: {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Remove SNPs in APOE region from GWAS summary statistics")
    parser.add_argument("--in", dest="input_file", required=True, help="Path to input GWAS summary stats (TSV/CSV)")
    parser.add_argument("--out", dest="output_file", required=True, help="Path to output file")
    parser.add_argument("--snp_col", default="SNP", help="Column name for SNP ID (default: SNP)")
    parser.add_argument("--chr_col", default="CHR", help="Column name for chromosome (default: CHR)")
    parser.add_argument("--bp_col", default="BP", help="Column name for base pair position (default: BP)")
    parser.add_argument("--sep", default="\t", help="Field separator (default: tab)")

    args = parser.parse_args()

    remove_apoe_region(
        input_file=args.input_file,
        output_file=args.output_file,
        snp_col=args.snp_col,
        chr_col=args.chr_col,
        bp_col=args.bp_col,
        sep=args.sep
    )

