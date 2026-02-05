#!/usr/bin/env python3

import argparse
import pandas as pd

# APOE region (hg19): chr19:44411941–46412079
APOE_CHR = '19'
APOE_START = 44411941
APOE_END = 46412079

def load_reference_map(mapping_file, sep2='\t'):
    # Expect columns: SNP (rsID), CHR, BP
    ref_df = pd.read_csv(mapping_file, sep=sep2, dtype=str)
    ref_df['BP'] = ref_df['BP'].astype(int)
    return ref_df

def remove_apoe_by_rsid(sumstats_file, mapping_file, output_file, sep1=' ', sep2='\t'):
    # Load data
    sumstats = pd.read_csv(sumstats_file, sep=sep1, dtype=str)
    ref_map = load_reference_map(mapping_file, sep2=sep2)

    # Merge to get CHR and BP
    merged = pd.merge(sumstats, ref_map, on='SNP', how='left')

    # Filter out SNPs in the APOE region
    mask = ~((merged['CHR'] == APOE_CHR) &
             (merged['BP'].between(APOE_START, APOE_END)))
    
    filtered = merged[mask]

    # Drop CHR and BP columns if they were not in the original data
    filtered = filtered[sumstats.columns]

    # Save
    filtered.to_csv(output_file, sep=sep2, index=False)
    print(f"Filtered summary statistics saved to: {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Remove SNPs near APOE from GWAS by rsID")
    parser.add_argument("--sumstats_file", required=True, help="Input GWAS summary stats file (with rsIDs)")
    parser.add_argument("--mapping_file", required=True, help="Reference file with SNP, CHR, BP columns")
    parser.add_argument("--output_file", required=True, help="Output file")
    parser.add_argument("--sep1", default='\t', help="Field separator (default: tab)")
    parser.add_argument("--sep2", default='\t', help="Field separator (default: tab)")

    args = parser.parse_args()

    remove_apoe_by_rsid(args.sumstats_file, args.mapping_file, args.output_file, args.sep1, args.sep2)
