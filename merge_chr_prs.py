#!/usr/bin/env python3

import pandas as pd
import argparse
import os

def main(pattern, chrom_range, output_file):
    prs_total = None
    chromosomes = parse_chrom_range(chrom_range)

    for chr in chromosomes:
        filename = pattern.format(chr=chr)
        if not os.path.exists(filename):
            print(f"Warning: File '{filename}' not found, skipping.")
            continue

        df = pd.read_csv(filename, delim_whitespace=True)
        df.columns = [col.lstrip('#') for col in df.columns]
        if 'SCORE1_AVG' not in df.columns:
            print(f"Error: 'SCORE1_AVG' column not found in {filename}. Skipping.")
            continue

        df = df[['FID', 'IID', 'SCORE1_AVG']]
        df.rename(columns={'SCORE1_AVG': f'PRS_chr{chr}'}, inplace=True)

        if prs_total is None:
            prs_total = df
        else:
            prs_total = prs_total.merge(df, on=['FID', 'IID'])

    if prs_total is None:
        print("No valid PRS files found. Exiting.")
        return

    # Sum across all PRS_chr columns
    score_columns = [col for col in prs_total.columns if col.startswith('PRS_chr')]
    prs_total['PRS_total'] = prs_total[score_columns].sum(axis=1)

    prs_total.to_csv(output_file, index=False)
    print(f"Saved merged PRS to '{output_file}'.")


def parse_chrom_range(chrom_range_str):
    try:
        if '-' in chrom_range_str:
            start, end = map(int, chrom_range_str.split('-'))
            return range(start, end + 1)
        else:
            return [int(c) for c in chrom_range_str.split(',')]
    except ValueError:
        raise ValueError("Chromosome range must be like '1-22' or '1,2,3'.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Merge per-chromosome PRS .sscore files into a total PRS score.')
    parser.add_argument('--pattern', type=str, default='chr{chr}_prs.sscore',
                        help="Input file pattern with '{chr}' as placeholder (default: chr{chr}_prs.sscore)")
    parser.add_argument('--chrom', type=str, default='1-22',
                        help="Chromosome range (e.g., '1-22' or '1,2,3') (default: 1-22)")
    parser.add_argument('--out', type=str, default='total_PRS.csv',
                        help='Output file name (default: total_PRS.csv)')
    args = parser.parse_args()

    main(args.pattern, args.chrom, args.out)
