#!/usr/bin/env python3

import pandas as pd
from glob import glob
import argparse
import os
from sklearn.preprocessing import StandardScaler

def merge_sscore_files(input_dir, output_file):
    # Step 1: Get all .sscore files in the input directory
    sscore_files = glob(os.path.join(input_dir, "*.sscore"))
    
    if not sscore_files:
        print(f"No .sscore files found in directory: {input_dir}")
        return
    
    # Step 2: Read and merge all .sscore files
    merged_scores = None

    for file in sscore_files:
        # Read the current .sscore file
        df = pd.read_csv(file, delim_whitespace=True)

        # Ensure the FID column doesn't have a hashtag (#) at the beginning
        if '#FID' in df.columns:
            df['FID'] = df['#FID'].astype(str).str.lstrip('#')  # Convert to string and remove leading hashtag
            df.drop(columns=['#FID'], inplace=True)  # Drop the old column with the hashtag

        # Identify columns to rename
        rename_columns = [
            "SCORE1_AVG",
            "ALLELE_CT",
            "NAMED_ALLELE_DOSAGE_SUM",
        ]
        for col in rename_columns:
            if col in df.columns:
                df.rename(columns={col: f"{col}_{os.path.basename(file)}"}, inplace=True)
        
        # Merge with existing data
        if merged_scores is None:
            merged_scores = df
        else:
            merged_scores = pd.merge(merged_scores, df, on=["FID", "IID", "PHENO1"], how="outer")

    # Step 3: Standardize the PRS columns to mean=0, SD=1
    prs_columns = [col for col in merged_scores.columns if 'SCORE1_AVG' in col]  # Adjust if PRS column names differ
    if prs_columns:
        scaler = StandardScaler()
        # Keep original PRS columns and standardize
        for prs_col in prs_columns:
            merged_scores[f"{prs_col}_std"] = scaler.fit_transform(merged_scores[[prs_col]])
            print(f"Standardized column {prs_col} to mean = 0, SD = 1.")
    else:
        print("No PRS columns found to standardize.")
    
    # Step 4: Drop duplicate PHENO1 columns (if any remain)
    merged_scores = remove_duplicate_columns(merged_scores, "PHENO1")

    # Step 5: Save the merged result to the output file
    merged_scores.to_csv(output_file, sep="\t", index=False)
    print(f"Merged .sscore files saved to: {output_file}")

def remove_duplicate_columns(df, column_name):
    """
    Ensures there is only one column with the specified name, dropping others if they exist.
    """
    # Identify columns matching the specified name
    matching_columns = [col for col in df.columns if col == column_name]
    
    # If duplicates exist, drop all but the first instance
    if len(matching_columns) > 1:
        df = df.loc[:, ~df.columns.duplicated()]  # Drop duplicate columns with identical names
    
    return df

def main():
    # Command-line argument parser
    parser = argparse.ArgumentParser(description="Merge PLINK2 .sscore files into a single file and standardize PRS.")
    parser.add_argument(
        "-i", "--input_dir", required=True, help="Directory containing .sscore files to merge."
    )
    parser.add_argument(
        "-o", "--output_file", required=True, help="Name of the output merged file."
    )
    
    # Parse the arguments
    args = parser.parse_args()
    
    # Call the merge function
    merge_sscore_files(args.input_dir, args.output_file)

if __name__ == "__main__":
    main()
