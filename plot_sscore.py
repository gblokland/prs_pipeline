#!/usr/bin/env python3
# Author: Gabriella Blokland, MHeNs, FHML, UM, 2024-11
# Generates two types of plots for each SCORE column:
# Overall Histogram: Shows the distribution of the scores without splitting by PHENO1.
# Histogram Split by PHENO1: Overlays histograms for each unique PHENO1 value.

import pandas as pd
import matplotlib.pyplot as plt
import argparse
import os

def plot_merged_sscore(input_file, output_dir):
    # Step 1: Read the merged .sscore file
    try:
        df = pd.read_csv(input_file, delim_whitespace=True)
    except Exception as e:
        print(f"Error reading file: {e}")
        return

    # Step 2: Identify SCORE columns to plot
    score_columns = [col for col in df.columns if col.startswith("SCORE")]

    if not score_columns:
        print("No SCORE columns found in the file to plot.")
        return

    # Step 3: Check for PHENO1 column
    if "PHENO1" not in df.columns:
        print("PHENO1 column not found. Only overall plots will be generated.")
        pheno_available = False
    else:
        pheno_available = True

    # Step 4: Plot each SCORE column
    for score_col in score_columns:
        # Overall Histogram
        plt.figure(figsize=(10, 6))
        plt.hist(df[score_col].dropna(), bins=50, alpha=0.7, color="blue", edgecolor="black")
        plt.title(f"Overall Distribution of {score_col}", fontsize=14)
        plt.xlabel("Score", fontsize=12)
        plt.ylabel("Frequency", fontsize=12)
        plt.grid(axis="y", alpha=0.75)

        # Save the overall plot
        overall_output_path = os.path.join(output_dir, f"{score_col}_overall_distribution.png")
        plt.savefig(overall_output_path)
        plt.close()
        print(f"Overall plot saved: {overall_output_path}")

        # Histogram Split by PHENO1
        if pheno_available:
            plt.figure(figsize=(10, 6))
            for pheno in sorted(df["PHENO1"].unique()):
                subset = df[df["PHENO1"] == pheno][score_col].dropna()
                plt.hist(
                    subset, bins=50, alpha=0.5, label=f"PHENO1 = {pheno}", edgecolor="black"
                )

            # Add labels and legend
            plt.title(f"Distribution of {score_col} by PHENO1", fontsize=14)
            plt.xlabel("Score", fontsize=12)
            plt.ylabel("Frequency", fontsize=12)
            plt.legend(title="PHENO1", fontsize=10)
            plt.grid(axis="y", alpha=0.75)

            # Save the split-by-PHENO1 plot
            pheno_output_path = os.path.join(output_dir, f"{score_col}_by_PHENO1_distribution.png")
            plt.savefig(pheno_output_path)
            plt.close()
            print(f"PHENO1-split plot saved: {pheno_output_path}")

def main():
    # Command-line argument parser
    parser = argparse.ArgumentParser(description="Plot distributions from merged .sscore file.")
    parser.add_argument(
        "-i", "--input_file", required=True, help="Path to the merged .sscore file."
    )
    parser.add_argument(
        "-o", "--output_dir", required=True, help="Directory to save the plots."
    )

    # Parse the arguments
    args = parser.parse_args()

    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)

    # Call the plot function
    plot_merged_sscore(args.input_file, args.output_dir)

if __name__ == "__main__":
    main()

