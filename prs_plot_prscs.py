#!/usr/bin/env python3
# Author: Gabriella Blokland, MHeNs, FHML, UM, 2024-11
# Scatter Plot of Betas Across Chromosomes
# with QQ Plot and Dynamic Output Naming

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import argparse
import os
import scipy.stats as stats
import datashader as ds
from datashader import transfer_functions as tf

def plot_prscs_betas(prscs_file, output_dir, has_headers):
    # Step 1: Define column names for PRScs output
    column_names = ["CHR", "BP", "SNP", "A1", "A2", "BETA"]

    # Step 2: Load PRScs output with or without headers
    try:
        if has_headers:
            df = pd.read_csv(prscs_file, delim_whitespace=True)
        else:
            df = pd.read_csv(prscs_file, delim_whitespace=True, header=None, names=column_names)
    except FileNotFoundError:
        print(f"Error: The file {prscs_file} was not found.")
        return
    except pd.errors.ParserError:
        print(f"Error: There was an issue parsing the file {prscs_file}.")
        return
    except Exception as e:
        print(f"Error reading file: {e}")
        return

    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Extract the base name of the input file for dynamic naming
    base_name = os.path.splitext(os.path.basename(prscs_file))[0]

    # Step 3: Plot histogram of betas
    plt.figure(figsize=(10, 6))
    sns.histplot(df["BETA"], bins=100, kde=True, color='blue', edgecolor='black')
    plt.title("Distribution of Betas (Effect Sizes)", fontsize=14)
    plt.xlabel("Beta", fontsize=12)
    plt.ylabel("Frequency", fontsize=12)
    plt.grid(axis="y", alpha=0.75)
    plt.tight_layout()
    hist_output_path = os.path.join(output_dir, f"{base_name}_beta_distribution.png")
    plt.savefig(hist_output_path)
    plt.close()
    print(f"Beta histogram saved: {hist_output_path}")

    # Step 4: QQ plot of the betas
    plt.figure(figsize=(10, 6))
    stats.probplot(df["BETA"].dropna(), dist="norm", plot=plt)
    plt.title("QQ Plot of Betas", fontsize=14)
    plt.grid(True)
    qq_output_path = os.path.join(output_dir, f"{base_name}_beta_qqplot.png")
    plt.savefig(qq_output_path)
    plt.close()
    print(f"QQ plot saved: {qq_output_path}")

    # Step 5: Scatter plot of betas across chromosomes using Datashader
    # Ensure numeric columns
    df['CHR'] = pd.to_numeric(df['CHR'], errors='coerce')
    df['BP'] = pd.to_numeric(df['BP'], errors='coerce')
    df['BETA'] = pd.to_numeric(df['BETA'], errors='coerce')

    # Drop rows with NaN values in key columns
    df = df.dropna(subset=['CHR', 'BP', 'BETA'])
    
    df['BETA'] = df['BETA'].astype(float)

    # Sort data by chromosome and position
    df = df.sort_values(by=['CHR', 'BP']).reset_index(drop=True)

    # Calculate cumulative offsets for genome positions
    chr_offsets = df.groupby('CHR')['BP'].max().cumsum().shift(fill_value=0)
    df['genome_pos'] = df['BP'] + df['CHR'].map(chr_offsets)

    # Debug: Print out some key data points
    print(f"Data points summary: {df[['genome_pos', 'BETA']].head()}")

    # Create a Datashader Canvas with appropriate dimensions
    cvs = ds.Canvas(plot_width=800, plot_height=100)

    # Aggregate data points using the canvas (without using 'Points')
    agg = cvs.points(df, 'genome_pos', 'BETA')  # Aggregating the data points based on genome_pos and BETA

    # Debug: Check the aggregated data to ensure it's being properly processed
    print(f"Aggregated data summary: {agg}")

    # Apply a transfer function to map the aggregate to an image
    img = tf.shade(agg, cmap='Viridis', how='linear')

    # Plot with matplotlib
    fig, ax = plt.subplots(figsize=(16, 2))
    ax.imshow(img.to_pil(), origin='lower')
    ax.set_title("Scatter Plot of Betas Across Chromosomes")
    ax.set_xlabel("Genome Position")
    ax.set_ylabel("Effect Size (Beta)")
    
    #ax.set_yscale('log')
    #ax.set_ylim(-0.03, 0.03)  # Modify the limits according to your data range

    # Save the plot
    scatter_output_path = os.path.join(output_dir, f"{base_name}_beta_scatter.png")
    plt.savefig(scatter_output_path)
    plt.close()
    print(f"Beta scatter plot saved: {scatter_output_path}")
    
def main():
    # Command-line argument parser
    parser = argparse.ArgumentParser(description="Plot PRScs Betas.")
    parser.add_argument(
        "-i", "--input_file", required=True, help="Path to the PRScs output file.")
    parser.add_argument(
        "-o", "--output_dir", required=True, help="Directory to save the plots.")
    parser.add_argument(
        "--has_headers", action="store_true", help="Specify if the input file has headers.")

    # Parse arguments
    args = parser.parse_args()

    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)

    # Call the plot function
    plot_prscs_betas(args.input_file, args.output_dir, args.has_headers)

if __name__ == "__main__":
    main()
