#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
from scipy.stats import pearsonr
from sklearn.linear_model import LinearRegression
import argparse

def plot_prs_vs_pc1(pca_file, prs_file, prs_column, output_dir):
    try:
        # Load the PCA and PRS data
        df_pca = pd.read_csv(pca_file, sep="\t")
        df_prs = pd.read_csv(prs_file, sep="\t")
    except Exception as e:
        print(f"Error reading files: {e}")
        return

    print("PCA file columns:", df_pca.columns)
    print("PRS file columns:", df_prs.columns)
    
    # Ensure both dataframes contain the necessary columns
    if 'PC1' not in df_pca.columns:
        print("Error: 'PC1' column not found in PCA file.")
        return
    if 'IID' not in df_pca.columns or 'IID' not in df_prs.columns:
        print("Error: 'IID' column not found in one or both files.")
        return
    if prs_column not in df_prs.columns:
        print(f"Error: '{prs_column}' column not found in PRS file.")
        return

    # Merge the PCA and PRS data on 'IID'
    df_merged = pd.merge(df_pca, df_prs, on='IID', how='inner')

    # Scatter plot of PRS vs PC1
    plt.figure(figsize=(10, 6))
    plt.scatter(df_merged['PC1'], df_merged[prs_column], alpha=0.5, label="Data points")
    
    # Add a regression line
    X = df_merged['PC1'].values.reshape(-1, 1)
    y = df_merged[prs_column].values
    model = LinearRegression()
    model.fit(X, y)
    y_pred = model.predict(X)
    plt.plot(df_merged['PC1'], y_pred, color='red', label="Regression line")
    
    # Calculate R² and p-value for the correlation
    r_value, p_value = pearsonr(df_merged['PC1'], df_merged[prs_column])
    r_squared = r_value ** 2

    # Annotate statistics on the plot
    stats_text = (
        f"R² = {r_squared:.4f}\n"
        f"p-value = {p_value:.4e}"
    )
    plt.text(0.05, 0.95, stats_text, transform=plt.gca().transAxes, 
             fontsize=10, color="blue", verticalalignment='top')

    plt.title("PRS vs PC1 (Population Stratification Check)")
    plt.xlabel("PC1")
    plt.ylabel(f"Polygenic Risk Score ({prs_column})")
    plt.legend()

    # Extract the base name of the PRS file without the directory and extension
    prs_base_name = os.path.splitext(os.path.basename(prs_file))[0]

    # Save the plot to the output directory with a dynamic filename based on PRS file name
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, f"{prs_base_name}_PRS_vs_PC1.png")
    plt.savefig(output_file)
    print(f"Plot saved to {output_file}")

def main():
    parser = argparse.ArgumentParser(description="Plot PRS vs PC1 after merging files on IID and save the plot.")
    
    # Define the arguments
    parser.add_argument(
        "-pca", "--pca_file", required=True, help="Path to the PCA file (must include 'PC1' and 'IID' columns)."
    )
    parser.add_argument(
        "-prs", "--prs_file", required=True, help="Path to the PRS file (must include 'IID' and specified PRS column)."
    )
    parser.add_argument(
        "-prs_col", "--prs_column", required=True, help="Name of the PRS column in the PRS file."
    )
    parser.add_argument(
        "-o", "--output_dir", required=True, help="Directory where the plot will be saved."
    )

    # Parse the arguments
    args = parser.parse_args()

    # Call the plot function with parsed arguments
    plot_prs_vs_pc1(args.pca_file, args.prs_file, args.prs_column, args.output_dir)

if __name__ == "__main__":
    main()
