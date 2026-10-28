#!/usr/bin/env python3

import pandas as pd
import argparse
import matplotlib.pyplot as plt
import os
import numpy as np
from scipy.stats import pearsonr, pointbiserialr
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import log_loss

def calculate_correlation_and_plot(prs_file, prs_column, phenotype_column, output_dir):
    # Load the PRS data
    try:
        df_prs = pd.read_csv(prs_file)
    except Exception as e:
        print(f"Error reading file: {e}")
        return

    # Check if the necessary columns exist in the dataframe
    if prs_column not in df_prs.columns or phenotype_column not in df_prs.columns:
        print(f"Error: '{prs_column}' and/or '{phenotype_column}' column not found.")
        return

    # Handle missing values in the phenotype column (-9 is treated as missing)
    df_prs[phenotype_column].replace(-9, pd.NA, inplace=True)
    df_prs.dropna(subset=[phenotype_column, prs_column], inplace=True)

    # Determine if the phenotype is binary or continuous
    phenotype_values = df_prs[phenotype_column].unique()
    if len(phenotype_values) == 2:
        is_binary = True
        print("Detected binary phenotype.")
    else:
        is_binary = False
        print("Detected continuous phenotype.")

    # Perform correlation or logistic regression
    if is_binary:
        # Convert binary phenotype to 0 and 1 if not already
        df_prs[phenotype_column] = pd.Categorical(df_prs[phenotype_column]).codes

        # Logistic regression for pseudo-R-squared
        model = LogisticRegression()
        model.fit(df_prs[[prs_column]], df_prs[phenotype_column])

        # Calculate pseudo-R-squared (McFadden's)
        log_likelihood_full = -log_loss(df_prs[phenotype_column], model.predict_proba(df_prs[[prs_column]]), normalize=False)
        log_likelihood_null = -log_loss(df_prs[phenotype_column], [df_prs[phenotype_column].mean()] * len(df_prs), normalize=False)
        pseudo_r_squared = 1 - (log_likelihood_full / log_likelihood_null)

        # Calculate point-biserial correlation
        corr, p_value = pointbiserialr(df_prs[prs_column], df_prs[phenotype_column])
        print(f"Point-biserial Correlation: {corr:.4f}, P-value: {p_value:.4e}")
        print(f"Pseudo-R-squared (McFadden's): {pseudo_r_squared:.4f}")
        correlation_metric = f"Point-biserial Correlation: {corr:.4f}\nPseudo-R-squared: {pseudo_r_squared:.4f}\nP-value: {p_value:.4e}"

    else:
        # Calculate Pearson correlation
        corr, p_value = pearsonr(df_prs[prs_column], df_prs[phenotype_column])
        r_squared = corr**2
        print(f"Correlation: {corr:.4f}, R-squared: {r_squared:.4f}, P-value: {p_value:.4e}")
        correlation_metric = f"R-squared: {r_squared:.4f}\nP-value: {p_value:.4e}"

    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Save correlation or pseudo R-squared result to a file
    correlation_file = os.path.join(output_dir, f"{prs_column}_vs_{phenotype_column}_correlation.txt")
    with open(correlation_file, 'w') as f:
        if is_binary:
            f.write(f"Correlation between {prs_column} and {phenotype_column}: Point-biserial Correlation: {corr:.4f}\n")
            f.write(f"Correlation between {prs_column} and {phenotype_column}: Pseudo-R-squared (McFadden's): {pseudo_r_squared:.4f}\n")
            f.write(f"Correlation between {prs_column} and {phenotype_column}: P-value: {p_value:.4e}\n")
        else:
            f.write(f"Correlation between {prs_column} and {phenotype_column}: Pearson Correlation: {corr:.4f}\n")
            f.write(f"Correlation between {prs_column} and {phenotype_column}: R-squared: {r_squared:.4f}\n")
            f.write(f"Correlation between {prs_column} and {phenotype_column}: P-value: {p_value:.4e}\n")
    print(f"Correlation result saved to: {correlation_file}")

    # Plot scatter plot
    plt.figure(figsize=(8, 6))
    plt.scatter(df_prs[prs_column], df_prs[phenotype_column], alpha=0.5, edgecolor='k', label="Data points")

    # Add a regression line for continuous phenotypes
    if not is_binary:
        x = df_prs[prs_column]
        y = df_prs[phenotype_column]
        m, b = np.polyfit(x, y, 1)  # Linear regression
        plt.plot(x, m * x + b, color='red', label=f"Regression line (y={m:.2f}x + {b:.2f})")

    # Add title and labels
    plt.title(f"{prs_column} vs {phenotype_column}\n{correlation_metric}")
    plt.xlabel(prs_column)
    plt.ylabel(phenotype_column)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.legend()

    # Save the scatter plot to a file in the output directory
    scatter_plot_file = os.path.join(output_dir, f"{prs_column}_vs_{phenotype_column}_scatterplot.png")
    plt.savefig(scatter_plot_file, dpi=300, bbox_inches='tight')
    print(f"Scatter plot saved to: {scatter_plot_file}")
    plt.show()

def main():
    # Command-line argument parser
    parser = argparse.ArgumentParser(description="Calculate correlation between a PRS column and a phenotype, and create a scatter plot.")
    
    # Define the argument for the PRS data file
    parser.add_argument(
        "-prs", "--prs_file", required=True, help="Path to the PRS data file (must include the PRS and phenotype columns)."
    )
    
    # Define the argument for the PRS column name
    parser.add_argument(
        "-prs_col", "--prs_column", required=True, help="Name of the column containing PRS values."
    )
    
    # Define the argument for the phenotype column name
    parser.add_argument(
        "-pheno", "--phenotype_column", required=True, help="Name of the phenotype column to correlate with the PRS column."
    )

    # Define the argument for the output directory
    parser.add_argument(
        "-out", "--output_dir", required=True, help="Directory to save the output files (correlation result and scatter plot)."
    )
    
    # Parse the arguments
    args = parser.parse_args()

    # Call the calculate_correlation_and_plot function with the parsed arguments
    calculate_correlation_and_plot(args.prs_file, args.prs_column, args.phenotype_column, args.output_dir)

if __name__ == "__main__":
    main()
