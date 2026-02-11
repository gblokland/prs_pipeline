#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.linear_model import LogisticRegression, LinearRegression
import argparse
import os
import numpy as np

def plot_and_save(prs_file, prs_column, pheno_column, output_dir):
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Load the PRS data
    try:
        df_prs = pd.read_csv(prs_file)
    except Exception as e:
        print(f"Error reading file: {e}")
        return

    # Set missing value code (-9) to NaN
    df_prs[prs_column] = df_prs[prs_column].replace(-9, np.nan)
    df_prs[pheno_column] = df_prs[pheno_column].replace(-9, np.nan)

    # Check if necessary columns exist
    if prs_column not in df_prs.columns or pheno_column not in df_prs.columns:
        print(f"Error: '{prs_column}' and/or '{pheno_column}' column not found.")
        return

    # Sanitize column names for filenames (remove spaces, special characters, etc.)
    prs_column_sanitized = prs_column.replace(" ", "_").replace("-", "_")
    pheno_column_sanitized = pheno_column.replace(" ", "_").replace("-", "_")

    # Drop rows where either PRS or phenotype is NaN
    df_prs_cleaned = df_prs.dropna(subset=[prs_column, pheno_column])

    # Check if phenotype is binary or continuous
    pheno_unique = df_prs_cleaned[pheno_column].nunique()
    output_text = ""

    if pheno_unique == 2:  # Binary phenotype
        print("Binary phenotype detected.")
        output_text += "Binary phenotype detected.\n"

        # Train logistic regression
        model = LogisticRegression()
        model.fit(df_prs_cleaned[[prs_column]], df_prs_cleaned[pheno_column])

        # Predict probabilities
        probs = model.predict_proba(df_prs_cleaned[[prs_column]])[:, 1]
        auc = roc_auc_score(df_prs_cleaned[pheno_column], probs)
        output_text += f"AUC-ROC: {auc}\n"
        print(f"AUC-ROC: {auc}")

        # Plot ROC curve
        fpr, tpr, thresholds = roc_curve(df_prs_cleaned[pheno_column], probs)
        plt.figure(figsize=(8, 6))
        plt.plot(fpr, tpr, color='b', label=f'ROC curve (AUC = {auc:.2f})')
        plt.plot([0, 1], [0, 1], color='gray', linestyle='--')  # Diagonal line
        plt.title(f"ROC Curve for {prs_column_sanitized} vs {pheno_column_sanitized}")
        plt.xlabel("False Positive Rate")
        plt.ylabel("True Positive Rate")
        plt.legend(loc="lower right")
        plt.grid(alpha=0.5)

        # Save ROC plot
        roc_plot_path = os.path.join(output_dir, f"roc_curve_{prs_column_sanitized}_{pheno_column_sanitized}.png")
        plt.savefig(roc_plot_path)
        plt.close()
        output_text += f"ROC curve saved to {roc_plot_path}\n"

    else:  # Continuous phenotype
        print("Continuous phenotype detected.")
        output_text += "Continuous phenotype detected.\n"

        # Train linear regression
        model = LinearRegression()
        model.fit(df_prs_cleaned[[prs_column]], df_prs_cleaned[pheno_column])

        # Predict values
        preds = model.predict(df_prs_cleaned[[prs_column]])
        r_squared = model.score(df_prs_cleaned[[prs_column]], df_prs_cleaned[pheno_column])
        output_text += f"R-squared: {r_squared:.4f}\n"
        print(f"R-squared: {r_squared:.4f}")

        # Scatter plot with regression line
        plt.figure(figsize=(8, 6))
        plt.scatter(df_prs_cleaned[prs_column], df_prs_cleaned[pheno_column], alpha=0.6, label="Data points")
        plt.plot(df_prs_cleaned[prs_column], preds, color="red", label=f"Regression line ($R^2$ = {r_squared:.2f})")
        plt.title(f"{prs_column_sanitized} vs {pheno_column_sanitized}")
        plt.xlabel(prs_column)
        plt.ylabel(pheno_column)
        plt.legend()
        plt.grid(alpha=0.5)

        # Save scatter plot
        scatter_plot_path = os.path.join(output_dir, f"scatter_plot_{prs_column_sanitized}_{pheno_column_sanitized}.png")
        plt.savefig(scatter_plot_path)
        plt.close()
        output_text += f"Scatter plot saved to {scatter_plot_path}\n"

    # Save output text to a file
    output_text_path = os.path.join(output_dir, f"results_{prs_column_sanitized}_{pheno_column_sanitized}.txt")
    with open(output_text_path, "w") as f:
        f.write(output_text)
    print(f"Results saved to {output_text_path}")

def main():
    # Command-line argument parser
    parser = argparse.ArgumentParser(description="Train a model and plot/save ROC-AUC or regression for PRS prediction.")
    
    # Define the argument for the PRS file
    parser.add_argument(
        "-file", "--prs_file", required=True, help="Path to the PRS data file."
    )
    
    # Define the argument for the PRS column name
    parser.add_argument(
        "-prs", "--prs_column", required=True, help="Name of the PRS column."
    )
    
    # Define the argument for the phenotype column name
    parser.add_argument(
        "-pheno", "--pheno_column", required=True, help="Name of the phenotype column."
    )

    # Define the argument for the output directory
    parser.add_argument(
        "-out", "--output_dir", required=True, help="Directory to save plots and results."
    )

    # Parse the arguments
    args = parser.parse_args()

    # Call the function with parsed arguments
    plot_and_save(args.prs_file, args.prs_column, args.pheno_column, args.output_dir)

if __name__ == "__main__":
    main()
