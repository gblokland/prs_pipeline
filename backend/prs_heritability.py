import pandas as pd
import statsmodels.api as sm
import argparse
import os

def calculate_prs_heritability(input_file, phenotype_col, prs_col, covariate_cols=None, delimiter=',', output_file=None):
    """
    Calculate heritability explained (R^2) by a Polygenic Risk Score (PRS).
    
    Args:
        input_file (str): Path to the input TSV/CSV file.
        phenotype_col (str): Column name for the phenotype (continuous).
        prs_col (str): Column name for the PRS.
        covariate_cols (list): List of column names for covariates (optional).
        delimiter (str): Delimiter for the input file (default: ',').
        output_file (str): Path to save results (optional, without extension).
    """
    # Load data
    print("Loading input file...")
    data = pd.read_csv(input_file, delimiter=delimiter)
    
    # Filter data: Keep only rows where Phase == 1
    if 'Phase' in data.columns:
        print("Filtering data to keep only rows where Phase == 1...")
        data = data[data['Phase'] == 1]
        print(f"Filtered data has {data.shape[0]} rows remaining.")
    else:
        print("Warning: Column 'Phase' not found. Skipping filtering.")

    # Check that required columns exist
    required_cols = [phenotype_col, prs_col]
    if covariate_cols:
        required_cols += covariate_cols
    missing_cols = [col for col in required_cols if col not in data.columns]
    if missing_cols:
        raise ValueError(f"Missing columns in input file: {', '.join(missing_cols)}")
    
    # Prepare PRS column
    X = data[[prs_col]]
    
    # Handle covariates: Convert integers to dummies
    if covariate_cols:
        print("Processing covariates...")
        covariates = data[covariate_cols]
        for col in covariate_cols:
            if pd.api.types.is_integer_dtype(data[col]):
                print(f"Converting integer covariate '{col}' to dummy variables...")
                dummies = pd.get_dummies(data[col], prefix=col, drop_first=True)
                covariates = covariates.drop(columns=col)
                covariates = pd.concat([covariates, dummies], axis=1)
        X = pd.concat([X, covariates], axis=1)

    # Add intercept
    X = sm.add_constant(X)

    # Prepare phenotype
    y = data[phenotype_col]

    # Fit linear regression
    print("Fitting linear regression model...")
    model = sm.OLS(y, X).fit()
    
    # Extract R-squared
    r2 = model.rsquared
    print(f"Heritability explained by PRS (R^2): {r2:.4f}")
    
    # Save results if output_file is specified
    if output_file:
        txt_output = output_file + ".txt"
        csv_output = output_file + ".csv"

        print(f"Saving results to {txt_output} and {csv_output}...")
        
        # Write to text file
        with open(txt_output, 'w') as f:
            f.write(f"Heritability explained (R^2): {r2:.4f}\n")
        
        # Write to CSV file
        output_data = {
            "Phenotype Column": phenotype_col,
            "PRS Column": prs_col,
            "Covariates": ', '.join(covariate_cols) if covariate_cols else "None",
            "Heritability Explained (R^2)": r2
        }
        output_df = pd.DataFrame([output_data])
        output_df.to_csv(csv_output, index=False)
        
        print("Results saved successfully.")

    # Print model summary (optional for debugging)
    print("\nModel Summary:")
    print(model.summary())
    
    return r2

if __name__ == "__main__":
    # Command-line argument parser
    parser = argparse.ArgumentParser(description="Calculate heritability explained (R^2) by PRS using linear regression.")
    
    parser.add_argument("input_file", help="Path to input file (TSV or CSV).")
    parser.add_argument("phenotype_col", help="Column name for the phenotype.")
    parser.add_argument("prs_col", help="Column name for the PRS.")
    parser.add_argument("--covariates", nargs='*', help="List of covariate column names (optional).", default=None)
    parser.add_argument("--delimiter", help="Delimiter for the input file (default: ','). Use '\\t' for TSV files.", default=',')
    parser.add_argument("--output_file", help="Base path to save results (without extension).", default=None)
    
    # Parse arguments
    args = parser.parse_args()
    
    # Run the function
    try:
        calculate_prs_heritability(
            input_file=args.input_file,
            phenotype_col=args.phenotype_col,
            prs_col=args.prs_col,
            covariate_cols=args.covariates,
            delimiter=args.delimiter,
            output_file=args.output_file
        )
    except Exception as e:
        print(f"Error: {e}")
