import pandas as pd
import argparse

def merge_prscs_with_covariates_and_phenotype(prscs_file, covariate_file, phenotype_file, output_file):
    try:
        # Load the PRScs, covariate, and phenotype files into pandas DataFrames
        df_prscs = pd.read_csv(prscs_file, sep='\t')
        df_covariates = pd.read_csv(covariate_file)
        df_phenotype = pd.read_csv(phenotype_file, sep='\t')
        
        # Print column names for debugging
        print("PRScs columns:", df_prscs.columns)
        print("Covariate columns:", df_covariates.columns)
        print("Phenotype columns:", df_phenotype.columns)

    except Exception as e:
        print(f"Error reading files: {e}")
        return

    # Ensure the required columns exist in all DataFrames
    required_columns = ['FID', 'IID']
    for df, name in zip([df_prscs, df_covariates, df_phenotype], ['PRScs', 'Covariates', 'Phenotype']):
        for col in required_columns:
            if col not in df.columns:
                raise KeyError(f"Column '{col}' is missing in {name} data.")

    # Ensure consistent data types for merging columns
    for df in [df_prscs, df_covariates, df_phenotype]:
        df['FID'] = df['FID'].astype(str)
        df['IID'] = df['IID'].astype(str)

    # Merge the PRScs data with the covariate data
    #df_merged = pd.merge(df_prscs, df_covariates, on=['FID', 'IID'], how='inner')
    df_merged = pd.merge(df_prscs, df_covariates, on=['IID'], how='inner')

    # Merge the result with the phenotype data
    #df_merged = pd.merge(df_merged, df_phenotype, on=['FID', 'IID'], how='inner')
    df_merged = pd.merge(df_merged, df_phenotype, on=['IID'], how='inner')

    # Rename SEX_y to SEX if it exists
    if 'SEX_y' in df_merged.columns:
        df_merged = df_merged.rename(columns={'SEX_y': 'SEX'})

    # Drop unnecessary columns
    columns_to_drop = ['FID_x', 'FID_y', 'PHENO1', 'PHENO', 'SEX_x']
    df_merged = df_merged.drop(columns=[col for col in columns_to_drop if col in df_merged.columns], errors='ignore')

    # Reorder columns to have FID and IID first
    columns_order = ['FID', 'IID', 'Phase'] + [col for col in df_merged.columns if col not in ['FID', 'IID', 'Phase']]
    df_merged = df_merged[columns_order]

    # Sort by FID, IID, and Phase
    sort_columns = ['FID', 'IID', 'Phase']
    df_merged = df_merged.sort_values(by=[col for col in sort_columns if col in df_merged.columns])

    # Save the merged file
    df_merged.to_csv(output_file, index=False) #sep='\t', 
    print(f"Merged data saved to: {output_file}")

def main():
    # Command-line argument parser
    parser = argparse.ArgumentParser(description="Merge PRScs data with covariates and phenotype data.")
    parser.add_argument('-prscs', '--prscs_file', required=True, help="Path to PRScs file")
    parser.add_argument('-cov', '--covariate_file', required=True, help="Path to covariate file")
    parser.add_argument('-pheno', '--phenotype_file', required=True, help="Path to phenotype file")
    parser.add_argument('-o', '--output_file', required=True, help="Path to save the merged output")
    #prscs_file = 'path/to/PRScs_merged_sscore.tsv'
    #covariate_file = 'path/to/covariate_data.csv'
    #phenotype_file = 'path/to/phenotype_data.csv'
    #output_file = 'path/to/output_merged_data.csv'

    args = parser.parse_args()

    # Call the merge function with command-line arguments
    merge_prscs_with_covariates_and_phenotype(args.prscs_file, args.covariate_file, args.phenotype_file, args.output_file)

if __name__ == "__main__":
    main()

