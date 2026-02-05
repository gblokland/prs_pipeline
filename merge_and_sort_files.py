import pandas as pd
import argparse

def merge_files(file1, file2, output_file):
    # Read the two files as dataframes
    df1 = pd.read_csv(file1, sep='\t')  # Assuming tab-separated input files
    df2 = pd.read_csv(file2, sep='\t')
    
    # Merge the dataframes on 'CHROM' and 'POS' columns
    merged_df = pd.merge(df1, df2, on=['CHROM', 'POS'], how='inner')  # 'inner' to keep common rows

    # Sort the merged dataframe by the 'SNP' column
    merged_df_sorted = merged_df.sort_values(by='SNP')

    # Save the sorted, merged dataframe to the output file
    merged_df_sorted.to_csv(output_file, sep='\t', index=False)

def main():
    # Set up command line argument parsing
    parser = argparse.ArgumentParser(description="Merge two files on CHROM and POS columns, and sort by SNP.")
    parser.add_argument("file1", help="Path to the first input file")
    parser.add_argument("file2", help="Path to the second input file")
    parser.add_argument("output_file", help="Path to the output merged and sorted file")

    # Parse the command line arguments
    args = parser.parse_args()
    
    # Merge the files, sort by SNP, and save the result
    merge_files(args.file1, args.file2, args.output_file)

if __name__ == "__main__":
    main()
