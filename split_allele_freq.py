import pandas as pd
import sys

def split_allele_freq(input_file, output_file):
    # Read the input file into a pandas DataFrame, assuming tab-separated values
    df = pd.read_csv(input_file, sep='\t')

    # Function to split {ALLELE:FREQ} pairs into individual columns
    def process_allele_freq(column_data):
        # Initialize default values for allele and frequency
        allele, freq = '', '0'

        # Process the allele:frequency pair if it exists
        if ':' in column_data:
            parts = column_data.split(':')
            if len(parts) >= 2:
                allele, freq = parts[0], parts[1]
            else:
                # Handle cases with unexpected format
                print(f"Unexpected data format: {column_data}")
        
        return pd.Series([allele, freq])

    # Apply the function to each {ALLELE:FREQ} column to extract alleles and frequencies
    df[['Allele1', 'Freq1']] = df['{ALLELE:FREQ}.1'].apply(process_allele_freq)
    df[['Allele2', 'Freq2']] = df['{ALLELE:FREQ}.2'].apply(process_allele_freq)

    # Convert frequencies to numeric, handle any non-numeric values by setting them to 0
    df['Freq1'] = pd.to_numeric(df['Freq1'], errors='coerce').fillna(0)
    df['Freq2'] = pd.to_numeric(df['Freq2'], errors='coerce').fillna(0)

    # Drop the original '{ALLELE:FREQ}' columns as they are now split
    df = df.drop(columns=['{ALLELE:FREQ}.1', '{ALLELE:FREQ}.2'])
    # Drop the columns 'N_ALLELES' and 'N_CHR'
    df = df.drop(columns=['N_ALLELES', 'N_CHR'])

    # Save the processed DataFrame to the output file
    df.to_csv(output_file, sep='\t', index=False)

# Example usage
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python split_allele_freq.py <input_file> <output_file>")
        sys.exit(1)
    input_file = sys.argv[1]  # First argument: input file
    output_file = sys.argv[2]  # Second argument: output file
    
    split_allele_freq(input_file, output_file)
