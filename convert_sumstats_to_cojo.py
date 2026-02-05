import pandas as pd
import sys

def convert_sumstats_to_cojo(gwas_file, freq_file, output_file, chr_col_options=None, bp_col_options=None, refall_col_options=None, altall_col_options=None, freq_col_options=None, beta_col_options=None, se_col_options=None, pval_col_options=None, n_col_options=None):
    # Load GWAS summary statistics
    gwas_data = pd.read_csv(gwas_file, sep='\t')
    
    # Load allele frequency data
    freq_data = pd.read_csv(freq_file, sep='\t')
    
    #######################################################
    # Default options for CHR column names if not provided
    if chr_col_options is None:
        chr_col_options = ['CHROM', 'Chr', 'CHR', 'chromosome']

    # Find the first matching frequency column from the options
    chr_col = None
    for col in chr_col_options:
        if col in gwas_data.columns:
            chr_col = col
            break

    if chr_col is None:
        raise ValueError(f"None of the provided chromosome columns {chr_col_options} were found in the gwas file.")

    #######################################################
    # Default options for BP column names if not provided
    if bp_col_options is None:
        bp_col_options = ['BP', 'POS', 'BPOS', 'base_pair_location']

    # Find the first matching frequency column from the options
    bp_col = None
    for col in bp_col_options:
        if col in gwas_data.columns:
            bp_col = col
            break

    if bp_col is None:
        raise ValueError(f"None of the provided BP columns {bp_col_options} were found in the gwas file.")

    #######################################################
    # Default options for A1 column names if not provided
    if refall_col_options is None:
        refall_col_options = ['A1', 'RefAllele', 'Allele1', 'EA', 'EffectAllele', 'effect_allele']

    # Find the first matching A1 column from the options
    refall_col = None
    for col in refall_col_options:
        if col in gwas_data.columns:
            refall_col = col
            break

    if refall_col is None:
        raise ValueError(f"None of the provided A1 columns {refall_col_options} were found in the gwas file.")

    #######################################################
    # Default options for A2 column names if not provided
    if altall_col_options is None:
        altall_col_options = ['A2', 'AltAllele', 'Allele2', 'other_allele']

    # Find the first matching A2 column from the options
    altall_col = None
    for col in altall_col_options:
        if col in gwas_data.columns:
            altall_col = col
            break

    if altall_col is None:
        raise ValueError(f"None of the provided A2 columns {altall_col_options} were found in the gwas file.")

    #######################################################
    # Default options for BETA column names if not provided
    if beta_col_options is None:
        beta_col_options = ['BETA', 'beta', 'B']

    # Find the first matching BETA column from the options
    beta_col = None
    for col in beta_col_options:
        if col in gwas_data.columns:
            beta_col = col
            break

    if beta_col is None:
        raise ValueError(f"None of the provided BETA columns {beta_col_options} were found in the gwas file.")

    #######################################################
    # Default options for SE column names if not provided
    if se_col_options is None:
        se_col_options = ['SE', 'STDERR', 'standard_error']

    # Find the first matching SE column from the options
    se_col = None
    for col in se_col_options:
        if col in gwas_data.columns:
            se_col = col
            break

    if se_col is None:
        raise ValueError(f"None of the provided SE columns {se_col_options} were found in the gwas file.")

    #######################################################
    # Default options for P column names if not provided
    if pval_col_options is None:
        pval_col_options = ['PVAL', 'PVALUE', 'P', 'p_value']

    # Find the first matching pval column from the options
    pval_col = None
    for col in pval_col_options:
        if col in gwas_data.columns:
            pval_col = col
            break

    if pval_col is None:
        raise ValueError(f"None of the provided PVALUE columns {pval_col_options} were found in the gwas file.")

    #######################################################
    # Default options for N column names if not provided
    if n_col_options is None:
        n_col_options = ['N', 'NMISS', 'SAMPLESIZE', 'SAMPLE_SIZE', 'sample_size']

    # Find the first matching N column from the options
    n_col = None
    for col in n_col_options:
        if col in gwas_data.columns:
            n_col = col
            break

    if n_col is None:
        raise ValueError(f"None of the provided N columns {n_col_options} were found in the gwas file.")

    # Rename columns to match frequency file format
    gwas_data.rename(columns={
        chr_col: 'CHROM',
        bp_col: 'POS',
    }, inplace=True)

    #######################################################
    # Default options for Freq column names if not provided
    if freq_col_options is None:
        freq_col_options = ['Freq1', 'Freq', 'FREQ', 'MAF', 'AlleleFreq', 'AlleleFreq1', 'AF', 'effect_allele_frequency', 'freq']

#    # Find the first matching frequency column from the options
#    freq_col = None
#    for col in freq_col_options:
#        if col in gwas_data.columns:
#            freq_col = col
#            break
#
#    if freq_col is None:
#        raise ValueError(f"None of the provided frequency columns {freq_col_options} were found in the gwas file.")


    #######################################################
    # Merge GWAS and frequency data on SNP ID or CHROM-POS
    # Check if any of the frequency columns exist in the GWAS file
    freq_col = None
    for col in gwas_data.columns:
        if col.lower() in [fc.lower() for fc in freq_col_options]:
            freq_col = col
            break

    # If no frequency column is found, merge with the freq_file
    if not freq_col:
        print("No frequency column found in the GWAS file. Merging with the frequency file...")

        # Perform the merge based on a common column (for example: 'CHR' and 'BP')
        # Adjust the column names ('CHR', 'BP') to match the actual columns in your files
        merged_data = pd.merge(gwas_data, freq_data, how='left', on=['CHROM','POS'])

        # After merging, you can use the frequency column from the freq file
        if 'Freq1' in merged_data.columns:
            freq_col = 'Freq1'
        else:
            raise ValueError("Frequency column not found in the merged file either.")

    else:
        print(f"Frequency column found in the GWAS file: {freq_col}")
        merged_data = gwas_data  # No need to merge, just use the original GWAS file
        # At this point, `merged_df` contains the necessary frequency data, either from the GWAS file or from the merged freq_file.

    #######################################################
    # Rename columns to match COJO format, using the identified frequency column
    merged_data.rename(columns={
        'SNP': 'SNP',
        chr_col: 'CHR',
        bp_col: 'BP',
        refall_col: 'A1',
        altall_col: 'A2',
        freq_col: 'freq',  # Use the found frequency column
        beta_col: 'b',
        se_col: 'se',
        pval_col: 'p',
        n_col: 'N'
    }, inplace=True)

    # Reorder columns and rename them for COJO format
    cojo_data = merged_data[['SNP', 'A1', 'A2', 'freq', 'b', 'se', 'p', 'N']]

    # Save to new file
    cojo_data.to_csv(output_file, sep='\t', index=False)
    print(f"COJO-formatted file saved to {output_file}")
    
if __name__ == "__main__":
    #if len(sys.argv) < 4:
    #    print("Usage: python convert_sumstats_to_cojo.py <gwas_file> <freq_file> <output_file> <chr_col_options> <bp_col_options> <refall_col_options> <altall_col_options> <freq_col_options> <beta_col_options> <se_col_options> <pval_col_options> <n_col_options>")
    #    sys.exit(1)
    
    #gwas_file = sys.argv[1]
    #freq_file = sys.argv[2]
    #output_file = sys.argv[3]
    ## Optionally pass a list of possible frequency column names
    #chr_col_options = sys.argv[4] #['Freq1', 'Frequency', 'AF']  # Add more as needed
    #bp_col_options = sys.argv[5] 
    #refall_col_options = sys.argv[6] 
    #altall_col_options = sys.argv[7] 
    #freq_col_options = sys.argv[8] 
    #beta_col_options = sys.argv[9] 
    #se_col_options = sys.argv[10] 
    #pval_col_options = sys.argv[11] 
    #n_col_options = sys.argv[12] 

    # Mandatory arguments (adjust according to your needs)
    required_args = 4  # You expect at least 4 arguments, e.g., sys.argv[0] to sys.argv[3]

    if len(sys.argv) < required_args:
        print(f"Usage: {sys.argv[0]} <gwas_file> <freq_file> <output_file> [<chr_col_options>] [<bp_col_options>] [<refall_col_options>] [<altall_col_options>] [<freq_col_options>] [<beta_col_options>] [<se_col_options>] [<pval_col_options>] [<n_col_options>]")
        sys.exit(1)

    # Mandatory arguments
    gwas_file = sys.argv[1]
    freq_file = sys.argv[2]
    output_file = sys.argv[3]

    # Handle optional arguments - Optionally pass a list of possible frequency column names
    optional_args = sys.argv[4:]  # Everything after sys.argv[3] is optional

    if len(optional_args) > 0:
        chr_col_options = optional_args[0]  # First optional argument
        print(f"Optional argument 1: {chr_col_options}")
    else:
        print("No optional arguments provided.")
        chr_col_options = None


    # You can also handle more optional arguments similarly
    if len(optional_args) > 1:
        bp_col_options = optional_args[1]
        print(f"Optional argument 2: {bp_col_options}")
    else:
        bp_col_options = None
        
    # You can also handle more optional arguments similarly
    if len(optional_args) > 2:
        refall_col_options = optional_args[2]
        print(f"Optional argument 3: {refall_col_options}")
    else:
        refall_col_options = None

    # You can also handle more optional arguments similarly
    if len(optional_args) > 3:
        altall_col_options = optional_args[3]
        print(f"Optional argument 4: {altall_col_options}")
    else:
        altall_col_options = None

    # You can also handle more optional arguments similarly
    if len(optional_args) > 4:
        freq_col_options = optional_args[4]
        print(f"Optional argument 5: {freq_col_options}")
    else:
        freq_col_options = None

    # You can also handle more optional arguments similarly
    if len(optional_args) > 5:
        beta_col_options = optional_args[5]
        print(f"Optional argument 6: {beta_col_options}")
    else:
        beta_col_options = None

    # You can also handle more optional arguments similarly
    if len(optional_args) > 6:
        se_col_options = optional_args[6]
        print(f"Optional argument 7: {se_col_options}")
    else:
        se_col_options = None

    # You can also handle more optional arguments similarly
    if len(optional_args) > 7:
        pval_col_options = optional_args[7]
        print(f"Optional argument 8: {pval_col_options}")
    else:
        pval_col_options = None

    # You can also handle more optional arguments similarly
    if len(optional_args) > 8:
        n_col_options = optional_args[8]
        print(f"Optional argument 9: {n_col_options}")
    else:
        n_col_options = None

    convert_sumstats_to_cojo(gwas_file, freq_file, output_file, chr_col_options, bp_col_options, refall_col_options, altall_col_options, freq_col_options, beta_col_options, se_col_options, pval_col_options, n_col_options)
