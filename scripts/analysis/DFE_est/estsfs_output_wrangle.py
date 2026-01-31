import pandas as pd
import sys


# Path to Est-sfs output

est_sfs_file = sys.argv[1] # pvalues output
outfile=sys.argv[2]

# Read the Est-sfs output file, skipping the first 8 header lines
est_sfs_data = pd.read_csv(est_sfs_file, sep=r'\s+', skiprows=8,
                           header=None,
                           names=['line_number', 'config_index', 'P_major_ancestral',
                                  'P_A', 'P_C', 'P_G', 'P_T'])


# Extract probability of major allele being ancestral
p_major_ancestral = est_sfs_data['P_major_ancestral']

# Extract the ancestral nucleotide with the highest probability
nucleotide_probs = est_sfs_data[['P_A', 'P_C', 'P_G', 'P_T']]
ancestral_base = nucleotide_probs.idxmax(axis=1).str[-1]  # Get the last character from 'P_X' => 'X'

# Combine into a new DataFrame
output_df = pd.DataFrame({
    'P_major_ancestral': p_major_ancestral,
    'Ancestral_base': ancestral_base
})


output_df.to_csv(outfile, sep='\t', index=False)