import pandas as pd
import sys

#### One Outgroup #####
# Path to Est-sfs output

est_sfs_file = sys.argv[1] # estsfs pvalues out
count_file = sys.argv[2] #  Population genetics (allele count) data
out_file = sys.argv[3] #
nsamp = sys.argv[4]
ref_col = sys.argv[5]  # name of outgroup column (e.g. Ref_rat, etc.)

# Read the Est-sfs output file, skipping the first 8 header lines
est_sfs_data = pd.read_csv(est_sfs_file, sep=r'\s+', skiprows=8,
                           header=None,
                           names=['line_number', 'config_index', 'P_major_ancestral',
                                  'P_trees'])

count_data = pd.read_csv(count_file, sep="\t")

# Extract probability of major allele being ancestral
count_data['P_major_ancestral'] = est_sfs_data['P_major_ancestral']
count_data['P_tress'] = est_sfs_data['P_trees']

# --- Remove sites with ambiguous ancestral probability or missing outgroup info ---
count_data = count_data[
    (count_data['P_major_ancestral'] <= 0.45) |
    (count_data['P_major_ancestral'] >= 0.55)
]

count_data = count_data[
    ~(count_data[ref_col].isna() | (count_data[ref_col] == "N"))
]


# --- Major_allele ---
bases = ['A', 'T', 'G', 'C']
count_cols = ['count_A', 'count_T', 'count_G', 'count_C']

count_data['Major_allele'] = count_data[count_cols].idxmax(axis=1).str.replace('count_', '')


# --- Derived_allele ---
def define_derived(row):
    if row['State'] == "M":
        return row['Major_allele']
    else:
        # polymorphic bases
        pol_bases = [base for base, count in
                     zip(bases, [row['count_A'], row['count_T'], row['count_G'], row['count_C']]) if count != 0]
        major = row['Major_allele']
        minor_candidates = [b for b in pol_bases if b != major]
        if not minor_candidates:
            return None  # in case only one base present
        minor = minor_candidates[0]

        if row['P_major_ancestral'] > 0.5:
            return minor
        else:
            return major


count_data['Derived_allele'] = count_data.apply(define_derived, axis=1)


# --- DAC ---
def get_dac(row, total_samples):
    if row['State'] == "P":
        return row[f'count_{row["Derived_allele"]}']
    elif row['State'] == "M":
        if pd.isna(row[ref_col]) or row[ref_col] == "N":
            return None  # missing ancestral info
        elif row[ref_col] == row['Major_allele']:
            return 0
        else:
            return total_samples
    else:
        return None


total_samples = nsamp  # same as M in your R function
count_data['DAC'] = count_data.apply(lambda r: get_dac(r, total_samples), axis=1)
count_data['DAF'] = count_data['DAC'] /total_samples   #

# Save
count_data.to_csv(out_file, sep="\t", index=False)

