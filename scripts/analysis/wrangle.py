import pandas as pd
from itertools import compress
import numpy as np

codon_aa_dict = {
    "GCT": "Ala", "GCC": "Ala", "GCA": "Ala", "GCG": "Ala",
    "CGT": "Arg", "CGC": "Arg", "CGA": "Arg", "CGG": "Arg", "AGA": "Arg", "AGG": "Arg",
    "GGT": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly",
    "TTA": "Leu", "TTG": "Leu", "CTT": "Leu", "CTC": "Leu", "CTA": "Leu", "CTG": "Leu",
    "CCT": "Pro", "CCC": "Pro", "CCA": "Pro", "CCG": "Pro",
    "TCT": "Ser", "TCC": "Ser", "TCA": "Ser", "TCG": "Ser", "AGT": "Ser", "AGC": "Ser",
    "ACT": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr",
    "GTT": "Val", "GTC": "Val", "GTA": "Val", "GTG": "Val",
    "TTT": "Phe", "TTC": "Phe",
    "TAT": "Tyr", "TAC": "Tyr",
    "CAT": "His", "CAC": "His",
    "CAA": "Gln", "CAG": "Gln",
    "AAT": "Asn", "AAC": "Asn",
    "AAA": "Lys", "AAG": "Lys",
    "GAT": "Asp", "GAC": "Asp",
    "GAA": "Glu", "GAG": "Glu",
    "TGT": "Cys", "TGC": "Cys",
    "ATT": "Ile", "ATC": "Ile", "ATA": "Ile",
    "ATG": "Met",
    "TGG": "Trp",
    "TAA": "Stop", "TGA": "Stop","TAG": "Stop"
}

#### No missing polymorphism #####

def state_founder(df, M):
    states = []
    for _, row in df.iterrows(): # df.iterrows() yields both the index and row (as Series)
        if row['count_N'] != 0:
            states.append("Miss")
        elif row['count_A'] == M or row['count_T'] == M or row['count_G'] == M or row['count_C'] == M:
            states.append("M")
        else:
            states.append("P")
    return states


##### Filter biallelic ####
def filter_multiallele(df):
    indices_to_remove = []
    for i, row in df.iterrows():
        if (row[2:6] == 0).sum() < 2:
            indices_to_remove.append(i)
    return df.drop(indices_to_remove)


##### Polymorphism type define ####
## Type I: Ts vs TV
def poltype_i(df):
    bases = ["A", "T", "G", "C"]
    transition = {"AG", "GA", "TC", "CT"}
    transversion = {"AC", "CA", "AT", "TA", "CG", "GC", "GT", "TG"}

    def determine_poltype(row):
        if row['State'] == "M":
            return "M"
        else:
            present_bases = [bases[i] for i, count in enumerate(row[['count_A', 'count_T', 'count_G', 'count_C']]) if
                             count != 0]
            if len(present_bases) >= 2:
                tmp = present_bases[0] + present_bases[1]
                if tmp in transition:
                    return "Ts"
                else:
                    return "Tv"
            return "M"

    df['Poltype_I'] = df.apply(determine_poltype, axis=1)
    return df

## Type II: GCn vs GCc

def poltype_ii(df):
    bases = ["A", "T", "G", "C"]
    GCneutral = {"AT", "TA", "CG", "GC"}
    GCchanging = {"AC","CA","AG","GA","CT","TC","GT","TG"}

    def determine_poltype(row):
        if row['State'] == "M":
            return "M"
        else:
            present_bases = [bases[i] for i, count in enumerate(row[['count_A', 'count_T', 'count_G', 'count_C']]) if
                             count != 0]
            if len(present_bases) >= 2:
                tmp = present_bases[0] + present_bases[1]
                if tmp in GCneutral:
                    return "GCn"
                else:
                    return "GCc"
            return "M"

    df['Poltype_II'] = df.apply(determine_poltype, axis=1)
    return df

## Type III: Nonsyn vs Syn


def poltype_iii(df, codon_aa_dict):
    bases = ["A", "T", "G", "C"]
    count_cols = ['count_A', 'count_T', 'count_G', 'count_C']

    PoltypeIII = pd.Series(["M"] * len(df), index=df.index)

    PoltypeIII[df[(df['Degeneracy'] == "0fold") & (df['State'] == "P")].index] = "Nonsyn"
    PoltypeIII[df[(df['Degeneracy'] == "4fold") & (df['State'] == "P")].index] = "Syn"
    PoltypeIII[df[(df['Degeneracy'] == "stop") & (df['State'] == "P")].index] = "Stop_pol"

    idx2fold = df[(df['Degeneracy'] == "2fold") & (df['State'] == "P")].index

    for i in idx2fold:
        codon_vec = list(df.at[i, 'codon'])
        codon_pos = df.at[i, 'codon_pos'] - 1  # Adjusting for zero-based indexing
        pols = [bases[j] for j, count in enumerate(df.loc[i, count_cols]) if count != 0]

        change1 = codon_vec.copy()
        change2 = codon_vec.copy()

        change1[codon_pos] = pols[0]
        change2[codon_pos] = pols[1]

        res = (codon_aa_dict["".join(change1)] == codon_aa_dict["".join(change2)])
        PoltypeIII[i] = "Syn" if res else "Nonsyn"


    idx3fold = df[(df['Degeneracy'] == "3fold") & (df['State'] == "P")].index

    for i in idx3fold:
        if df.at[i, 'count_G'] == 0:
            PoltypeIII[i] = "Syn"
        else:
            PoltypeIII[i] = "Nonsyn"

    df['Poltype_III'] = PoltypeIII
    return df


### Codon change polymorphism ####

def find_codon_change(df):
    bases = ["A", "T", "G", "C"]
    count_cols = ['count_A', 'count_T', 'count_G', 'count_C']

    codon_change = df['codon'].copy()
    poly = df[df['State'] == "P"].index  # only polymorphic rows

    for i in poly:
        codon_vec = list(df.at[i, 'codon'])
        codon_pos = df.at[i, 'codon_pos'] - 1  # convert 1-based to 0-based
        ref = codon_vec[codon_pos]

        # Get bases with nonzero counts
        pols = list(compress(bases, df.loc[i, count_cols] != 0))

        # Derived = polymorphic bases that are not the reference base
        der = [base for base in pols if base != ref]

        if len(der) == 1:
            codon_vec[codon_pos] = der[0]
            codon_change[i] = ''.join(codon_vec)
        elif len(der) > 1:
            codon_vec_c1 = codon_vec.copy()
            codon_vec_c2 = codon_vec.copy()
            codon_vec_c1[codon_pos] = der[0]
            codon_vec_c2[codon_pos] = der[1]
            codon_change[i] = ','.join([''.join(codon_vec_c1), ''.join(codon_vec_c2)])
        # if no derived allele found, leave codon_change[i] as original codon

    df['codon_change'] = codon_change
    return df


#### Amino acid change polymorphism ####

def find_aa_change(df, codon_aa_dict):
    mult_idx = df[df['codon_change'].str.contains(",")].index
    sing_idx = df[~df['codon_change'].str.contains(",")].index

    aa_change = df['aa'].copy()

    for i in sing_idx:
        tmp = df.at[i, 'codon_change']
        aa_change.at[i] = codon_aa_dict[tmp]

    for i in mult_idx:
        tmp = df.at[i, 'codon_change'].split(',')
        aa_change.at[i] = ','.join([codon_aa_dict[codon] for codon in tmp])

    df['aa_change'] = aa_change
    return df

##### Create matrices ####

def create_codon_matrix_nocomplex(df):
    codons = ["AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT",
              "ATA", "ATC", "ATG", "ATT", "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT",
              "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT", "GAA", "GAC", "GAG", "GAT",
              "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT",
              "TAA", "TAC", "TAG", "TAT", "TCA", "TCC", "TCG", "TCT", "TGA", "TGC", "TGG", "TGT",
              "TTA", "TTC", "TTG", "TTT"]

    codon_matrix = pd.DataFrame(0, index=codons, columns=codons)

    for i in range(len(df)):
        if "," not in df.at[i, 'codon_change']:
            ref_codon = df.at[i, 'codon']
            change_codon = df.at[i, 'codon_change']
            codon_matrix.at[ref_codon, change_codon] += 1

    return codon_matrix


def create_aa_matrix_nocomplex(df):
    aas_3 = ["Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile",
             "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val", "Stop"]

    aa_matrix = pd.DataFrame(0, index=aas_3, columns=aas_3)

    for i in range(len(df)):
        if "," not in df.at[i, 'aa_change']:
            ref_aa = df.at[i, 'aa']
            change_aa = df.at[i, 'aa_change']
            aa_matrix.at[ref_aa, change_aa] += 1

    return aa_matrix

##### Add MAC and MAF info #### Not major, minor

def MAC_info(df):
    # Calculate the total count for each row
    total_counts = df[['count_A', 'count_T', 'count_G', 'count_C']].sum(axis=1)

    # Initialize the MAC list
    MAC = []

    # Calculate MAC for each row
    for i in range(len(df)):
        # Get the counts for A, T, G, C
        counts = df.loc[i, ['count_A', 'count_T', 'count_G', 'count_C']].values

        # Check if there are more than two zero counts
        if sum(counts == 0) > 2:
            MAC.append(0)
        else:
            # Sort the counts and take the third value
            MAC.append(sorted(counts)[2])

    # Add MAC and MAF to the DataFrame
    df['MAC'] = MAC
    df['MAF'] = df['MAC'] / total_counts

    return df


######## SFS build functions ######- different from stats.py: return df not list

def unfolded_sfs(df, M):
    """
    Build unfolded SFS from dataframe with column 'DAC' (derived allele count).

    Parameters:
        df (pd.DataFrame): must contain column 'DAC'
        M (int): sample size

    Returns:
        pd.DataFrame: index = 0..M, column = 'count'
    """
    spectra = pd.DataFrame(0, index=np.arange(0, M + 1), columns=["count"])
    for i in range(0, M + 1):
        spectra.loc[i, "count"] = (df["DAC"] == i).sum()
    return spectra


def folded_sfs(df, M):
    """
    Build folded SFS from dataframe with column 'MAC' (minor allele count).
    Handles both even and odd M.

    Parameters:
        df (pd.DataFrame): must contain column 'MAC'
        M (int): sample size

    Returns:
        pd.DataFrame: index = 0..floor(M/2), column = 'count'
    """
    if M % 2 == 0:  # even
        max_i = M // 2
    else:  # odd
        max_i = (M - 1) // 2

    spectra = pd.DataFrame(0, index=np.arange(0, max_i + 1), columns=["count"])
    for i in range(0, max_i + 1):
        spectra.loc[i, "count"] = ((df["MAC"] == i) | (df["MAC"] == (M - i))).sum()
    return spectra