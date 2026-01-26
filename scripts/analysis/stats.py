import pandas as pd
import numpy as np

########## Spectra build
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
    spectra = spectra["count"].tolist()
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

    spectra = spectra["count"].tolist()
    return spectra

############## PopGen stats

def PiCalculate_uSFS(SFS, ngenomes):
    numerator = 0.0
    for i in range(1, ngenomes):  # skip i = 0 and i = n (monomorphic)
        count_i = float(SFS[i])
        numerator += count_i * i * (ngenomes - i)

    total_sites = sum(float(x) for x in SFS) # to be sure, in case not integers
    pi = (2 * numerator) / (ngenomes * (ngenomes - 1) * total_sites)
    return pi

def PiCalculate_fSFS(SFS, ngenomes):
    numerator = 0.0
    for i in range(1, len(SFS)):  # Start from 1 to exclude the monomorphic sites
        count_i = float(SFS[i])
        numerator += count_i * i * (ngenomes - i)

    total_sites = sum(float(x) for x in SFS)
    pi = (2 * numerator) / (ngenomes * (ngenomes - 1) * total_sites)
    return pi

#############

def ThetaCalculate_uSFS(SFS, ngenomes):
    theta_sum = 0.0
    for i in range(1, ngenomes):  # skip monomorphic bins at i=0 and i=n
        count_i = float(SFS[i])
        theta_sum += count_i / i
    total_sites = sum(float(x) for x in SFS)
    return theta_sum / total_sites

def ThetaCalculate_fSFS(SFS):
    theta_sum = 0.0
    for i in range(1, len(SFS)):  # skip monomorphic bin at index 0
        count_i = float(SFS[i])
        theta_sum += count_i / i
    total_sites = sum(float(x) for x in SFS)
    return theta_sum / total_sites


##############  Bootstrap: same sample size (sample sites)
def bootstrap_pi(df, M, ngenomes, n_bootstrap=1000, ci=95,folded=False):
    """
    df: DataFrame with your full data
    M: sample size
    ngenomes: number of haploid genomes (basically M for Dmel)
    folded: whether to compute folded SFS (True) or unfolded SFS (False)
    n_bootstrap: number of bootstrap replicates
    ci: confidence interval (e.g., 95 for 95% CI)
    """
    pi_values = []

    for _ in range(n_bootstrap):
        resampled_df = df.sample(n=len(df), replace=True)
        if folded:
            sfs = folded_sfs(resampled_df, M)
            pi = PiCalculate_fSFS(sfs, ngenomes)
        else:
            sfs = unfolded_sfs(resampled_df, M)
            pi = PiCalculate_uSFS(sfs, ngenomes)
        pi_values.append(pi)

    pi_array = np.array(pi_values)
    lower = np.percentile(pi_array, (100 - ci) / 2)
    upper = np.percentile(pi_array, 100 - (100 - ci) / 2)
    mean_pi = np.mean(pi_array)

    return mean_pi, lower, upper

##############  Bootstrap: sample genes ??




##############  Bootstrap: same sample size (sample sites)
def bootstrap_pnps(df_neut,df_sel, M, ngenomes, n_bootstrap=1000, ci=95,folded=False):
    """
    df_neut and df_sel: DataFrames with your full data
    M: sample size
    ngenomes: number of haploid genomes (basically M for Dmel)
    folded: whether to compute folded SFS (True) or unfolded SFS (False)
    n_bootstrap: number of bootstrap replicates
    ci: confidence interval (e.g., 95 for 95% CI)
    """
    pnps_values = []

    for _ in range(n_bootstrap):
        resampled_neut = df_neut.sample(n=len(df_neut), replace=True)
        resampled_sel =  df_sel.sample(n=len(df_sel), replace=True)
        if folded:
            sfs_neut = folded_sfs(resampled_neut, M)
            sfs_sel = folded_sfs(resampled_sel, M)

            pi_neut = PiCalculate_fSFS(sfs_neut, ngenomes)
            pi_sel = PiCalculate_uSFS(sfs_sel, ngenomes)
            pinps=pi_sel / pi_neut
        else:
            sfs_neut = unfolded_sfs(resampled_neut, M)
            sfs_sel = unfolded_sfs(resampled_sel, M)

            pi_neut = PiCalculate_uSFS(sfs_neut, ngenomes)
            pi_sel = PiCalculate_uSFS(sfs_sel, ngenomes)
            pinps = pi_sel / pi_neut

        pnps_values.append(pinps)

    pnps_array = np.array(pnps_values)
    lower = np.percentile(pnps_array, (100 - ci) / 2)
    upper = np.percentile(pnps_array, 100 - (100 - ci) / 2)
    mean = np.mean(pnps_array)

    return mean, lower, upper
