import pandas as pd
import os
import sys
from pathlib import Path
analysis_dir = Path(__file__).resolve().parents[1]
sys.path.append(str(analysis_dir))
import stats
import wrangle

count_file=sys.argv[1] # Population genetics (allele count) data

### Gene names in each bin
bin1=sys.argv[2]
bin2=sys.argv[3]
bin3=sys.argv[4]
bin4=sys.argv[5]

out_file1=sys.argv[6] ## Pi values
out_file2=sys.argv[7] ## P0/P4 values
nsamples=int(sys.argv[8])

##### Read and prepare data
count_df=pd.read_csv(count_file, delimiter="\t")

bin1_genes=pd.read_csv(bin1, sep="\t", header=None)
bin2_genes=pd.read_csv(bin2, sep="\t", header=None)
bin3_genes=pd.read_csv(bin3, sep="\t", header=None)
bin4_genes=pd.read_csv(bin41, sep="\t", header=None)

bin1_df = count_df[count_df["Gene_name"].isin(bin1_genes[0])]
fourfold_bin1 = bin1_df[bin1_df["Degeneracy"] == "4fold"]
zerofold_bin1 = bin1_df[bin1_df["Degeneracy"] == "0fold"]

bin2_df = count_df[count_df["Gene_name"].isin(bin2_genes[0])]
fourfold_bin2 = bin2_df[bin2_df["Degeneracy"] == "4fold"]
zerofold_bin2 = bin2_df[bin2_df["Degeneracy"] == "0fold"]

bin3_df = count_df[count_df["Gene_name"].isin(bin3_genes[0])]
fourfold_bin3 = bin3_df[bin3_df["Degeneracy"] == "4fold"]
zerofold_bin3 = bin3_df[bin3_df["Degeneracy"] == "0fold"]

bin4_df = count_df[count_df["Gene_name"].isin(bin3_genes[0])]
fourfold_bin4 = bin4_df[bin4_df["Degeneracy"] == "4fold"]
zerofold_bin4 = bin4_df[bin4_df["Degeneracy"] == "0fold"]

##### Create SFS
df_neut1= wrangle.unfolded_sfs(fourfold_bin1, nsamples)['count'].tolist()
df_neut2= wrangle.unfolded_sfs(fourfold_bin2, nsamples)['count'].tolist()
df_neut3= wrangle.unfolded_sfs(fourfold_bin3, nsamples)['count'].tolist()
df_neut4= wrangle.unfolded_sfs(fourfold_bin4, nsamples)['count'].tolist()

df_sel1 = wrangle.unfolded_sfs(zerofold_bin1, nsamples)['count'].tolist()
df_sel2 = wrangle.unfolded_sfs(zerofold_bin2, nsamples)['count'].tolist()
df_sel3 = wrangle.unfolded_sfs(zerofold_bin3, nsamples)['count'].tolist()
df_sel4 = wrangle.unfolded_sfs(zerofold_bin4, nsamples)['count'].tolist()

####  Calculate Pi

pi_4fold_bin1=stats.PiCalculate_uSFS(df_neut1, nsamples)
pi_4fold_bin2=stats.PiCalculate_uSFS(df_neut2, nsamples)
pi_4fold_bin3=stats.PiCalculate_uSFS(df_neut3, nsamples)
pi_4fold_bin4=stats.PiCalculate_uSFS(df_neut4, nsamples)

pi_0fold_bin1=stats.PiCalculate_uSFS(df_sel1, nsamples)
pi_0fold_bin2=stats.PiCalculate_uSFS(df_sel2, nsamples)
pi_0fold_bin3=stats.PiCalculate_uSFS(df_sel3, nsamples)
pi_0fold_bin4=stats.PiCalculate_uSFS(df_sel4, nsamples)

##### ##### ##### ##### ##### ##### ##### #####
##### #####  Bootstrap CI ### ##### ##### #####

##### Bin1

bs_ci_4fold_bin1=stats.bootstrap_pi(fourfold_bin1,  M=nsamples, ngenomes=nsamples)
bs_ci_0fold_bin1=stats.bootstrap_pi(zerofold_bin1,  M=nsamples, ngenomes=nsamples)

final_bin1 = pd.DataFrame({
    "Shetbin": ["bin1"] * 2,
    "Degeneracy": ["4fold", "0fold"],
    "Pi": [pi_4fold_bin1, pi_0fold_bin1],
    "BS_mean": [bs_ci_4fold_bin1[0], bs_ci_0fold_bin1[0]],
    "BS_LB": [bs_ci_4fold_bin1[1], bs_ci_0fold_bin1[1]],
    "BS_UB": [bs_ci_4fold_bin1[2], bs_ci_0fold_bin1[2]]
})

#### Bin 2

bs_ci_4fold_bin2=stats.bootstrap_pi(fourfold_bin2,  M=nsamples, ngenomes=nsamples)
bs_ci_0fold_bin2=stats.bootstrap_pi(zerofold_bin2,  M=nsamples, ngenomes=nsamples)

final_bin2 = pd.DataFrame({
    "Shetbin": ["bin2"] * 2,
    "Degeneracy": ["4fold", "0fold"],
    "Pi": [pi_4fold_bin2, pi_0fold_bin2],
    "BS_mean": [bs_ci_4fold_bin2[0], bs_ci_0fold_bin2[0]],
    "BS_LB": [bs_ci_4fold_bin2[1], bs_ci_0fold_bin2[1]],
    "BS_UB": [bs_ci_4fold_bin2[2], bs_ci_0fold_bin2[2]]
})

#### Bin3

bs_ci_4fold_bin3=stats.bootstrap_pi(fourfold_bin3,  M=nsamples, ngenomes=nsamples)
bs_ci_0fold_bin3=stats.bootstrap_pi(zerofold_bin3,  M=nsamples, ngenomes=nsamples)

final_bin3 = pd.DataFrame({
    "Shetbin": ["bin3"] * 2,
    "Degeneracy": ["4fold", "0fold"],
    "Pi": [pi_4fold_bin3, pi_0fold_bin3],
    "BS_mean": [bs_ci_4fold_bin3[0], bs_ci_0fold_bin3[0]],
    "BS_LB": [bs_ci_4fold_bin3[1], bs_ci_0fold_bin3[1]],
    "BS_UB": [bs_ci_4fold_bin3[2], bs_ci_0fold_bin3[2]]
})


#### Bin4

bs_ci_4fold_bin4=stats.bootstrap_pi(fourfold_bin4,  M=nsamples, ngenomes=nsamples)
bs_ci_0fold_bin4=stats.bootstrap_pi(zerofold_bin4,  M=nsamples, ngenomes=nsamples)

final_bin4 = pd.DataFrame({
    "Shetbin": ["bin4"] * 2,
    "Degeneracy": ["4fold", "0fold"],
    "Pi": [pi_4fold_bin4, pi_0fold_bin4],
    "BS_mean": [bs_ci_4fold_bin4[0], bs_ci_0fold_bin4[0]],
    "BS_LB": [bs_ci_4fold_bin4[1], bs_ci_0fold_bin4[1]],
    "BS_UB": [bs_ci_4fold_bin4[2], bs_ci_0fold_bin4[2]]
})

combined_df = pd.concat([final_bin1, final_bin2, final_bin3, final_bin4], ignore_index=True)
combined_df.to_csv(out_file1, sep="\t",index=False)

######### Bootstrap Pin/Ps #######

bs_ci_pinps_bin1 = stats.bootstrap_pnps(fourfold_bin1, zerofold_bin1, M=nsamples, ngenomes=nsamples)
bs_ci_pinps_bin2 = stats.bootstrap_pnps(fourfold_bin2, zerofold_bin2, M=nsamples, ngenomes=nsamples)
bs_ci_pinps_bin3 = stats.bootstrap_pnps(fourfold_bin3, zerofold_bin3, M=nsamples, ngenomes=nsamples)
bs_ci_pinps_bin4 = stats.bootstrap_pnps(fourfold_bin4, zerofold_bin4, M=nsamples, ngenomes=nsamples)

PnPs_all = pd.DataFrame({
    "Shetbin": ["bin1", "bin2", "bin3", "bin4"],
    "PnPs": [pi_0fold_bin1/pi_4fold_bin1, pi_0fold_bin2/pi_4fold_bin2,
             pi_0fold_bin3/pi_4fold_bin3,pi_0fold_bin4/pi_4fold_bin4],
    "BS_mean": [bs_ci_pinps_bin1[0], bs_ci_pinps_bin2[0],bs_ci_pinps_bin3[0],bs_ci_pinps_bin4[0]],
    "BS_LB": [bs_ci_pinps_bin1[1], bs_ci_pinps_bin2[1],bs_ci_pinps_bin3[1],bs_ci_pinps_bin4[1]],
    "BS_UB": [bs_ci_pinps_bin1[2], bs_ci_pinps_bin2[2],bs_ci_pinps_bin3[2],bs_ci_pinps_bin4[2]]
})

PnPs_all.to_csv(out_file2, sep="\t",index=False)
