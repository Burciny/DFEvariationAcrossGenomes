import fastdfe as fd
import pandas as pd
import numpy as np
import sys
analysis_dir = Path(__file__).resolve().parents[1]
sys.path.append(str(analysis_dir))
import stats
import wrangle

############## Activate this for alpha calculation with modified threshold ##############

# ---- custom alpha calculation ----
# def custom_get_alpha(self, model, params, s_min=1):
#     # original code copied IDENTICALLY except threshold changed
#     y = model._discretize(params, self.bins) * fd.discretization.H_fixed_regularized(self.s)
#     return np.sum(y[self.s > s_min]) / np.sum(y)
#
# # ---- monkey patch ----
# fd.discretization.Discretization.get_alpha = custom_get_alpha

############## ############## ############## ############## ##############

# ===============================
# FUNCTIONS
# ===============================

PARAMS_OF_INTEREST = ['S_d', 'b', 'p_b', 'S_b', 'eps', 'alpha']


def get_inference_parameters(inf_obj, inf_name):
    # MLE parameters (excluding 'h' and anything else unwanted)
    params_mle = {k: v for k, v in inf_obj.params_mle.items()
                  if k in PARAMS_OF_INTEREST}
    params_df = pd.DataFrame(list(params_mle.items()), columns=["Parameter", "MLE"])

    # Add alpha and theta
    params_df = pd.concat([params_df, pd.DataFrame({
        "Parameter": ["alpha", "theta"],
        "MLE": [inf_obj.alpha, inf_obj.theta]
    })], ignore_index=True)

    # CIs from bootstraps — only for columns in PARAMS_OF_INTEREST
    if inf_obj.bootstraps is not None:
        bs = inf_obj.bootstraps
        ci_rows = []
        for param in PARAMS_OF_INTEREST:
            if param in bs.columns:
                ci_rows.append({
                    "Parameter": param,
                    "Lower_Bound": bs[param].quantile(0.025),
                    "Upper_Bound": bs[param].quantile(0.975)
                })
            else:
                ci_rows.append({
                    "Parameter": param,
                    "Lower_Bound": np.nan,
                    "Upper_Bound": np.nan
                })
        cis_df = pd.DataFrame(ci_rows)
    else:
        cis_df = pd.DataFrame({
            "Parameter": PARAMS_OF_INTEREST,
            "Lower_Bound": [np.nan] * len(PARAMS_OF_INTEREST),
            "Upper_Bound": [np.nan] * len(PARAMS_OF_INTEREST)
        })

    # theta has no bootstrap CI
    cis_df = pd.concat([cis_df, pd.DataFrame({
        "Parameter": ["theta"],
        "Lower_Bound": [np.nan],
        "Upper_Bound": [np.nan]
    })], ignore_index=True)

    full_df = pd.merge(params_df, cis_df, on="Parameter", how="left")
    full_df["inference"] = inf_name

    return full_df



def get_discretized_dfe(inf_obj,inf_name):
    """
    Save discretized DFE values and their standard deviations from a FastDFE inference object.

    Parameters:
        inf_obj (BaseInference): The inference object.
        output_prefix (str): Path prefix for the output file.
        output_suffix (str): Suffix for the output file. for model differentiation
    """
    dfe_values, dfe_deviations = inf_obj.get_discretized(confidence_intervals=True)

    dfe_df = pd.DataFrame({
        "Bin": [f"bin_{i+1}" for i in range(len(dfe_values))],
        "Value": dfe_values,
        "Std_Lower": dfe_deviations[0],
        "Std_Upper": dfe_deviations[1]
    })
    dfe_df["inference"] = [inf_name] * len(dfe_df)

    return dfe_df

    #dfe_df.to_csv(f"{output_prefix}_{output_suffix}", sep="\t", index=False)


#########
########################################################################
######### Polarized Data unfolded inference ################################

count_file=sys.argv[1] # Population genetics (allele count) data

### Gene names in each bin
bin1=sys.argv[2]
bin2=sys.argv[3]
bin3=sys.argv[4]
bin4=sys.argv[5]

nsamples=int(sys.argv[6])
out= sys.argv[7] # output_directory/prefix

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

sfs_neut = fd.Spectra(dict(bin1=df_neut1, bin2=df_neut2, bin3=df_neut3, bin4=df_neut4))
sfs_sel = fd.Spectra(dict(bin1=df_sel1, bin2=df_sel2, bin3=df_sel3, bin4=df_sel4))

#### Run inference

# Do not share eps

basic_joint_vareps_inf= fd.JointInference(sfs_neut=sfs_neut, sfs_sel=sfs_sel, n_runs=10, n=nsamples, do_bootstrap=True,
                                   shared_params = [fd.SharedParams(types=["bin1", "bin2", "bin3", "bin4"],
                                   params=['S_d', 'b', 'p_b', 'S_b'])],
                                          bounds={'S_b': (0.00001, 10000)}, fixed_params=dict(all=dict(h=0.5)))

basic_joint_vareps_inf.run()

#### SAVE

basic_joint_vareps_inf.to_file(file=out+"_full")
jointinferred_vareps_full=fd.JointInference.from_file(file=out+"_full")


### Marginal and joint inferences
bin1_marginal_vareps=jointinferred_vareps_full.marginal_inferences["bin1"]
bin1_joint_vareps=jointinferred_vareps_full.joint_inferences["bin1"]

bin2_marginal_vareps=jointinferred_vareps_full.marginal_inferences["bin2"]
bin2_joint_vareps=jointinferred_vareps_full.joint_inferences["bin2"]

bin3_marginal_vareps=jointinferred_vareps_full.marginal_inferences["bin3"]
bin3_joint_vareps=jointinferred_vareps_full.joint_inferences["bin3"]

bin4_marginal_vareps=jointinferred_vareps_full.marginal_inferences["bin4"]
bin4_joint_vareps=jointinferred_vareps_full.joint_inferences["bin4"]


#### Save info

joint_distDFE=get_discretized_dfe(bin1_joint_vareps, "joint")
bin1_distDFE=get_discretized_dfe(bin1_marginal_vareps,"shet_bin1")
bin2_distDFE=get_discretized_dfe(bin2_marginal_vareps,"shet_bin2")
bin3_distDFE=get_discretized_dfe(bin3_marginal_vareps,"shet_bin3")
bin4_distDFE=get_discretized_dfe(bin4_marginal_vareps,"shet_bin4")

jointinferred_full_distDFE= pd.concat([joint_distDFE,bin1_distDFE,bin2_distDFE, bin3_distDFE, bin4_distDFE])
jointinferred_full_distDFE.to_csv(file=out+"_discretizedDFE", sep="\t", index=False)

joint_params=get_inference_parameters(bin1_joint_vareps, "joint")
bin1_params=get_inference_parameters(bin1_marginal_vareps, "shet_bin1")
bin2_params=get_inference_parameters(bin2_marginal_vareps, "shet_bin2")
bin3_params=get_inference_parameters(bin3_marginal_vareps, "shet_bin3")
bin4_params=get_inference_parameters(bin4_marginal_vareps, "shet_bin4")

jointinferred_full_params= pd.concat([joint_params,bin1_params,bin2_params, bin3_params, bin4_params])
jointinferred_full_params.to_csv(file=out+ "_paramsMLE", sep="\t", index=False)
