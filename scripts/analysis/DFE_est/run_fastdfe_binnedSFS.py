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

#### helper functions ##

def get_inference_parameters(inf_obj,inf_name):
    """
    Save MLE parameter estimates, alpha and theta, and their confidence intervals (if available)
    from a FastDFE BaseInference object.

    Parameters:
        inf_obj (BaseInference): The inference object.
        output_prefix (str): Path prefix for the output file.
        output_suffix (str): Suffix for the output file. for model differentiation
    """
    # MLE parameters
    params_mle = inf_obj.params_mle
    params_df = pd.DataFrame(list(params_mle.items()), columns=["Parameter", "MLE"])

    # Add alpha and theta
    extra_params_df = pd.DataFrame({
        "Parameter": ["alpha", "theta"],
        "MLE": [inf_obj.alpha, inf_obj.theta]
    })
    params_df = pd.concat([params_df, extra_params_df], ignore_index=True)

    # Confidence intervals
    cis_dict = inf_obj.get_cis_params_mle()
    cis_df = pd.DataFrame(list(cis_dict.items()), columns=["Parameter", "CI"])
    cis_df["Lower_Bound"] = cis_df["CI"].apply(lambda x: x[0])
    cis_df["Upper_Bound"] = cis_df["CI"].apply(lambda x: x[1])
    cis_df.drop("CI", axis=1, inplace=True)

    # Add placeholder CI for theta (if not in original CI dict)
    if "theta" not in cis_df["Parameter"].values:
        cis_df = pd.concat([cis_df, pd.DataFrame({
            "Parameter": ["theta"],
            "Lower_Bound": [np.nan],
            "Upper_Bound": [np.nan]
        })], ignore_index=True)

    # Merge everything
    full_df = pd.merge(params_df, cis_df, on="Parameter", how="left")
    full_df["inference"] = [inf_name] * len(full_df)
    return full_df

    # Save
#    full_df.to_csv(f"{output_prefix}_{output_suffix}", sep="\t", index=False)


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

def save_pval_matrix(nested_inf_obj, out_prefix, model_labels= ["full.no_anc", "full.anc", "dele.no_anc", "dele.anc"]):
    """
    Save a nested model comparison p-value matrix to a tab-separated file.

    Parameters:
    - pval_matrix: 2D array of p-values (dtype=object, may include None)
    - model_labels: list of strings used for row/column labels
    - out_prefix: string prefix for output file
    """
    # Create DataFrame from the matrix
    pval_matrix = nested_inf_obj[0]

    pval_df = pd.DataFrame(pval_matrix, index=model_labels, columns=model_labels)

    # Save as tab-separated file without .csv extension
    pval_df.to_csv(out_prefix + "_nested_model_pvals", sep="\t", index=True)


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
                                   params=['S_d', 'b', 'p_b', 'S_b'])],bounds={'S_b': (0.00001, 10000)})

basic_joint_inf.run()
basic_joint_vareps_inf.run()


#### SAVE

basic_joint_vareps_inf.to_file(file=out+"_full")
jointinferred_vareps_full=fd.JointInference.from_file(file=out+"_full")

#### Plots all shared

jointinferred_full.plot_inferred_parameters()
jointinferred_full.plot_discretized()

### Plots vareps

jointinferred_vareps_full.plot_inferred_parameters()
jointinferred_vareps_full.plot_discretized()

### LRT vareps
jointinferred_vareps_full.joint_inference_makes_sense()
jointinferred_vareps_full.perform_lrt_shared(do_bootstrap=True)

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


###################################################################################################################################################
##################### Nested models for each marginal inference ##################

#### Bin 1 ####
nested_models_bin1=bin1_marginal_vareps.compare_nested_models()
save_pval_matrix(nested_models_bin1, out_prefix=out+ "_bin1_nested_model_pvals")

base_inference_dict=nested_models_bin1[1]
full_no_anc_inference = base_inference_dict['full.no_anc']
full_no_anc_inference.to_file(file=out+"_full_no_anc")

dele_no_anc_inference = base_inference_dict['dele.no_anc']
dele_no_anc_inference.to_file(file=out+"_dele_no_anc")

dele_anc_inference = base_inference_dict['dele.anc']
dele_anc_inference.to_file(file=out+"_dele_anc")

bin1_marginal_vareps.plot_nested_models(file=out+"_bin1_nested_model_comps")

##### Bin 2 #####
nested_models_bin2=bin2_marginal_vareps.compare_nested_models()
save_pval_matrix(nested_models_bin2, out_prefix=out+ "_bin2_nested_model_pvals")

base_inference_dict=nested_models_bin2[1]
full_no_anc_inference = base_inference_dict['full.no_anc']
full_no_anc_inference.to_file(file=out+"_full_no_anc")

dele_no_anc_inference = base_inference_dict['dele.no_anc']
dele_no_anc_inference.to_file(file=out+"_dele_no_anc")

dele_anc_inference = base_inference_dict['dele.anc']
dele_anc_inference.to_file(file=out+"_dele_anc")

bin2_marginal_vareps.plot_nested_models(file=out+"_bin2_nested_model_comps")


##### Bin 3 #####
nested_models_bin3=bin3_marginal_vareps.compare_nested_models()
save_pval_matrix(nested_models_bin3, out_prefix=out+ "_bin3_nested_model_pvals")

base_inference_dict=nested_models_bin3[1]
full_no_anc_inference = base_inference_dict['full.no_anc']
full_no_anc_inference.to_file(file=out+"_full_no_anc")

dele_no_anc_inference = base_inference_dict['dele.no_anc']
dele_no_anc_inference.to_file(file=out+"_dele_no_anc")

dele_anc_inference = base_inference_dict['dele.anc']
dele_anc_inference.to_file(file=out+"_dele_anc")

bin3_marginal_vareps.plot_nested_models(file=out+"_bin3_nested_model_comps")


##### Bin 4 #####
nested_models_bin4=bin4_marginal_vareps.compare_nested_models()
out="/Users/burcinyildirim/Documents/Uppsala_ICM/Data_process/Mouse/CDS/FastDFE/BObinnedSFS/Nested_models/Mouse_marginal_BObin4"
save_pval_matrix(nested_models_bin4, out_prefix=out+ "_bin4_nested_model_pvals")

base_inference_dict=nested_models_bin4[1]
full_no_anc_inference = base_inference_dict['full.no_anc']
full_no_anc_inference.to_file(file=out+"_full_no_anc")

dele_no_anc_inference = base_inference_dict['dele.no_anc']
dele_no_anc_inference.to_file(file=out+"_dele_no_anc")

dele_anc_inference = base_inference_dict['dele.anc']
dele_anc_inference.to_file(file=out+"_dele_anc")

bin4_marginal_vareps.plot_nested_models(file=out+"_bin4_nested_model_comps")

