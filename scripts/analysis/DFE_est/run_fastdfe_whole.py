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

def save_inference_parameters(inf_obj, output_prefix, output_suffix="full_params_mle"):
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

    # Save
    full_df.to_csv(f"{output_prefix}_{output_suffix}", sep="\t", index=False)


def save_discretized_dfe(inf_obj, output_prefix, output_suffix="full_discretized_dfe"):
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

    dfe_df.to_csv(f"{output_prefix}_{output_suffix}", sep="\t", index=False)


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


#########
count_file=sys.argv[1] # Population genetics (allele count) data
out= sys.argv[2] # out_directory/prefix
nsamp= sys.argv[3] ## Sample size Mouse=20, Yeast=16, Dmel=69

### Read data

count_df=pd.read_csv(count_file, delimiter="\t")

## Build SFS
fourfold = bin1_df[bin1_df["Degeneracy"] == "4fold"]
zerofold = bin1_df[bin1_df["Degeneracy"] == "0fold"]

df_neut= wrangle.unfolded_sfs(fourfold, nsamples)['count'].tolist()
df_sel= wrangle.unfolded_sfs(zerofold, nsamples)['count'].tolist()

sfs_neut = fd.Spectrum(df_neut)
sfs_sel = fd.Spectrum(df_sel)

## inference with full model (Gamma Exponential parametrization)

basic_inf= fd.BaseInference(sfs_neut=sfs_neut, sfs_sel=sfs_sel, n_runs=10, n=nsamp, do_bootstrap=True)
basic_inf.run()

### Save whole inference object
basic_inf.to_file(file=out+"_full")
# inferred_full=fd.BaseInference.from_file(file=out+"_full")

### Save plots ###
basic_inf.plot_discretized(file=out+"_full_discrete_DFE")
basic_inf.plot_sfs_comparison(file=out+"_full_SFScomp")

### Save parameters ###
save_inference_parameters(basic_inf, output_prefix= out, output_suffix="full_params_mle")
### Save discretized DFE ###
save_discretized_dfe(basic_inf, output_prefix= out, output_suffix="full_discretized_dfe")


######## Compare_nested_models ############

nested_models=basic_inf.compare_nested_models()
save_pval_matrix(nested_models, out_prefix=out)

###
base_inference_dict=nested_models[1]

full_no_anc_inference = base_inference_dict['full.no_anc']
full_no_anc_inference.to_file(file=out+"_full_no_anc")

dele_no_anc_inference = base_inference_dict['dele.no_anc']
dele_no_anc_inference.to_file(file=out+"_dele_no_anc")

dele_anc_inference = base_inference_dict['dele.anc']
dele_anc_inference.to_file(file=out+"_dele_anc")

### Plot of nested model comparison
basic_inf.plot_nested_models(file=out+"_nested_model_comps")
