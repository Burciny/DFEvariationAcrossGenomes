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

################################################################################################################################################
#################################### GET ALL INFERENCES FOR ALL ########################################################################

# ===============================
# PARAMETER COUNTS (for AIC)
# ===============================

k_dict = {
    "full_anc": 5,
    "full_noanc": 4,
    "dele_anc": 3,
    "dele_noanc": 2
}

# ===============================
# MODEL RUNNER
# ===============================

def run_all_models(sfs_neut, sfs_sel):

    models = {}

    # K=5
    full_anc = fd.BaseInference(
        sfs_neut=sfs_neut,
        sfs_sel=sfs_sel,
        n_runs=10,
        n=nsamples,
        do_bootstrap=True,
        fixed_params=dict(all=dict(h=0.5)),
        bounds={'S_b': (0.00001, 10000)}
    )
    full_anc.run()
    models["full_anc"] = full_anc

    # K=4
    full_noanc = fd.BaseInference(
        sfs_neut=sfs_neut,
        sfs_sel=sfs_sel,
        n_runs=10,
        n=nsamples,
        fixed_params=dict(all=dict(eps=0, h=0.5)),
        do_bootstrap=True,
        bounds={'S_b': (0.00001, 10000)}
    )
    full_noanc.run()
    models["full_noanc"] = full_noanc

    # K=3
    dele_anc = fd.BaseInference(
        sfs_neut=sfs_neut,
        sfs_sel=sfs_sel,
        n_runs=10,
        n=nsamples,
        fixed_params=dict(all=dict(S_b=1, p_b=0, h=0.5)),
        do_bootstrap=True
    )
    dele_anc.run()
    models["dele_anc"] = dele_anc

    # K=2
    dele_noanc = fd.BaseInference(
        sfs_neut=sfs_neut,
        sfs_sel=sfs_sel,
        n_runs=10,
        n=nsamples,
        fixed_params=dict(all=dict(S_b=1, p_b=0, eps=0, h=0.5)),
        do_bootstrap=True
    )
    dele_noanc.run()
    models["dele_noanc"] = dele_noanc

    return models

# ===============================
# AIC CALCULATION
# ===============================

def calculate_aic_table(models, k_dict, bin_name):

    rows = []

    for name, inf in models.items():

        logL = inf.likelihood
        K = k_dict[name]
        AIC = 2*K - 2*logL

        rows.append({
            "bin": bin_name,
            "model": name,
            "logL": logL,
            "K": K,
            "AIC": AIC
        })

    df = pd.DataFrame(rows)

    df["delta_AIC"] = df["AIC"] - df["AIC"].min()
    df["weight_raw"] = np.exp(-0.5 * df["delta_AIC"])
    df["AIC_weight"] = df["weight_raw"] / df["weight_raw"].sum()

    return df.sort_values("AIC").reset_index(drop=True)


# ===============================
# SETTINGS
# ===============================

sys.argv[1]=nsamples
sys.argv[2]=outpath
sys.argv[3]=out_prefix # species name
sys.argv[4]= out_suffix # suffix for inference


count_file=sys.argv[5] # Population genetics (allele count) data

### Gene names in each bin
bin1=sys.argv[6]
bin2=sys.argv[7]
bin3=sys.argv[8]
bin4=sys.argv[9]

##### Read and prepare data
count_df=pd.read_csv(count_file, delimiter="\t")
fourfold = count_df[count_df["Degeneracy"] == "4fold"]
zerofold = count_df[count_df["Degeneracy"] == "0fold"]

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

df_neut_whole= wrangle.unfolded_sfs(fourfold, nsamples)['count'].tolist()
df_sel_whole = wrangle.unfolded_sfs(zerofold, nsamples)['count'].tolist()


sfs_dict = {
    "full": (df_neut_whole, df_sel_whole
    ),
    "shet_bin1": (df_neut1, df_sel1
    ),
    "shet_bin2": (df_neut2, df_sel2
    ),
    "shet_bin3": (df_neut3, df_sel3
    ),
    "shet_bin4": (df_neut4, df_sel4
    ),
}

# ===============================
# MAIN LOOP
# ===============================

all_bin_results = {}
all_aic_tables = []

for name, (df_neut, df_sel) in sfs_dict.items():

    print(f"\nRunning: {name}")

    if sum(df_neut) == 0 or sum(df_sel) == 0:
        print("Skipping zero SFS")
        continue

    sfs_neut = fd.Spectrum(df_neut)
    sfs_sel  = fd.Spectrum(df_sel)

    models = run_all_models(sfs_neut, sfs_sel)
    all_bin_results[name] = models

    aic_table = calculate_aic_table(models, k_dict, name)
    all_aic_tables.append(aic_table)


# ===============================
# SAVE
# ===============================

final_aic_table = pd.concat(all_aic_tables, ignore_index=True)
final_aic_table.to_csv(f"{outpath}/{out_prefix}_{out_suffix}_AICresults", sep="\t", index=False)

with open(f"{outpath}/{out_prefix}_{out_suffix}_DFE_inference_objects.pkl", "wb") as f:
    pickle.dump(all_bin_results, f)


################################################################################################################################################
#################################### SAVE BEST AND MODEL AVERAGE PARAMETERS ########################################################################

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


def get_model_averaged_parameters(models, aic_table, inf_name, selected_models):

    models_used = ",".join(selected_models)

    top_aic = aic_table[aic_table["model"].isin(selected_models)].copy()
    top_aic["AIC_weight"] = top_aic["AIC_weight"] / top_aic["AIC_weight"].sum()

    weight_dict = dict(zip(top_aic["model"], top_aic["AIC_weight"]))

    rows = []

    for param in PARAMS_OF_INTEREST:

        val_dict = {}
        var_dict = {}

        for model_name in weight_dict.keys():

            inf = models[model_name]

            # Get MLE value
            if param in inf.params_mle:
                val = inf.params_mle[param]
            elif param == "alpha":
                val = inf.alpha
            elif param == "theta":
                val = inf.theta
            else:
                val = 0.0  # param not relevant for this model (e.g. S_b in dele models)

            # Get bootstrap variance
            if param == "S_b" and model_name in ["dele_anc", "dele_noanc"]:
                var_model = 0.0
            elif param in ("theta",):
                var_model = 0.0  # theta not in bootstraps
            elif inf.bootstraps is not None and param in inf.bootstraps.columns:
                var_model = inf.bootstraps[param].var()
            else:
                var_model = 0.0

            val_dict[model_name] = val
            var_dict[model_name] = var_model

        avg_est = sum(weight_dict[m] * val_dict[m] for m in weight_dict)

        var_est = sum(
            weight_dict[m] * (var_dict[m] + (val_dict[m] - avg_est) ** 2)
            for m in weight_dict
        )

        se = np.sqrt(var_est)

        rows.append({
            "Parameter": param,
            "MLE": avg_est,
            "Lower_Bound": avg_est - 1.96 * se,
            "Upper_Bound": avg_est + 1.96 * se,
            "inference": inf_name,
            "models_used": models_used
        })

    return pd.DataFrame(rows)

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



# ===============================
# MAIN LOOP
# ===============================

##########################################################################################

aic = pd.read_csv(f"{outpath}/{out_prefix}_{out_suffix}_AICresults", delimiter="\t")

with open(f"{outpath}/{out_prefix}_{out_suffix}_DFE_inference_objects.pkl", "rb") as f:
    DFE_infs = pickle.load(f)


all_bestfit_results = []
all_dfe_results = []
all_modelavg_results = []

for b in aic["bin"].unique():

    print("Processing bin:", b)

    models = DFE_infs[b]
    aic_bin = aic[aic["bin"] == b].sort_values("AIC").copy()
    aic_bin["AIC_weight"] = aic_bin["AIC_weight"]  # already computed
    aic_bin["cum_weight"] = aic_bin["AIC_weight"].cumsum()

    # --- Best fit ---
    best_model_name = aic_bin["model"].iloc[0]
    best_model = models[best_model_name]
    all_bestfit_results.append(get_inference_parameters(best_model, b))
    all_dfe_results.append(get_discretized_dfe(best_model, b))

    # --- Model averaging ---
    selected_models = aic_bin.loc[aic_bin["cum_weight"] <= 0.95, "model"].tolist()
    if len(aic_bin.loc[aic_bin["cum_weight"] > 0.95]) > 0:
        selected_models.append(aic_bin.loc[aic_bin["cum_weight"] > 0.95].iloc[0]["model"])
    all_modelavg_results.append(get_model_averaged_parameters(models, aic_bin, b, selected_models))



pd.concat(all_bestfit_results, ignore_index=True).to_csv(
    f"{outpath}/{out_prefix}_{out_suffix}_bestfit_parameters", sep="\t", index=False
)

pd.concat(all_dfe_results, ignore_index=True).to_csv(
    f"{outpath}/{out_prefix}_{out_suffix}_bestfit_dfe", sep="\t", index=False
)

pd.concat(all_modelavg_results, ignore_index=True).to_csv(
    f"{outpath}/{out_prefix}_{out_suffix}_model_averaged_parameters", sep="\t", index=False
)
