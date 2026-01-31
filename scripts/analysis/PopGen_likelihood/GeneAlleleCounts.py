import pandas as pd
import torch
import numpy as np
import sys
import pickle


info_full= sys.argv[1] # Population genetics count data
outpath=sys.argv[2] # Output directory


df_full = pd.read_csv(info_full, sep="\t", low_memory=False)
df_full["Gene_name"] = df_full["Gene_name"].str.split("_").str[0]

# Only keep 0fold and 4fold sites (regardless of polymorphism)
deg_df = df_full[df_full["Degeneracy"].isin(["0fold", "4fold"])].copy()

zerofold_df = deg_df[deg_df["Degeneracy"] == "0fold"].copy()

gene_to_counts_full_0fold_all = (
    zerofold_df.groupby("Gene_name")["MAC"]
    .apply(list)
    .to_dict()
)

######### Full data- 0 fold sites, polymorphic

zerofold_poly_df = df_full[(df_full["State"] == "P") & (df_full["Degeneracy"] == "0fold")].copy()

gene_to_counts_full_0fold = (
    zerofold_poly_df.groupby("Gene_name")["MAC"]
    .apply(list)
    .to_dict()
)

# 1. Full data — All sites
with open(f"{out}/gene_to_counts_full_0fold_all.pkl", "wb") as f:
    pickle.dump(gene_to_counts_full_0fold_all, f)

# 2. Full data — poly sites (MAC)
with open(f"{out}/gene_to_counts_full_0fold.pkl", "wb") as f:
     pickle.dump(gene_to_counts_full_0fold, f)

