import pandas as pd
import sys
from pathlib import Path
analysis_dir = Path(__file__).resolve().parents[1]
sys.path.append(str(analysis_dir))
import stats

count_file=sys.argv[1] # Population genetics (allele count) data
out_file=sys.argv[2]
status=sys.argv[3] # "unfolded" or "folded" SFS
nsamp=sys.argv[4] # Sample size Mouse=20, Yeast=16, Dmel=69

count_df=pd.read_csv(count_file, delimiter="\t")

# Unique genes
all_genes = count_df["Gene_name"].unique()

# Compute pi per gene
results = []
for gene in all_genes:
    gene_df = count_df[count_df["Gene_name"] == gene]

    # Separate 4fold and 0fold
    fourfold = gene_df[gene_df["Degeneracy"] == "4fold"]
    zerofold = gene_df[gene_df["Degeneracy"] == "0fold"]

    if (status=="unfolded"):
        # Build SFS
        sfs_4fold = stats.unfolded_sfs(fourfold, M=nsamp)
        sfs_0fold = stats.unfolded_sfs(zerofold, M=nsamp)

        # Calculate pi
        pi_4fold = stats.PiCalculate_uSFS(sfs_4fold, nsamp) if sum(sfs_4fold) > 0 else float("nan")
        pi_0fold = stats.PiCalculate_uSFS(sfs_0fold, nsamp) if sum(sfs_0fold) > 0 else float("nan")
    else:
        # Build SFS
        sfs_4fold = stats.folded_sfs(fourfold, M=nsamp)
        sfs_0fold = stats.folded_sfs(zerofold, M=nsamp)

        # Calculate pi
        pi_4fold = stats.PiCalculate_fSFS(sfs_4fold, nsamp) if sum(sfs_4fold) > 0 else float("nan")
        pi_0fold = stats.PiCalculate_fSFS(sfs_0fold, nsamp) if sum(sfs_0fold) > 0 else float("nan")

    # Store
    results.append({
        "Gene": gene,
        "pi_4fold": pi_4fold,
        "pi_0fold": pi_0fold
    })

# Convert to DataFrame
pi_df = pd.DataFrame(results)
pi_df.to_csv(out_file, sep="\t",index=False)

