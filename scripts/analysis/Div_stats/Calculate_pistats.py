import pandas as pd
import sys
from pathlib import Path
analysis_dir = Path(__file__).resolve().parents[1]
sys.path.append(str(analysis_dir))
import stats


count_file=sys.argv[1] # Population genetics (allele count) data
out_file=sys.argv[2]
nsamp=sys.argv[3] # Sample size Mouse=20, Yeast=16, Dmel=69

count_df=pd.read_csv(count_file, delimiter="\t")

fourfold = count_df[count_df["Degeneracy"] == "4fold"]
zerofold = count_df[count_df["Degeneracy"] == "0fold"]

fourfold_spectra= stats.unfolded_sfs(fourfold, M=nsamp)
zerofold_spectra= stats.unfolded_sfs(zerofold, M=nsamp)

pi_4fold=stats.PiCalculate_uSFS(fourfold_spectra, nsamp)
pi_0fold=stats.PiCalculate_uSFS(zerofold_spectra, nsamp)

#######

bs_ci_4fold=stats.bootstrap_pi(fourfold,  M=nsamp, ngenomes=nsamp)
bs_ci_0fold=stats.bootstrap_pi(zerofold,  M=nsamp, ngenomes=nsamp)

final = pd.DataFrame({
    "Degeneracy": ["4fold", "0fold"],
    "Pi": [pi_4fold, pi_0fold],
    "BS_mean": [bs_ci_4fold[0], bs_ci_0fold[0]],
    "BS_LB": [bs_ci_4fold[1], bs_ci_0fold[1]],
    "BS_UB": [bs_ci_4fold[2], bs_ci_0fold[2]]
})

final.to_csv(out_file, sep="\t",index=False)


