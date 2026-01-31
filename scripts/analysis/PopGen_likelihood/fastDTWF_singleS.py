import fastDTWF
import torch
import sys


idx=sys.argv[1]
outpath=sys.argv[2]
s_grid_file=sys.argv[3] # file with all 101 selection coefficients
thr=sys.argv[4]

torch.set_num_threads(int(thr))


with open(s_grid_file) as f:
    s_val = float(f.readlines()[int(idx)].strip())

### Drosophila demographic and genomic parameters (Johri et al. 2020, Ragsdale and Gutenkunst 2017)
## change for other species
Nanc=1225393
Ncur=1357760
Tgrowth=500000
mu = torch.tensor(3e-9, dtype=torch.float64)
s_het = torch.tensor(s_val, dtype=torch.float64)

###

freq_probs = fastDTWF.get_likelihood(
    pop_size_list=[Nanc, Ncur],    # population sizes for two epoch (haploids)
    switch_points=[Tgrowth,0],        # time for popsize change
    sample_size=69,        # sample from pop (haploids)
    s_het=s_het,
    mu_0_to_1=mu,
    mu_1_to_0=mu,
    dtwf_tv_sd=0.1,
    dtwf_row_eps=1e-8,
    sampling_tv_sd=0.05,
    sampling_row_eps=1e-8,
    no_fix=True,             # Whether to condition on non-fixation
    sfs=False,                # Whether to use an infinite sites model
    injection_rate=0.,        # Mutation rate under infinite sites model
    use_condensed=True,
    refresh_gens=500
)


#np.savetxt(f"{outpath}/likelihood_s_{s_val:.1e}", freq_probs.numpy())
torch.save(freq_probs, f"{outpath}/likelihood_s_{s_val:.1e}")

