import pandas as pd
import torch
import numpy as np
import sys
import pickle

outpath=sys.argv[1] # output directiory
prefix=sys.argv[2] # prefix of FastDTWF file
AFfile=sys.argv[3] # gene allele counts
out_prefix=sys.argv[4] #

with open(f"{AFfile}", "rb") as f:
     gene_to_counts = pickle.load(f)


s_grid = np.logspace(-8, 0, 100)
s_grid = np.append(0, s_grid)

likelihoods = []

for s_val in s_grid:
    filename = f"{outpath}/{prefix}_s_{s_val:.1e}"
    # Load the tensor
    freq_probs = torch.load(filename)
    likelihoods.append(freq_probs)

# Stack all tensors (shape) [101, nsamp+1] 101 Shet values, samp+1 frequency class
likelihood_tensor = torch.stack(likelihoods)


gene_loglikelihoods = {}
log_likelihood_tensor = torch.log(likelihood_tensor)  # shape [101, nsamp+1]
#log_likelihood_tensor = log_likelihood_tensor.masked_fill(torch.isinf(log_likelihood_tensor), -1e10)

for gene, counts in gene_to_counts.items():
    gene_logL = torch.zeros(log_likelihood_tensor.shape[0])
    for k in counts:
        gene_logL += log_likelihood_tensor[:, k] # All shet values for category k
    gene_loglikelihoods[gene] = gene_logL


torch.save(gene_loglikelihoods, f"{outpath}/{out_prefix}_gene_loglikelihoods")

with open(f"{outpath}/{out_prefix}_gene_loglikelihoods.pkl", "wb") as f:
    pickle.dump(gene_loglikelihoods, f)


#torch.load(f"{outpath}/{out_prefix}_gene_loglikelihoods")
