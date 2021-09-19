################## RNA velocity analysis with scvelo PART2: perform velocity estimations ################################

# written by Michael Stadler and Charlotte Soneson, FMI, Basel, Switzerland (2020)


import sys
import os
import matplotlib
import numpy as np
import pandas as pd
import scanpy as sc
import anndata2ri
import anndata
import scvelo as scv
from rpy2.robjects import r

# get input arguments from command line
n = len(sys.argv)
if (n == 3):
	scefile = sys.argv[1]
	adfileout = sys.argv[2]
else: 
	raise ValueError("must have two arguments (scefile, adfileout)")

print("input sce: ", scefile, flush = True)
print("output adfile: ", adfileout)

matplotlib.use('AGG')
sc.settings.verbosity = 3
scv.settings.verbosity = 3  

anndata2ri.activate()

# read sce and convert it to annData object
adata = r(f'sce <- readRDS("{scefile}")')
scv.utils.show_proportions(adata)

# preprocess the data
adata2 = scv.pp.filter_genes(adata, min_shared_counts = 30, copy = True)
scv.pp.normalize_per_cell(adata2, enforce = True)
scv.pp.filter_genes_dispersion(adata2, n_top_genes = 2000)
scv.pp.log1p(adata2)

## Convert unspliced matrix to dense format
adata2.layers['spliced'] = adata2.layers['spliced'].todense()
adata2.layers['unspliced'] = adata2.layers['unspliced'].todense()
scv.pp.moments(adata2, n_pcs = 30, n_neighbors = 30)

# compute velocity and velocity graph
scv.tl.recover_dynamics(adata2)
scv.tl.velocity(adata2, mode = 'dynamical')
scv.tl.velocity_graph(adata2)
## Maximal cosine correlation for each cell (used by scVelo to calculate self-transition probabilities)
adata2.obs['max_cosine_corr'] = adata2.uns['velocity_graph'].max(1).A.flatten()

scv.tl.recover_latent_time(adata2)
scv.tl.velocity_confidence(adata2)

# session info
scv.logging.print_version()
sc.logging.print_versions()

# save 
adata2.write(filename = adfileout)