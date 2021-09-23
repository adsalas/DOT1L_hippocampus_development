# cellrank analysis for computing cell fate probabilities of the mature hippocampal cell lineages

## Import the required libraries


```python
import sys
import os
import matplotlib
import numpy as np
import pandas as pd
import scanpy as scp
import scvelo as scv
import cellrank as cr
```


```python
# Define some settings
scv.settings.verbosity = 3
scv.settings.set_figure_params('scvelo')
cr.settings.verbosity = 2
```


```python
# Load the sc object generated after running scvelo
sc = scp.read_h5ad("/path/to/the/saved/sc_scvelo_2")
```


```python
# Import some metadata contained in the original RaceID object and that will be required later
sc_Final_colors = pd.read_table("/path/to/the/cells_barcode_clusters_colors_tsne.tsv", sep= " ")
```


```python
# Create 2 new columns in the sc object for storing the tsne1 and tsne2 coordinates
sc.obs['tnse1'] = 0.0
sc.obs['tnse2'] = 0.0
```


```python
# Add the original t-SNE coordinates individually to the sc object in the .obs layer
for i in range(len(sc.obs)): 
    for j in range(len(sc_Final_colors)): 
        if sc.obs.cell_name[i] == sc_Final_colors.Exper_barcode[j]:
            sc.obs.tnse1[i] = sc_Final_colors.tsne1[j]
            sc.obs.tnse2[i] = sc_Final_colors.tsne2[j]
            break
```


```python
# Combine the 'tsne' coordinates and store them in the .obsm layer
A = sc.obs['tnse1']
B = sc.obs['tnse2']
d = {'tsne1': A, 'tsne2': B}
df = pd.DataFrame(d)
sc.obsm['tSNE2'] = np.array(df[['tsne1','tsne2']])
sc.obsm
```


```python
# Create a column to store the new cluster name and populate with 'unknown'
sc.obs['new_cluster'] = "unknown"
```


```python
# Add the corresponding cluster collapsed "new name" to each cell matching the original data set
for i in range(len(sc.obs)): 
    for j in range(len(sc_Final_colors)): 
        if sc.obs.cell_name[i] == sc_Final_colors.Exper_barcode[j]:
            sc.obs.new_cluster[i] = sc_Final_colors.new_cluster[j]
            break
```


```python
# Filter the cells to be kept for the cellrank computations
clusterNames=sc.obs.new_cluster.unique()
clusterNamesFTD = clusterNames[[0,1,2,4,5,6,7,10,11,12,13,14,16]]
clusterNamesFTD

sc1 = sc[(sc.obs['new_cluster'] == clusterNamesFTD[0]) | # Granule
              (sc.obs['new_cluster'] == clusterNamesFTD[1]) | # IPC2
              (sc.obs['new_cluster'] == clusterNamesFTD[2]) | # CA1pyr
              (sc.obs['new_cluster'] == clusterNamesFTD[3]) | # NSC1
              (sc.obs['new_cluster'] == clusterNamesFTD[4]) | # NSC4cyc
              (sc.obs['new_cluster'] == clusterNamesFTD[5]) | # CA3pyr
              (sc.obs['new_cluster'] == clusterNamesFTD[6]) | # Subic1
              (sc.obs['new_cluster'] == clusterNamesFTD[7]) | # NSC2cyc
              (sc.obs['new_cluster'] == clusterNamesFTD[8]) | # Subic2
              (sc.obs['new_cluster'] == clusterNamesFTD[9]) | # IPC1           
              (sc.obs['new_cluster'] == clusterNamesFTD[10]) | # NSC5
              (sc.obs['new_cluster'] == clusterNamesFTD[11]) | # CortHem
              (sc.obs['new_cluster'] == clusterNamesFTD[12])   # NSC3           
        ]

# Make a copy of the sc1 object 
sc_2 = sc1.copy()
```


```python
# Normalize per cell
scv.pp.normalize_per_cell(sc_2, enforce = False)
```


```python
# Extract the HVGs
scv.pp.filter_genes_dispersion(sc_2)
```


```python
# Perform the preprocessing steps

scv.pp.log1p(sc_2)
scv.pp.neighbors(sc_2, n_pcs=30, n_neighbors=30)
scv.pp.moments(sc_2, n_pcs = 30, n_neighbors =30, use_highly_variable= True)

# compute velocity and velocity graph
scv.tl.recover_dynamics(sc_2)
scv.tl.velocity(sc_2, mode = 'dynamical', enforce= True)
scv.tl.velocity_graph(sc_2)

# Maximal cosine correlation for each cell (used by scVelo to calculate self-transition probabilities)
sc_2.obs['max_cosine_corr'] = sc_2.uns['velocity_graph'].max(1).A.flatten()

scv.tl.recover_latent_time(sc_2)
scv.tl.velocity_confidence(sc_2)

scp.tl.pca(sc_2)
```


```python
# Define a variable named "cluster_colors_Final" to use for plotting
colors = np.array(sc_2.obs["colorFinal"])
colors
sc_2.uns["cluster_colors_Final"] = colors
```


```python
# Plot the scvelo estimated velocity vectors on the t-SNE representation obtained from RaceID analysis
scv.pl.velocity_embedding_stream(sc_2, basis='tSNE2', color= sc_2.uns["cluster_colors_Final"], legend_fontsize=12, title='', smooth=.8, min_mass=4)
```

## Computing Macrostates


```python
k = cr.tl.transition_matrix(sc_2, 
    weight_connectivities=0.2, 
    scheme= 'correlation',
    show_progress_bar=False                 
)
g = cr.tl.estimators.GPCCA(k)
```


```python
g.compute_schur(n_components=5)
```


```python
g.compute_macrostates(
    n_states=2, 
    cluster_key="new_cluster")
```


```python
# Define the terminal states (more mature hippocampal cell lineages)
g.set_terminal_states({"CA3pyr": sc_2[sc_2.obs["new_cluster"] == "CA3pyr"].obs_names, 
                       "CA1pyr": sc_2[sc_2.obs["new_cluster"] == "CA1pyr"].obs_names,
                      "Subic1": sc_2[sc_2.obs["new_cluster"] == "Subic1"].obs_names,
                      "Subic2": sc_2[sc_2.obs["new_cluster"] == "Subic2"].obs_names,
                      "Granule": sc_2[sc_2.obs["new_cluster"] == "Granule"].obs_names})
```


```python
# Compute lineage probabilities towards terminal states
cr.tl.lineages(sc_2, backward=False)
```

## Plotting the cell fate probabilities for the pre-defined terminal states


```python
cr.pl.lineages(sc_2, same_plot=False, basis = 'tSNE2')
```
