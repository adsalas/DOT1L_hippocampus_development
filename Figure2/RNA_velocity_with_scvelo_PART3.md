# scvelo analysis for the hippocampal cell lineages


```python
# Import the required libraries
import sys
import os
import matplotlib
import numpy as np
import pandas as pd
import scanpy as scp
import scvelo as scv
import re
```


```python
# Define some settings
scv.settings.figdir = '/path/to/your/saving/directory/'
scv.settings.plot_prefix = 'scvelo_'
scv.settings.set_figure_params('scvelo')  # for beautified visualization
```


```python
# Read and import the object (scvelo object generated from Alevin's output pre-processed as a SingleCellExperiment)
sc = scp.read_h5ad("/path/to/the/saved/sce-scvelo.h5ad")
```


```python
# Import the metadata file
sc_metadata = pd.read_csv("/path/to/the/saved/metadata_sc_subcluster.csv")
```


```python
# Import other information contained in the original RaceID object that will be required later
sc_Final_colors = pd.read_table("/path/to/the/saved/cells_barcode_clusters_colors_tsne.tsv", sep= " ")
```


```python
# Create a new column to store the colors and populate with grey color code "#808080"
sc.obs['colorFinal'] = "#808080"
```


```python
# Changing the 'experiment' column information to upper case so that it matches the information of my original object
sc.obs.experiment = [x.upper() for x in sc.obs.experiment]
```


```python
# Create a cloumn for the 'samples' 
sc.obs['samples']= 'unknown'

for i in range(len(sc.obs)): 
    sc.obs.samples[i] = sc.obs.experiment[i][0:4]
```


```python
# Create a new column to store the cell name (experiment + barcode) and populate with 'unknown'
sc.obs['cell_name'] = "unknown"

# Populate the cell_name with the concatenation of the 'experiment' and the 'barcode'
for i in range(len(sc.obs)): 
    sc.obs.cell_name[i] = sc.obs.experiment[i] + '_' + sc.obs.barcode[i]
```


```python
# Create a new column to store the corrected cluster annotation
sc.obs['cluster_correc'] = "unknown"
```


```python
# Add the corresponding cluster annotations to each cell matching the original data set
for i in range(len(sc.obs)): 
    for j in range(len(sc_metadata)): 
        if sc.obs.cell_name[i] == sc_metadata.Exper_barcode[j]:
            sc.obs.cluster_correc[i] = sc_metadata.new_cluster[j]
        
```


```python
# Add a column for cluster names without sub-cluster information
sc.obs['cluster_collapsed'] = "unknown"
```


```python
# Populate the column
for i in range(len(sc.obs)): 
    # Correct the information for the CA clusters
    if sc.obs.cluster_correc[i][0:2] == "CA": 
        sc.obs.cluster_collapsed[i] = sc.obs.cluster_correc[i][0:6]
             
    else:
        sc.obs.cluster_collapsed[i] = sc.obs.cluster_correc[i]    
        
    # Correct a typo in the originally used label from 'NCS2_8' to 'NSC2_8'
    if sc.obs.cluster_correc[i] == 'NCS2_8':
        sc.obs.cluster_collapsed[i] = 'NSC2_8'
```


```python
# Add the corresponding cluster annotations to each cell matching the original data set
for i in range(len(sc.obs)): 
    for j in range(len(sc_Final_colors)): 
        if sc.obs.cell_name[i] == sc_Final_colors.Exper_barcode[j]:
            sc.obs.colorFinal[i] = sc_Final_colors.colors[j]
            break
        
```


```python
colors = np.array(sc.obs["colorFinal"])
```


```python
sc.uns["cluster_colors_Final"] = colors
```


```python
# Plot velocity embedding as stream
scv.pl.velocity_embedding_stream(sc, color= sc.uns["cluster_colors_Final"], figsize=(20,20))
```


```python
# Plot velocity embedding as grid
scv.pl.velocity_embedding_grid(sc, color= sc.uns["cluster_colors_Final"], figsize=(20,20))
```


```python
# Make a copy of the sc object containing just CA1 subclusters
sc1 = sc[(sc.obs['cluster_correc'].str.contains('CA1pyr1',regex = True)) | (sc.obs['cluster_correc'].str.contains('CA1pyr2',regex = True)), :].copy()
```


```python
# Checking the velocity vectors for the cells in the CA1 subclusters 1 and 2
scv.pl.velocity_embedding_stream(sc1, color= 'cluster_correc', figsize=(10,10), basis='umap')
```


```python
# Make a copy of the sc object containing just CA3 subclusters
sc2 = sc[(sc.obs['cluster_correc'].str.contains('CA3pyr1',regex = True)) | (sc.obs['cluster_correc'].str.contains('CA3pyr2',regex = True)), :].copy()
```


```python
# Checking the velocity vectors for the cells in the CA3 subclusters 1 and 2
scv.pl.velocity_embedding_stream(sc2, color= 'cluster_correc', figsize=(10,10), basis='umap')
```
