# import library
import pandas as pd
import numpy as np
import stlearn as st
import scanpy as sc
import warnings
import matplotlib.pyplot as plt
import seaborn as sns

warnings.filterwarnings("ignore")


### During
adata = st.ReadXenium(feature_cell_matrix_file="cell_feature_matrix.h5",
                     cell_summary_file="cells.csv.gz",
                     library_id="ESCC_During",
                     image_path="21_ESCC_CaseX1_During.png",
                     scale=1,
                     spot_diameter_fullres=15 # Recommend
                     )
                     
# Filter genes and cells with at least 10 counts
sc.pp.filter_genes(adata, min_counts=10)
sc.pp.filter_cells(adata,min_counts=10)

# Store the raw data for using PSTS
adata.raw = adata

# Normalization data
sc.pp.normalize_total(adata)

# Squareroot normalize transcript counts. We need to deal with sparse matrix of .X
from scipy.sparse import csr_array
adata.X = np.sqrt(adata.X.toarray()) + np.sqrt(adata.X.toarray() + 1)

# Run PCA
st.em.run_pca(adata,n_comps=50,random_state=0)

# Compute neighborhood graph of cells using the PCA representation
st.pp.neighbors(adata,n_neighbors=25,use_rep='X_pca',random_state=0)

st.tl.clustering.louvain(adata,random_state=0)

st.pl.cluster_plot(adata, use_label="louvain", image_alpha=0, size=4, figsize=(5, 5))

### Calculating the number of grid spots we will generate
n_ = 125
print(f'{n_} by {n_} has this many spots:\n', n_*n_)

### Gridding.
grid = st.tl.cci.grid(adata, n_row=n_, n_col=n_, use_label = 'louvain')
print( grid.shape ) # Slightly less than the above calculation, since we filter out spots with 0 cells.

import matplotlib.pyplot as plt
fig, axes = plt.subplots(ncols=2, figsize=(10,5))
st.pl.cluster_plot(grid, use_label='louvain', size=10, ax=axes[0], show_plot=False, show_image=False)
st.pl.cluster_plot(adata, use_label='louvain', ax=axes[1], show_plot=False, show_image=False)
axes[0].set_title(f'Grid louvain dominant spots')
axes[1].set_title(f'Cell louvain labels')
plt.show()

groups = list(grid.obs['louvain'].cat.categories)
for group in groups[0:2]:
    fig, axes = plt.subplots(ncols=3, figsize=(16,5))
    group_props = grid.uns['louvain'][group].values
    grid.obs['group'] = group_props
    st.pl.feat_plot(grid, feature='group', ax=axes[0], show_plot=False, vmax=1, show_color_bar=False, show_image=False)
    st.pl.cluster_plot(grid, use_label='louvain', list_clusters=[group], ax=axes[1], show_plot=False, show_image=False)
    st.pl.cluster_plot(adata, use_label='louvain', list_clusters=[group], ax=axes[2], show_plot=False, show_image=False)
    axes[0].set_title(f'Grid {group} proportions (max = 1)')
    axes[1].set_title(f'Grid {group} max spots')
    axes[2].set_title(f'Individual cell {group}')
    plt.show()

fig, axes = plt.subplots(ncols=2, figsize=(12,6))
st.pl.gene_plot(grid, gene_symbols='IFNG', ax=axes[0],show_color_bar=False,show_plot=False, show_image=False)
st.pl.gene_plot(adata, gene_symbols='IFNG', ax=axes[1], show_color_bar=False, show_plot=False, show_image=False)
axes[0].set_title(f'Grid IFNG expression')
axes[1].set_title(f'Cell IFNG expression')
plt.show()

fig, axes = plt.subplots(ncols=2, figsize=(12,6))
st.pl.gene_plot(grid, gene_symbols='IFNGR2', ax=axes[0],show_color_bar=False,show_plot=False, show_image=False)
st.pl.gene_plot(adata, gene_symbols='IFNGR2', ax=axes[1], show_color_bar=False, show_plot=False, show_image=False)
axes[0].set_title(f'Grid IFNGR2 expression')
axes[1].set_title(f'Cell IFNGR2 expression')
plt.show()

# Loading the LR databases available within stlearn (from NATMI)
lrs = st.tl.cci.load_lrs(['connectomeDB2020_lit'], species='human')
print(len(lrs))

# Running the analysis #
st.tl.cci.run(grid, lrs,
                  min_spots = 20, #Filter out any LR pairs with no scores for less than min_spots
                  distance=None, # None defaults to spot+immediate neighbours; distance=0 for within-spot mode
                  n_pairs=1000, # Number of random pairs to generate; low as example, recommend ~10,000
                  n_cpus=None, # Number of CPUs for parallel. If None, detects & use all available.
                  )
                  
### Can adjust significance thresholds.
st.tl.cci.adj_pvals(grid, correct_axis='spot',
                   pval_adj_cutoff=0.05, adj_method='fdr_bh')

best_lr = grid.uns['lr_summary'].index.values[0] # Just choosing one of the top from lr_summary

# Showing the rankings of the LR from a global and local perspective.
# Ranking based on number of significant hotspots.
plt.figure()
st.pl.lr_summary(grid, n_top=500)
st.pl.lr_summary(grid, n_top=50, figsize=(10,3))
plt.savefig("lr_rank_plot.png", dpi=500)

### Can adjust significance thresholds.
st.tl.cci.adj_pvals(grid, correct_axis='spot',
                   pval_adj_cutoff=0.05, adj_method='fdr_bh')

best_lr = grid.uns['lr_summary'].index.values[0] # Just choosing one of the top from lr_summary

stats = ['lr_scores', '-log10(p_adjs)', 'lr_sig_scores']
fig, axes = plt.subplots(ncols=len(stats), figsize=(18,6))
for i, stat in enumerate(stats):
    st.pl.lr_result_plot(grid, use_result=stat, use_lr=best_lr, show_color_bar=False, ax=axes[i], show_image=False)
    axes[i].set_title(f'{best_lr} {stat}')
    
stats = ['lr_scores', '-log10(p_adjs)', 'lr_sig_scores']
fig, axes = plt.subplots(ncols=len(stats), figsize=(18,6))
for i, stat in enumerate(stats):
    st.pl.lr_result_plot(grid, use_result=stat, use_lr="IFNG_IFNGR1", show_color_bar=False, ax=axes[i], show_image=False)
    axes[i].set_title(f'{"IFNG_IFNGR1"} {stat}')
    
stats = ['lr_scores', '-log10(p_adjs)', 'lr_sig_scores']
fig, axes = plt.subplots(ncols=len(stats), figsize=(18,6))
for i, stat in enumerate(stats):
    st.pl.lr_result_plot(grid, use_result=stat, use_lr="IFNG_IFNGR2", show_color_bar=False, ax=axes[i], show_image=False)
    axes[i].set_title(f'{"IFNG_IFNGR2"} {stat}')
    
stats = ['lr_scores', '-log10(p_adjs)', 'lr_sig_scores']
fig, axes = plt.subplots(ncols=len(stats), figsize=(18,6))
for i, stat in enumerate(stats):
    st.pl.lr_result_plot(grid, use_result=stat, use_lr="CD274_CD80", show_color_bar=False, ax=axes[i], show_image=False)
    axes[i].set_title(f'{"CD274_CD80"} {stat}')
    
stats = ['lr_scores', '-log10(p_adjs)', 'lr_sig_scores']
fig, axes = plt.subplots(ncols=len(stats), figsize=(18,6))
for i, stat in enumerate(stats):
    st.pl.lr_result_plot(grid, use_result=stat, use_lr="TGFB1_TGFBR2", show_color_bar=False, ax=axes[i], show_image=False)
    axes[i].set_title(f'{"TGFB1_TGFBR2"} {stat}')
    
### Pre
adata = st.ReadXenium(feature_cell_matrix_file="cell_feature_matrix.h5",
                     cell_summary_file="cells.csv.gz",
                     library_id="ESCC_Pre",
                     image_path="20_ESCC_CaseX1_Pre.png",
                     scale=1,
                     spot_diameter_fullres=15 # Recommend
                     )
                     
# Filter genes and cells with at least 10 counts
sc.pp.filter_genes(adata, min_counts=10)
sc.pp.filter_cells(adata,min_counts=10)

# Store the raw data for using PSTS
adata.raw = adata

# Normalization data
sc.pp.normalize_total(adata)

# Squareroot normalize transcript counts. We need to deal with sparse matrix of .X
from scipy.sparse import csr_array
adata.X = np.sqrt(adata.X.toarray()) + np.sqrt(adata.X.toarray() + 1)

# Run PCA
st.em.run_pca(adata,n_comps=50,random_state=0)

# Compute neighborhood graph of cells using the PCA representation
st.pp.neighbors(adata,n_neighbors=25,use_rep='X_pca',random_state=0)

st.tl.clustering.louvain(adata,random_state=0)

st.pl.cluster_plot(adata, use_label="louvain", image_alpha=0, size=4, figsize=(5, 5))

### Calculating the number of grid spots we will generate
n_ = 125
print(f'{n_} by {n_} has this many spots:\n', n_*n_)

### Gridding.
grid = st.tl.cci.grid(adata, n_row=n_, n_col=n_, use_label = 'louvain')
print( grid.shape ) # Slightly less than the above calculation, since we filter out spots with 0 cells.

import matplotlib.pyplot as plt
fig, axes = plt.subplots(ncols=2, figsize=(10,5))
st.pl.cluster_plot(grid, use_label='louvain', size=10, ax=axes[0], show_plot=False, show_image=False)
st.pl.cluster_plot(adata, use_label='louvain', ax=axes[1], show_plot=False, show_image=False)
axes[0].set_title(f'Grid louvain dominant spots')
axes[1].set_title(f'Cell louvain labels')
plt.show()

groups = list(grid.obs['louvain'].cat.categories)
for group in groups[0:2]:
    fig, axes = plt.subplots(ncols=3, figsize=(16,5))
    group_props = grid.uns['louvain'][group].values
    grid.obs['group'] = group_props
    st.pl.feat_plot(grid, feature='group', ax=axes[0], show_plot=False, vmax=1, show_color_bar=False, show_image=False)
    st.pl.cluster_plot(grid, use_label='louvain', list_clusters=[group], ax=axes[1], show_plot=False, show_image=False)
    st.pl.cluster_plot(adata, use_label='louvain', list_clusters=[group], ax=axes[2], show_plot=False, show_image=False)
    axes[0].set_title(f'Grid {group} proportions (max = 1)')
    axes[1].set_title(f'Grid {group} max spots')
    axes[2].set_title(f'Individual cell {group}')
    plt.show()

fig, axes = plt.subplots(ncols=2, figsize=(12,6))
st.pl.gene_plot(grid, gene_symbols='IFNG', ax=axes[0],show_color_bar=False,show_plot=False, show_image=False)
st.pl.gene_plot(adata, gene_symbols='IFNG', ax=axes[1], show_color_bar=False, show_plot=False, show_image=False)
axes[0].set_title(f'Grid IFNG expression')
axes[1].set_title(f'Cell IFNG expression')
plt.show()

fig, axes = plt.subplots(ncols=2, figsize=(12,6))
st.pl.gene_plot(grid, gene_symbols='IFNGR2', ax=axes[0],show_color_bar=False,show_plot=False, show_image=False)
st.pl.gene_plot(adata, gene_symbols='IFNGR2', ax=axes[1], show_color_bar=False, show_plot=False, show_image=False)
axes[0].set_title(f'Grid IFNGR2 expression')
axes[1].set_title(f'Cell IFNGR2 expression')
plt.show()

# Loading the LR databases available within stlearn (from NATMI)
lrs = st.tl.cci.load_lrs(['connectomeDB2020_lit'], species='human')
print(len(lrs))

# Running the analysis #
st.tl.cci.run(grid, lrs,
                  min_spots = 20, #Filter out any LR pairs with no scores for less than min_spots
                  distance=None, # None defaults to spot+immediate neighbours; distance=0 for within-spot mode
                  n_pairs=1000, # Number of random pairs to generate; low as example, recommend ~10,000
                  n_cpus=None, # Number of CPUs for parallel. If None, detects & use all available.
                  )
                  
lr_info = grid.uns['lr_summary'] # A dataframe detailing the LR pairs ranked by number of significant spots.
print(lr_info.shape)
print(lr_info.head(50))

# Showing the rankings of the LR from a global and local perspective.
# Ranking based on number of significant hotspots.
st.pl.lr_summary(grid, n_top=500)
st.pl.lr_summary(grid, n_top=50, figsize=(10,3))### Can adjust significance thresholds.
st.tl.cci.adj_pvals(grid, correct_axis='spot',
                   pval_adj_cutoff=0.05, adj_method='fdr_bh')

best_lr = grid.uns['lr_summary'].index.values[0] # Just choosing one of the top from lr_summary

stats = ['lr_scores', '-log10(p_adjs)', 'lr_sig_scores']
fig, axes = plt.subplots(ncols=len(stats), figsize=(18,6))
for i, stat in enumerate(stats):
    st.pl.lr_result_plot(grid, use_result=stat, use_lr=best_lr, show_color_bar=False, ax=axes[i], show_image=False)
    axes[i].set_title(f'{best_lr} {stat}')

stats = ['lr_scores', '-log10(p_adjs)', 'lr_sig_scores']
fig, axes = plt.subplots(ncols=len(stats), figsize=(18,6))
for i, stat in enumerate(stats):
    st.pl.lr_result_plot(grid, use_result=stat, use_lr="IFNG_IFNGR1", show_color_bar=False, ax=axes[i], show_image=False)
    axes[i].set_title(f'{"IFNG_IFNGR1"} {stat}')

stats = ['lr_scores', '-log10(p_adjs)', 'lr_sig_scores']
fig, axes = plt.subplots(ncols=len(stats), figsize=(18,6))
for i, stat in enumerate(stats):
    st.pl.lr_result_plot(grid, use_result=stat, use_lr="IFNG_IFNGR1", show_color_bar=False, ax=axes[i], show_image=False)
    axes[i].set_title(f'{"IFNG_IFNGR2"} {stat}')

stats = ['lr_scores', '-log10(p_adjs)', 'lr_sig_scores']
fig, axes = plt.subplots(ncols=len(stats), figsize=(18,6))
for i, stat in enumerate(stats):
    st.pl.lr_result_plot(grid, use_result=stat, use_lr="CD274_CD80", show_color_bar=False, ax=axes[i], show_image=False)
    axes[i].set_title(f'{"CD274_CD80"} {stat}')

stats = ['lr_scores', '-log10(p_adjs)', 'lr_sig_scores']
fig, axes = plt.subplots(ncols=len(stats), figsize=(18,6))
for i, stat in enumerate(stats):
    st.pl.lr_result_plot(grid, use_result=stat, use_lr="TGFB1_TGFBR2", show_color_bar=False, ax=axes[i], show_image=False)
    axes[i].set_title(f'{"TGFB1_TGFBR2"} {stat}')
    
### After
adata = st.ReadXenium(feature_cell_matrix_file="cell_feature_matrix.h5",
                     cell_summary_file="cells.csv.gz",
                     library_id="ESCC_After",
                     image_path="22_ESCC_CaseX1_After.png",
                     scale=1,
                     spot_diameter_fullres=15 # Recommend
                     )
                     
# Filter genes and cells with at least 10 counts
sc.pp.filter_genes(adata, min_counts=10)
sc.pp.filter_cells(adata,min_counts=10)

# Store the raw data for using PSTS
adata.raw = adata

# Normalization data
sc.pp.normalize_total(adata)

# Squareroot normalize transcript counts. We need to deal with sparse matrix of .X
from scipy.sparse import csr_array
adata.X = np.sqrt(adata.X.toarray()) + np.sqrt(adata.X.toarray() + 1)

# Run PCA
st.em.run_pca(adata,n_comps=50,random_state=0)

# Compute neighborhood graph of cells using the PCA representation
st.pp.neighbors(adata,n_neighbors=25,use_rep='X_pca',random_state=0)

st.tl.clustering.louvain(adata,random_state=0)

st.pl.cluster_plot(adata, use_label="louvain", image_alpha=0, size=4, figsize=(5, 5))

### Calculating the number of grid spots we will generate
n_ = 125
print(f'{n_} by {n_} has this many spots:\n', n_*n_)

### Gridding.
grid = st.tl.cci.grid(adata, n_row=n_, n_col=n_, use_label = 'louvain')
print( grid.shape ) # Slightly less than the above calculation, since we filter out spots with 0 cells.

import matplotlib.pyplot as plt
fig, axes = plt.subplots(ncols=2, figsize=(10,5))
st.pl.cluster_plot(grid, use_label='louvain', size=10, ax=axes[0], show_plot=False, show_image=False)
st.pl.cluster_plot(adata, use_label='louvain', ax=axes[1], show_plot=False, show_image=False)
axes[0].set_title(f'Grid louvain dominant spots')
axes[1].set_title(f'Cell louvain labels')
plt.show()

groups = list(grid.obs['louvain'].cat.categories)
for group in groups[0:2]:
    fig, axes = plt.subplots(ncols=3, figsize=(16,5))
    group_props = grid.uns['louvain'][group].values
    grid.obs['group'] = group_props
    st.pl.feat_plot(grid, feature='group', ax=axes[0], show_plot=False, vmax=1, show_color_bar=False, show_image=False)
    st.pl.cluster_plot(grid, use_label='louvain', list_clusters=[group], ax=axes[1], show_plot=False, show_image=False)
    st.pl.cluster_plot(adata, use_label='louvain', list_clusters=[group], ax=axes[2], show_plot=False, show_image=False)
    axes[0].set_title(f'Grid {group} proportions (max = 1)')
    axes[1].set_title(f'Grid {group} max spots')
    axes[2].set_title(f'Individual cell {group}')
    plt.show()

fig, axes = plt.subplots(ncols=2, figsize=(12,6))
st.pl.gene_plot(grid, gene_symbols='IFNG', ax=axes[0],show_color_bar=False,show_plot=False, show_image=False)
st.pl.gene_plot(adata, gene_symbols='IFNG', ax=axes[1], show_color_bar=False, show_plot=False, show_image=False)
axes[0].set_title(f'Grid IFNG expression')
axes[1].set_title(f'Cell IFNG expression')
plt.show()

fig, axes = plt.subplots(ncols=2, figsize=(12,6))
st.pl.gene_plot(grid, gene_symbols='IFNGR2', ax=axes[0],show_color_bar=False,show_plot=False, show_image=False)
st.pl.gene_plot(adata, gene_symbols='IFNGR2', ax=axes[1], show_color_bar=False, show_plot=False, show_image=False)
axes[0].set_title(f'Grid IFNGR2 expression')
axes[1].set_title(f'Cell IFNGR2 expression')
plt.show()

# Loading the LR databases available within stlearn (from NATMI)
lrs = st.tl.cci.load_lrs(['connectomeDB2020_lit'], species='human')
print(len(lrs))

# Running the analysis #
st.tl.cci.run(grid, lrs,
                  min_spots = 20, #Filter out any LR pairs with no scores for less than min_spots
                  distance=None, # None defaults to spot+immediate neighbours; distance=0 for within-spot mode
                  n_pairs=1000, # Number of random pairs to generate; low as example, recommend ~10,000
                  n_cpus=None, # Number of CPUs for parallel. If None, detects & use all available.
                  )
                  
lr_info = grid.uns['lr_summary'] # A dataframe detailing the LR pairs ranked by number of significant spots.
print(lr_info.shape)
print(lr_info.head(50))

# Showing the rankings of the LR from a global and local perspective.
# Ranking based on number of significant hotspots.
plt.figure()
st.pl.lr_summary(grid, n_top=500)
st.pl.lr_summary(grid, n_top=50, figsize=(10,3))
plt.savefig("lr_rank_plot.png", dpi=500)

### Can adjust significance thresholds.
st.tl.cci.adj_pvals(grid, correct_axis='spot',
                   pval_adj_cutoff=0.05, adj_method='fdr_bh')

best_lr = grid.uns['lr_summary'].index.values[0] # Just choosing one of the top from lr_summary

stats = ['lr_scores', '-log10(p_adjs)', 'lr_sig_scores']
fig, axes = plt.subplots(ncols=len(stats), figsize=(18,6))
for i, stat in enumerate(stats):
    st.pl.lr_result_plot(grid, use_result=stat, use_lr=best_lr, show_color_bar=False, ax=axes[i], show_image=False)
    axes[i].set_title(f'{best_lr} {stat}')

stats = ['lr_scores', '-log10(p_adjs)', 'lr_sig_scores']
fig, axes = plt.subplots(ncols=len(stats), figsize=(18,6))
for i, stat in enumerate(stats):
    st.pl.lr_result_plot(grid, use_result=stat, use_lr="IFNG_IFNGR1", show_color_bar=False, ax=axes[i], show_image=False)
    axes[i].set_title(f'{"IFNG_IFNGR1"} {stat}')

stats = ['lr_scores', '-log10(p_adjs)', 'lr_sig_scores']
fig, axes = plt.subplots(ncols=len(stats), figsize=(18,6))
for i, stat in enumerate(stats):
    st.pl.lr_result_plot(grid, use_result=stat, use_lr="IFNG_IFNGR1", show_color_bar=False, ax=axes[i], show_image=False)
    axes[i].set_title(f'{"IFNG_IFNGR2"} {stat}')

stats = ['lr_scores', '-log10(p_adjs)', 'lr_sig_scores']
fig, axes = plt.subplots(ncols=len(stats), figsize=(18,6))
for i, stat in enumerate(stats):
    st.pl.lr_result_plot(grid, use_result=stat, use_lr="CD274_CD80", show_color_bar=False, ax=axes[i], show_image=False)
    axes[i].set_title(f'{"CD274_CD80"} {stat}')

stats = ['lr_scores', '-log10(p_adjs)', 'lr_sig_scores']
fig, axes = plt.subplots(ncols=len(stats), figsize=(18,6))
for i, stat in enumerate(stats):
    st.pl.lr_result_plot(grid, use_result=stat, use_lr="TGFB1_TGFBR2", show_color_bar=False, ax=axes[i], show_image=False)
    axes[i].set_title(f'{"TGFB1_TGFBR2"} {stat}')
