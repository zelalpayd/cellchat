import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

# 1. Load the h5ad file
adata = sc.read_h5ad('brainpart.h5ad')

# 2. Basic Preprocessing / Normalization
# Store raw counts in a separate layer before normalization
adata.layers["counts"] = adata.X.copy() 

# Normalize to 10,000 reads per cell
sc.pp.normalize_total(adata, target_sum=1e4)

# Logarithmize the data to stabilize variance and make it more suitable for downstream analysis
sc.pp.log1p(adata)


# Box plot showing median and quartiles
genes_to_plot = ['MMP15', 'ADAMTS7']  

fig, axes = plt.subplots(1, len(genes_to_plot), figsize=(7*len(genes_to_plot), 6))
if len(genes_to_plot) == 1:
    axes = [axes]

for idx, gene in enumerate(genes_to_plot):
    gene_matches = adata.var[adata.var['Gene'] == gene]
    
    if len(gene_matches) > 0:
        gene_idx = gene_matches.index[0]
        expr = adata[:, gene_idx].X
        if hasattr(expr, "toarray"):
            expr = expr.toarray().flatten()
        else:
            expr = np.array(expr).flatten()
        
        plot_data = pd.DataFrame({
            'Expression': expr,
            'Development Stage': adata.obs['development_stage'].values
        })
        
        # Box plot
        sns.boxplot(data=plot_data, x='Development Stage', y='Expression', ax=axes[idx])
        axes[idx].set_title(f'{gene} Expression by Development Stage')
        axes[idx].set_xlabel('Development Stage')
        axes[idx].set_ylabel('Normalized Expression')
        axes[idx].tick_params(axis='x', rotation=45)
        # Align x-axis labels horizontally to the right
        for label in axes[idx].get_xticklabels():
            label.set_ha('right')
    else:
        axes[idx].text(0.5, 0.5, f'{gene} not found', ha='center', va='center')

plt.show()
plt.tight_layout()
