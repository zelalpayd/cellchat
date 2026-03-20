import scanpy as sc
adata_full = sc.read_h5ad("../brainpart.h5ad")

#subsetting selected development stages; 9th, 12th and 15th week post-fertilization stages
adata_age = adata_full[adata_full.obs['development_stage'].isin(['12th week post-fertilization stage', '15th week post-fertilization stage', '9th week post-fertilization stage'])].copy()  

#filtering genes expressed in at least 10 cells
sc.pp.filter_genes(adata_age, min_cells=10)

#Changing the EnsembleIDs to gene names (data includes the gene names)
adata_age.var_names = adata_age.var["Gene"]

adata_age.write("brainpart_subset.h5ad")
