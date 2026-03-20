#the tutorial gave errors so Im following chatgpt suggestions
library(zellkonverter)
library(Seurat)

#Reading the subsetted data
sce <- readH5AD("../brainpart_subset.h5ad")

# Converting SCE to Seurat 
S <- SeuratObject::as.Seurat(sce, counts = "X", data = NULL)

# Checking and renaming the default assay if needed
if (DefaultAssay(S) == "originalexp") {
  S <- RenameAssays(S, originalexp = "RNA")
}

#Normalizing the data 
S <- NormalizeData(S, normalization.method = "LogNormalize", scale.factor = 10000)

#Saving the file
saveRDS(S, file = "../cerebellum_subset.rds")
