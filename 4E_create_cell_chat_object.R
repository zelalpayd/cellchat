library(Seurat)
library(CellChat)

create_cellchat_stage_safe <- function(seurat_obj, stage_name,
                                       group.by = "CellClass",
                                       min.cells = 10,
                                       CellChatDB) {

  # 1. Filter metadata for cell types with enough cells
  meta.full <- seurat_obj@meta.data
  valid.types <- names(table(meta.full[[group.by]])[table(meta.full[[group.by]]) >= min.cells])
  
  if(length(valid.types) == 0){
    stop(paste0("No cell types with >= ", min.cells, " cells found for stage ", stage_name))
  }
  
  meta.filtered <- meta.full[meta.full[[group.by]] %in% valid.types, ]
  meta.filtered[[group.by]] <- droplevels(meta.filtered[[group.by]])
  
  cat(paste0("Stage ", stage_name, ": ", nrow(meta.filtered), " cells, ", 
             length(unique(meta.filtered[[group.by]])), " cell types\n"))


  # 2. Extract expression matrix
  expr.full <- GetAssayData(seurat_obj, layer = "data")
  
  # Only keep cells present in both metadata and matrix
  cells.keep <- rownames(meta.filtered)
  cells.keep <- cells.keep[cells.keep %in% colnames(expr.full)]
  
  if(length(cells.keep) == 0){
    stop(paste0("No cells left for stage ", stage_name, " after filtering."))
  }
  
  expr.filtered <- expr.full[, cells.keep, drop = FALSE]
  meta.filtered <- meta.filtered[cells.keep, , drop = FALSE]  # match metadata to matrix
  
  # 3. Filter genes to CellChat ligands/receptors
  genes.use <- unique(c(CellChatDB$interaction$ligand, CellChatDB$interaction$receptor))
  genes.keep <- intersect(rownames(expr.filtered), genes.use)
  
  if(length(genes.keep) == 0){
    stop(paste0("No CellChat-relevant genes found in stage ", stage_name))
  }
  
  expr.filtered <- expr.filtered[genes.keep, , drop = FALSE]
  
  # 4. Create CellChat object
  cellchat <- createCellChat(object = expr.filtered,
                             meta = meta.filtered,
                             group.by = group.by)
  cellchat@DB <- CellChatDB
  
  # 5. Preprocessing
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat)
  

  # 6. Compute pathway-level communication
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  
  message(paste0("CellChat object for stage ", stage_name, " created successfully."))
  
  return(cellchat)
}
