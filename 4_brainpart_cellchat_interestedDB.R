library(CellChat)
library(Seurat)
library(patchwork)
library(ggplot2)
source("4E_create_cell_chat_object.R")

## Uploading CellChat Database
CellChatDB <- CellChatDB.human
#Subsetting the database to interested database, X can be "Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact" or "Non-protein Signaling"
CellChatDB.X<- subsetDB(CellChatDB, search = "X")
dplyr::glimpse(CellChatDB.cell_cell$interaction)

## Splitting the data into pcws
seurat_object <- readRDS("brainpart_subset.rds")

#dividing the seurat object into a list of seurat objects based on the development stage
seurat.list <- SplitObject(seurat_object, split.by = "development_stage")

## Creating CellChat objects 
# Loading the function and creating all three CellChat objects for cell-cell interactions
source("4E_create_cell_chat_object.R")

# Create 9pcw
cellchat_9pcw <- create_cellchat_stage_safe(
  seurat_obj = seurat.list[["9th week post-fertilization stage"]],
  stage_name = "9pcw",
  group.by = "CellClass",
  min.cells = 10,
  CellChatDB = CellChatDB.X
)

# 12pcw
cellchat_12pcw <- create_cellchat_stage_safe(
  seurat_obj = seurat.list[["12th week post-fertilization stage"]],
  stage_name = "12pcw",
  group.by = "CellClass",
  min.cells = 10,
  CellChatDB = CellChatDB.X
)

#15pcw
cellchat_15pcw <- create_cellchat_stage_safe(
  seurat_obj = seurat.list[["15th week post-fertilization stage"]],
  stage_name = "15pcw",
  group.by = "CellClass",
  min.cells = 10,
  CellChatDB = CellChatDB.X
)

# Merging the three CellChat objects and computing centrality
cellchat.list <- list(
  "9pcw"  = cellchat_9pcw.sub,
  "12pcw" = cellchat_12pcw.sub,
  "15pcw" = cellchat_15pcw.sub
)

cellchat.list <- lapply(cellchat.list, function(x) {
  netAnalysis_computeCentrality(x, slot.name = "netP")
})


# Merging 9 pcw and 12 pcw CellChat objects and computing centrality

cellchat_9_12.list <- list(
    "9pcw"  = cellchat_9pcw.sub,
    "12pcw" = cellchat_12pcw.sub
)

cellchat_9_12.list <- lapply(cellchat_9_12.list, function(x) {
  netAnalysis_computeCentrality(x, slot.name = "netP")

  
# Merging 12 pcw and 15 pcwCellChat objects and computing centrality

  
cellchat_12_15.list <- list(
    "12pcw" = cellchat_12pcw.sub,
    "15pcw"  = cellchat_15pcw.sub
  
)
cellchat_12_15.list <- lapply(cellchat_12_15.list, function(x) {
  netAnalysis_computeCentrality(x, slot.name = "netP")
})


#CellChat Analysis
library(cowplot)
# Compare interactions across developmental stages
gg1 <- compareInteractions(cellchat_merged, show.legend = F, group = c(1,2,3))
gg2 <- compareInteractions(cellchat_merged, show.legend = F, group = c(1,2,3), measure = "weight")
gg1 + gg2

save_plot("brainregion_comparison_interactions_all_stages_X.png", gg1 + gg2, base_width = 8, base_height = 4)
## Number of interactions and Interaction Strenghts for each time point
  #9pcw
groupSize <- as.numeric(table(cellchat_9pcw@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat_9pcw@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat_9pcw@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

#12pcw
groupSize <- as.numeric(table(cellchat_12pcw@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat_12pcw@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat_12pcw@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")


#15pcw
groupSize <- as.numeric(table(cellchat_15pcw@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat_15pcw@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat_15pcw@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")


#9_12pcw
#red increased, blue decreased interactions
netVisual_diffInteraction(cellchat_9_12, weight.scale = T)
netVisual_diffInteraction(cellchat_9_12, weight.scale = T, measure = "weight")


#12_15pcw
netVisual_diffInteraction(cellchat_12_15, weight.scale = T)
netVisual_diffInteraction(cellchat_12_15, weight.scale = T, measure = "weight")

  
## Heatmap showing differential number of interactions or interaction strength among different cell populations across two datasets


# 3. Create the plots for 12_15pcw
gg1_obj <- netVisual_heatmap(cellchat_12_15,
                                 measure = "count", title = "12-15 PCW (Count)")

gg2_obj <- netVisual_heatmap(cellchat_12_15, 
                                 measure = "weight", title = "12-15 PCW (Weight)")


#12_15pcw Count heatmap
library(grid)
library(cowplot)

p1 <- grid.grabExpr(ComplexHeatmap::draw(gg1_obj))

plot_grid(p1)
save_plot("brainregion_heatmap_12_15_count_X.png", p1, base_width = 6, base_height = 5)



#12_15pcw Weight heatmap
p2 <- grid.grabExpr(ComplexHeatmap::draw(gg2_obj))
plot_grid(p2)
save_plot(".brainregion_heatmap_12_15_weight_X.png", p2, base_width = 6, base_height = 5)


  
# 3. Create the plots for 9_12 pcw
gg3_obj <- netVisual_heatmap(cellchat_9_12, measure = "count", title = "9-12 PCW (Count)")

gg4_obj <- netVisual_heatmap(cellchat_9_12, measure = "weight", title = "9-12 PCW (Weight)")


#9_12pcw Count heatmap
library(grid)
library(cowplot)

p3 <- grid.grabExpr(ComplexHeatmap::draw(gg3_obj))

plot_grid(p3)
save_plot("brainregion_heatmap_9_12_count_X.png", p3, base_width = 6, base_height = 5)

#9_12pcw Weight heatmap
p4 <- grid.grabExpr(ComplexHeatmap::draw(gg4_obj))
plot_grid(p4)
save_plot("brainregion_heatmap_9_12_weight_X.png", p4, base_width = 6, base_height = 5)
  

## Compare the major sources and targets in a 2D space

library(scales)
all_idents <- levels(cellchat.list[[1]]@idents)


node.colors <- hue_pal()(length(all_idents))
names(node.colors) <- all_idents

num.link <- sapply(cellchat.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
gg <- list()
for (i in 1:length(cellchat.list)) {
  # Generate the plot
  p <- netAnalysis_signalingRole_scatter(
    cellchat.list[[i]], 
    title = names(cellchat.list)[i], 
    weight.MinMax = weight.MinMax,
    color.use = node.colors
  )
  
  # Add theme adjustments to fix the "tight and long" look
  gg[[i]] <- p + 
    theme(
      aspect.ratio = 1, # Forces a square shape
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
}
combined_plot <- patchwork::wrap_plots(plots = gg, ncol = 2) + 
  plot_layout(guides = "collect")

combined_plot
save_plot("brainregion_signaling_role_scatter_X.png", combined_plot, base_width = 8, base_height = 6)

### Identify the signaling changes of specific cell populations
  
# Scatter plot
gg1 <- netAnalysis_signalingChanges_scatter(
  cellchat_9_12.list, 
  idents.use = "Vascular",
  signaling.exclude = "MIF"
)


gg1
save_plot("brainregion_signaling_changes_vascular_X.png", gg1, base_width = 6, base_height = 5)


## Identify altered sigaling with distinct interaction strength
### Compare the overall information flow of each signaling pathway or ligand-receptor pair for 9_12 pcw
gg1 <- rankNet(cellchat_9_12, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat_9_12, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = F, do.stat = TRUE)
save_plot("brainregion_ranking_9_12_X.png", gg1 + gg2, base_width = 8, base_height = 4)

gg1 + gg2


### Compare the overall information flow of each signaling pathway or ligand-receptor pair for 12_15 pcw   
gg1 <- rankNet(cellchat_12_15, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat_12_15, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = F, do.stat = TRUE)

gg1 + gg2
save_plot("brainregion_ranking_12_15_X.png", gg1 + gg2, base_width = 8, base_height = 4)


### Compare outgoing (or incoming) signaling patterns associated with each cell population
library(ComplexHeatmap)

i = 1
pathway.union <- union(cellchat_9_12.list[[i]]@netP$pathways, 
                       cellchat_9_12.list[[i+1]]@netP$pathways)

# Generate the heatmaps for outgoing signaling of 9_12
ht1 = netAnalysis_signalingRole_heatmap(cellchat_9_12.list[[i]], 
                                        pattern = "outgoing", 
                                        signaling = pathway.union, 
                                        title = names(cellchat_9_12.list)[i], 
                                        width = 5, height = 6)

ht2 = netAnalysis_signalingRole_heatmap(cellchat_9_12.list[[i+1]], 
                                        pattern = "outgoing", 
                                        signaling = pathway.union, 
                                        title = names(cellchat_9_12.list)[i+1], 
                                        width = 5, height = 6)

library(grid)
library(cowplot)

# Convert heatmaps to grobs (graphical objects)
gb1 = grid.grabExpr(draw(ht1))
gb2 = grid.grabExpr(draw(ht2))

# Plot side-by-side
plot_grid(gb1, gb2, ncol = 2)
save_plot("brainregion_outgoing_signaling_role_heatmap_9_12_X.png", plot_grid(gb1, gb2, ncol = 2), base_width = 10, base_height = 6)

  
# Generate the heatmaps for incoming signaling of 9_12
ht5 = netAnalysis_signalingRole_heatmap(cellchat_9_12.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(cellchat_9_12.list)[i], width = 5, height = 6, color.heatmap = "GnBu")
ht6 = netAnalysis_signalingRole_heatmap(cellchat_9_12.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(cellchat_9_12.list)[i+1], width = 5, height = 6, color.heatmap = "GnBu")



# Convert heatmaps to grobs (graphical objects)
gb5 = grid.grabExpr(draw(ht5))
gb6 = grid.grabExpr(draw(ht6))

# Plot side-by-side
plot_grid(gb5, gb6, ncol = 2)
save_plot("brainpart_signaling_role_heatmap_incoming_9_12_X.png", plot_grid(gb5, gb6, ncol = 2), base_width = 10, base_height = 6)


#12_15
library(ComplexHeatmap)

i = 1
pathway.union <- union(cellchat_12_15.list[[i]]@netP$pathways, 
                       cellchat_12_15.list[[i+1]]@netP$pathways)

# Generate the heatmaps for outgoing signaling of 12_15pcw
ht3 = netAnalysis_signalingRole_heatmap(cellchat_12_15.list[[i]], 
                                        pattern = "outgoing", 
                                        signaling = pathway.union, 
                                        title = names(cellchat_12_15.list)[i], 
                                        width = 5, height = 6)

ht4 = netAnalysis_signalingRole_heatmap(cellchat_12_15.list[[i+1]], 
                                        pattern = "outgoing", 
                                        signaling = pathway.union, 
                                        title = names(cellchat_12_15.list)[i+1], 
                                        width = 5, height = 6)


# Convert heatmaps to grobs (graphical objects)
gb3 = grid.grabExpr(draw(ht3))
gb4 = grid.grabExpr(draw(ht4))

# Plot side-by-side
plot_grid(gb3, gb4, ncol = 2)
save_plot("brainpart_signaling_role_heatmap_12_15_X.png", plot_grid(gb3, gb4, ncol = 2), base_width = 10, base_height = 6)

# Generate the heatmaps for incoming signaling of 12_15pcw
ht7 = netAnalysis_signalingRole_heatmap(cellchat_12_15.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(cellchat_12_15.list)[i], width = 5, height = 6, color.heatmap = "GnBu")
ht8 = netAnalysis_signalingRole_heatmap(cellchat_12_15.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(cellchat_12_15.list)[i+1], width = 5, height = 6, color.heatmap = "GnBu")


library(grid)
library(cowplot)

# Convert heatmaps to grobs (graphical objects)
gb7 = grid.grabExpr(draw(ht7))
gb8 = grid.grabExpr(draw(ht8))

# Plot side-by-side
plot_grid(gb7, gb8, ncol = 2)
save_plot("brainpart_signaling_role_heatmap_incoming_12_15_X.png", plot_grid(gb7, gb8, ncol = 2), base_width = 10, base_height = 6)



# Compare the signaling gene expression distribution between different datasets

cellchat_merged@meta$datasets = factor(cellchat_merged@meta$datasets, levels = c("9pcw", "12pcw","15pcw")) # set factor level
plotGeneExpression(cellchat_merged, signaling = "LAMININ", split.by = "datasets", colors.ggplot = T, type = "violin")

#> Adding another scale for y, which will replace the existing scale.
save_plot("brainpart_laminin_expression_X.png", last_plot(), base_width = 6, base_height = 5)


#Save the data
save(cellchat.list, file = "cellchat.list_midbrain.RData")
save(cellchat_merged, file = "cellchat_merged_midbrain.RData")

save(cellchat_9_12, file = "cellchat_9_12_midbrain.RData")
save(cellchat_12_15, file = "cellchat_12_15_midbrain.RData")









  
