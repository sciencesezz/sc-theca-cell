getwd()

library(Seurat)
library(patchwork)
library(CellChat)

setwd("/mnt/data/home/sarahsczelecki/single-cell/seurat/outputs")

load("2w-stroma-subset_seurat.RData")
View(subset_seurat)

#Step 1: prepare required input data for CellChat Analysis (from Seurat Obj)

data.input <- subset_seurat[["RNA"]]$data
#change cluster # to String
Idents(subset_seurat) <- "RNA_snn_res.0.6"
old_ids <- Idents(subset_seurat)
new_ids <- paste0("Cluster", old_ids)
Idents(subset_seurat) <- new_ids
labels <- Idents(subset_seurat)
meta <- data.frame(labels = labels, row.names = names(labels))


#Step 2: Create a CellChat object
cellchat <- createCellChat(object = subset_seurat, group.by = "ident", assay = "RNA")
cellchat@idents <- factor(cellchat@idents,
                          levels = c("Cluster0", "Cluster1", "Cluster2", "Cluster3", "Cluster4", 
                                     "Cluster5", "Cluster6", "Cluster7", "Cluster8", "Cluster9"))

#Step 3: Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.mouse # use CellChatDB.human if running on human data
showDatabaseCategory(CellChatDB)

# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") # use Secreted Signaling
#CellChatDB.use <- subsetDB(CellChatDB, search = "ECM-Receptor", key = "annotation") # use ECM-Receptor
#CellChatDB.use <- subsetDB(CellChatDB, search = "Cell-Cell Contact", key = "annotation") # use Cell-Cell

# set the used database in the object
cellchat@DB <- CellChatDB.use

#Step 4: Preprocessing the expression data for cell-cell communication analysis
# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

#Step 5: Inference of cell-cell communication network
ptm = Sys.time()
cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = 10)

# All inferred ligand–receptor interactions
df.net <- subsetCommunication(cellchat, slot.name = "netP")
#Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)
#Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)
execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))

#visualize the aggregated cell-cell communication network
ptm = Sys.time()
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

#visualise e signaling sent from each cell group
mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

pathways.show <- c("WNT") #choose a pathway, to find available ones goes df.net
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(6,3) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

#chord
# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")

#heatmap
# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")

netAnalysis_contribution(cellchat, signaling = pathways.show)

netVisual_chord_gene(cellchat, sources.use = 4, targets.use = c(4,6,7), lab.cex = 0.5,legend.pos.y = 30)


cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
