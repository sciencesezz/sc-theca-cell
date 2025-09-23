setwd("/mnt/data/home/sarahsczelecki/single-cell/seurat/outputs")

#install.packages("Seurat")
library(Seurat)
#install.packages("ggplot2")
library(ggplot2)
library(sctransform)
library(tidyr)
library(dplyr)
#install.packages("hdf5r")
#install.packages("remotes")
library(hdf5r)
library(remotes)
#remotes::install_github("10xGenomics/loupeR")
#loupeR::setup()
library(loupeR)
library(patchwork)


#import R object
#load("2w3w1y-WT-merged_seurat.Robj")
#load("2w3w1yWT-stroma-subset_seurat.RData")
load("2w-stroma-subset_seurat.RData") #this should be the singlet one, I wrote over the old one
View(subset_seurat)
head(subset_seurat)
tail(subset_seurat)

###########################slingshot - single trajectory#########################
#if (!requireNamespace("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")
#BiocManager::install("slingshot")
library(slingshot)
library(SingleCellExperiment)
library(ggplot2)


# Get the necessary data from Seurat
counts <- GetAssayData(subset_seurat, layer = "counts")
metadata <- subset_seurat@meta.data
Idents(subset_seurat) <- "RNA_snn_res.0.6"
clusters <- Idents(subset_seurat)

# Get dimensional reduction (e.g., UMAP)
# You can also use PCA or tSNE

dimred <- Embeddings(subset_seurat, reduction = "umap")
# Create SingleCellExperiment object
sce <- SingleCellExperiment(
  assays = list(counts = counts),
  colData = metadata,
  reducedDims = list(UMAP = dimred)
)

# Add cluster information
sce$cluster <- clusters

# Specify start.clus if you know your starting cluster
sce <- slingshot(sce, 
                 clusterLabels = 'cluster',
                 reducedDim = 'UMAP',
                 start.clus = '6')  # Replace 'X' with your starting cluster

# Create plotting data
plot_data <- data.frame(
  UMAP1 = reducedDims(sce)$UMAP[,1],
  UMAP2 = reducedDims(sce)$UMAP[,2],
  Pseudotime = slingPseudotime(sce)[,1],
  Cluster = sce$cluster
)

# Create color palette for pseudotime
colors <- colorRampPalette(c("#FFFF9F", "#FF9F9F", "#9F9FFF"))(100)

# Extract and format curve data
curves <- slingCurves(sce)
curve_data <- data.frame(
  UMAP1 = curves[[1]]$s[,1],  # 's' contains the curve coordinates
  UMAP2 = curves[[1]]$s[,2]
)

# Plot with ggplot2
ggplot(plot_data, aes(x = UMAP1, y = UMAP2)) +
  geom_point(aes(color = Pseudotime), size = 1) +
  scale_color_gradientn(colors = colors, na.value = "grey50") +
  geom_path(data = curve_data, 
            aes(x = UMAP1, y = UMAP2),
            color = "black", 
            linewidth = 1) +
  theme_minimal() +
  labs(title = "Pseudotime Trajectory Analysis",
       x = "UMAP1",
       y = "UMAP2") +
  theme(legend.position = "right")

######################pseudotime slingshot multiple trajectories####################################
# Load required libraries
library(grDevices)
library(graphics)

# Previous setup code remains the same
sce <- as.SingleCellExperiment(subset_seurat)
reducedDims(sce)$UMAP <- subset_seurat[["umap"]]@cell.embeddings
clusters <- as.factor(Idents(subset_seurat))
colData(sce)$GMM <- clusters


set.seed(123)
sds <- slingshot(sce, 
                 reducedDim = 'UMAP',
                 clusterLabels = colData(sce)$GMM,
                 start.clus = '6',
                 stretch = 0.8,
                 smoother = 'smooth.spline',
                 shrink = 1.0)

pt <- slingPseudotime(sds)
n_lineages <- ncol(pt)


colors_pt <- colorRampPalette(c("grey", "lightblue", "navyblue"))(100)

# Individual lineage plots
for(i in 1:n_lineages) {
  
  # Open a square PNG device
  png(filename = paste0("Lineage_", i, ".png"), width = 1200, height = 1200, res = 150, family = "Arial")
  
  # Set margins: bottom, left, top, right
  par(mar = c(6,6,5,10))
  par(xpd = NA)  # allow legend outside plot region
  
  current_pt <- pt[,i]
  pt_range <- range(current_pt, na.rm = TRUE)
  
  # Plot cells colored by pseudotime
  plot(reducedDims(sce)$UMAP,
       col = colors_pt[cut(current_pt, breaks = 100)],
       pch = 16,
       cex = 0.8,
       xlab = "UMAP 1",
       ylab = "UMAP 2",
       main = paste("Lineage", i),
       cex.lab = 2,        # 20 pt equivalent
       cex.axis = 1.9,     # slightly smaller than labels
       cex.main = 2,       # 20 pt equivalent
       family = "Arial")
  
  # Add trajectories
  lines(SlingshotDataSet(sds), lwd = 2, col = "black", type = "curves")
  lines(SlingshotDataSet(sds), lwd = 3, col = "red", type = "curves", linInd = i)
  
  # Add pseudotime legend
  legend_image <- as.raster(matrix(rev(colors_pt), ncol = 1))
  rasterImage(legend_image,
              xleft = par("usr")[2] + 0.1 * diff(par("usr")[1:2]),
              ybottom = par("usr")[3],
              xright = par("usr")[2] + 0.2 * diff(par("usr")[1:2]),
              ytop = par("usr")[4])
  
  # Numeric labels for legend — shifted slightly right
  text(x = par("usr")[2] + 0.27 * diff(par("usr")[1:2]),   # shifted from 0.22 → 0.24
       y = par("usr")[4],
       labels = round(pt_range[2],1),
       cex = 2,
       family = "Arial")
  text(x = par("usr")[2] + 0.27 * diff(par("usr")[1:2]),
       y = par("usr")[3],
       labels = round(pt_range[1],1),
       cex = 2,
       family = "Arial")
  text(x = par("usr")[2] + 0.27 * diff(par("usr")[1:2]),
       y = mean(par("usr")[3:4]),
       labels = "Pseudotime",
       srt = 90,
       cex = 2,
       family = "Arial")
  
  dev.off()
}
################TEST STARTS
# Define your cluster colors manually (order should match levels(clusters))
cluster_colors <- c(
  "0" = "#F8766D",  
  "1" = "#D98900",  
  "2" = "#A99C00",  
  "3" = "#64B200",  
  "4" = "#00BD5B",
  "5" = "#00BFA4", 
  "6" = "#00BADE", 
  "7" = "#00A0FF", 
  "8" = "#B385FF",
  "9" = "#EE61EA")

# Make sure clusters are factors with correct levels
clusters <- factor(colData(sce)$GMM)

# Summary plot — square, matching styling
png(filename = "All_Trajectories.png", width = 1200, height = 1200, res = 150, family = "Arial")

par(mar = c(6,6,5,10))
par(xpd = NA)

# Plot using your colors
plot(reducedDims(sce)$UMAP,
     col = cluster_colors[as.character(clusters)],  # <- explicit mapping
     pch = 16,
     cex = 0.8,
     xlab = "UMAP 1",
     ylab = "UMAP 2",
     main = "All Trajectories",
     cex.lab = 2,
     cex.axis = 1.9,
     cex.main = 2,
     family = "Arial")

# Add trajectories
lines(SlingshotDataSet(sds), lwd = 2, col = "black", type = "curves")

# Add legend with same colors
legend(
  x = par("usr")[2] + 0.1 * diff(par("usr")[1:2]),
  y = par("usr")[4] - 0.05 * diff(par("usr")[3:4]),
  legend = levels(clusters),
  col = cluster_colors[levels(clusters)],  # <- same mapping
  pch = 16,
  cex = 2,
  title = "Clusters",
  title.adj = 0
)

dev.off()


##this code now exports each as an individual file, so I can make the figures better

###########################Trajectories monocle3#########################
library(monocle3)
library(Matrix)
#remotes::install_github('satijalab/seurat-wrappers')
library(SeuratWrappers)

#convert seurat object to monocle3 celldataset - seurat has a wrapper function
subset_cds <- as.cell_data_set(subset_seurat)

#####################Using Seurat Clusters######################################
#to get cell metadata, different than seurat to access, but same data
colData(subset_cds)
rowData(subset_cds)
rowData(subset_cds) <- DataFrame(gene_short_name = rownames(subset_cds))

#use clustering I did in seurat, for monocle3
#saving cells to 1 partition
#saving clusters perfomed in seurat into cds
#saving UMAP embeddings/coordinates performed in seurat into cds

#assign partitions - assign all cells to one partition (not sure why...I think
#monocle3 makes partitions and clusters, and partitions are just mega clusters, so
#assigning everything to 1 partition arbitrarily)

recreate.partition <- c(rep(1, length(subset_cds@colData@rownames)))
names(recreate.partition) <- subset_cds@colData@rownames
recreate.partition <- as.factor(recreate.partition)
subset_cds@clusters$UMAP$partitions <- recreate.partition

#assign the cluster information
#these two steps save the cluster information from the seurat object into the cluster
#and partition slots of the cds object, instead of using monocle3 clustering
#this weird syntax is how to access the different partition layers

list_cluster <- subset_seurat@active.ident
subset_cds@clusters$UMAP$clusters <- list_cluster

#save UMAP coordinates from seurat object inside the cds object
subset_cds@int_colData@listData$reducedDims$UMAP <- subset_seurat@reductions$umap@cell.embeddings

#plot to visualise - looks good, looks like seurat object

cluster.before.trajectory_sub <- plot_cells(subset_cds,
                                            color_cells_by = "cluster", 
                                            label_groups_by_cluster = FALSE,
                                            group_label_size = 5) + 
  theme(legend.position = "right") + 
  ggtitle("UMAP - Leiden")

cluster.before.trajectory_sub

#learn trajectory graph
subset_cds <- learn_graph(subset_cds, use_partition = FALSE)

plot_cells(subset_cds, 
           color_cells_by = "cluster",
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE, 
           label_roots = FALSE, 
           label_leaves = FALSE,
           group_label_size = 5) + 
  theme(legend.position = "right",  # Control legend position
        legend.title = element_text(size = 12),  # Control legend title size
        legend.text = element_text(size = 10))   # Control legend text size


#order the cells in pseudotime (earlier state, smaller pseudotime)
#based on my loupebrowser guess...clusterXXX is preTC??...starts as stroma mixed and then
#bifucates as it differentiates?

#root nodes - computationally?
# specifying root cells: `root_pr_nodes` argument - check the principal points
#plot_cells(cds,
#          color_cells_by = "cluster",
#         label_cell_groups=FALSE,
#        label_groups_by_cluster=FALSE,
#       label_leaves=FALSE,
#      label_branch_points=FALSE,
#     label_principal_points = TRUE,       # set this to TRUE
#    graph_label_size=3)


#this makes it the correct structure for the order_cells - guessing the root cells as a cluster
root_cells <- colnames(subset_cds)[clusters(subset_cds) == 6]

subset_cds <- order_cells(subset_cds, reduction_method = "UMAP", 
                          root_cells = root_cells)

library(ggplot2)

# Create the plot object
p <- plot_cells(subset_cds, 
                color_cells_by = "pseudotime", 
                label_groups_by_cluster = FALSE,
                show_trajectory_graph = TRUE,
                trajectory_graph_color = "grey80", 
                trajectory_graph_segment_size = 1.5,
                label_branch_points = FALSE, 
                label_roots = FALSE, 
                label_leaves = FALSE, 
                cell_size = 0.7, 
                graph_label_size = 6) 

# Optionally adjust font sizes with theme()
p <- p + theme(
  text = element_text(size = 20, family = "Arial"),
  axis.title = element_text(size = 20),
  axis.text = element_text(size = 20),
  legend.title = element_text(size = 20),
  legend.text = element_text(size = 20)
)
print(p)

# Save the plot to a file
ggsave(filename = "pseudotime_plot_monocle3.png",   # file name
       plot = p,                          # plot object
       width = 8, height = 8,             # width and height in inches
       dpi = 800)                          # resolution

###################genes in pseudotime######################################
library(monocle3)
library(ggplot2)
library(reshape2)

#extracting pseduotime from the root note from above...can change to investigate

# Genes of interest
#genes_of_interest <- c("Enpep", "Anpep", "Tcf21", "Wt1", "Gli1")
genes_of_interest <- c("Sod2", "Ppargc1a", "Sod1", "Sod3", "Ppargc1b")
# Extract pseudotime values (named vector, names = cells)
pt <- pseudotime(subset_cds)

# Extract expression matrix (genes x cells)
expr_mat <- log1p(exprs(subset_cds[genes_of_interest, ]))

# Transpose to cells x genes
expr_data <- as.data.frame(t(expr_mat))

# Add pseudotime and metadata directly (rownames are cell IDs)
expr_data$pseudotime <- pt[rownames(expr_data)]

meta <- as.data.frame(colData(subset_cds))
expr_data$Age <- meta[rownames(expr_data), "Age"]

# Reshape to long format
expr_long <- melt(expr_data,
                  id.vars = c("pseudotime", "Age"),
                  variable.name = "gene",
                  value.name = "expression")

# Plot by Age
p1 <- ggplot(expr_long, aes(x = pseudotime, y = expression, color = gene)) +
  geom_point(size = 1.0, alpha = 0.25) +
  geom_smooth(se = FALSE, method = "loess", size = 1.2) +
  facet_wrap(~ Age) +
  theme_classic() +
  theme(
    strip.text = element_text(size = 18, face = "bold"),
    axis.text.y  = element_text(size = 18),
    axis.text.x  = element_text(size = 18),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    plot.title   = element_text(hjust = 0.5, face = "bold", size = 18),
    legend.text  = element_text(size = 18),
    legend.title = element_text(size = 18, face = "bold")) +
  labs(title = "Gene expression across pseudotime by Age")

print(p1)

ggsave("genes-in-time-by-age.png", plot = p1, width = 9, height = 6, dpi = 800)
