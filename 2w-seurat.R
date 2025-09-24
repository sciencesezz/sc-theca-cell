###script for making filtered and normalised Robjects#####
getwd()
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
library(DoubletFinder)

setwd("/mnt/data/home/sarahsczelecki/single-cell/seurat/outputs")

#import R object
#load("2w1y1yadeno-final-mito45.RData")
load("2w-final-mito10.RData")
View(merged_seurat)
head(merged_seurat)
tail(merged_seurat)

########################Dimensonality/Normalisation############################
# Extract metadata
metadata <- merged_seurat@meta.data  

#levels <- c("2week_g3", "1year_g2", "1year_g3")
levels <- c("2week_g1", "2week_g2", "2week_g3")

#levels_age <- c("2week", "1year", "1year_adeno")

metadata$sample <- factor(metadata$sample, levels = levels)


#check number of cells per sample

#basic
#metadata %>% 
# ggplot(aes(x=sample, fill=sample)) + 
#geom_bar() +
#theme_classic() +
#theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
#theme(plot.title = element_text(hjust=0.5, face="bold")) +
#ggtitle("NCells")

#larger
metadata %>% 
  ggplot(aes(x = sample, fill = sample)) + 
  geom_bar() +
  theme_classic() +
  theme(
    axis.text.x  = element_text(size = 14, angle = 45, vjust = 1, hjust = 1),
    axis.text.y  = element_text(size = 14),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    plot.title   = element_text(hjust = 0.5, face = "bold", size = 18),
    legend.text  = element_text(size = 14),   # legend item labels
    legend.title = element_text(size = 16),   # legend title
    legend.key.size = unit(1.5, "cm")         # box size for legend keys
  ) +
  ggtitle("NCells") +
  xlab("Sample") +
  ylab("Cell count")

#check number of cells per age group
#metadata %>% 
# ggplot(aes(x=Age, fill=Age)) + 
#geom_bar() +
#theme_classic() +
#theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
#theme(plot.title = element_text(hjust=0.5, face="bold")) +
#ggtitle("NCells")

#larger
metadata %>% 
  ggplot(aes(x = Age, fill = Age)) + 
  geom_bar() +
  theme_classic() +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(), 
    axis.text.y  = element_text(size = 20),
    axis.title.x = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 20, face = "bold"),
    plot.title   = element_text(hjust = 0.5, face = "bold", size = 20),
    legend.text  = element_text(size = 20),   # legend item labels
    legend.title = element_text(size = 20, face = "bold"),   # legend title
    legend.key.size = unit(1.5, "cm")         # box size for legend keys
  ) +
  ggtitle("Number of Cells") +
  xlab("Sample") +
  ylab("Cell Count")


# Visualize the number UMIs/transcripts per cell
metadata %>% 
  ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(), 
    axis.text.y  = element_text(size = 20),
    axis.title.x = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 20, face = "bold"),
    plot.title   = element_text(hjust = 0.5, face = "bold", size = 20),
    legend.text  = element_text(size = 20),   # legend item labels
    legend.title = element_text(size = 20, face = "bold"),   # legend title
    legend.key.size = unit(1.5, "cm")         # box size for legend keys
  ) +
  ylab("Cell density") +
  geom_vline(xintercept = 500) + 
  ggtitle("UMI counts (transcripts) per cell")



# Visualize the distribution of genes detected per cell via histogram
metadata %>% 
  ggplot(aes(color= sample, x=nGene, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(), 
    axis.text.y  = element_text(size = 20),
    axis.title.x = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 20, face = "bold"),
    plot.title   = element_text(hjust = 0.5, face = "bold", size = 20),
    legend.text  = element_text(size = 20),   # legend item labels
    legend.title = element_text(size = 20, face = "bold"),   # legend title
    legend.key.size = unit(1.5, "cm")         # box size for legend keys
  ) +
  scale_x_log10() + 
  geom_vline(xintercept = 300) + 
  ggtitle("Genes detected per cell")

# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI (novelty score)
metadata %>%
  ggplot(aes(x = log10GenesPerUMI, color = sample, fill = sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(), 
    axis.text.y  = element_text(size = 20),
    axis.title.x = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 20, face = "bold"),
    plot.title   = element_text(hjust = 0.5, face = "bold", size = 20),
    legend.text  = element_text(size = 20),   # legend item labels
    legend.title = element_text(size = 20, face = "bold"),   # legend title
    legend.key.size = unit(1.5, "cm")         # box size for legend keys
  ) +
  geom_vline(xintercept = 0.8) + 
  ggtitle("Genes detected per UMI\n(novelty score)")

# Visualize the distribution of mitochondrial gene expression detected per cell
metadata %>% 
  ggplot(aes(color = sample, x = percent.mt, fill = sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(), 
    axis.text.y  = element_text(size = 20),
    axis.title.x = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 20, face = "bold"),
    plot.title   = element_text(hjust = 0.5, face = "bold", size = 20),
    legend.text  = element_text(size = 20),   # legend item labels
    legend.title = element_text(size = 20, face = "bold"),   # legend title
    legend.key.size = unit(1.5, "cm")         # box size for legend keys
  ) +
  geom_vline(xintercept = 0.2) + 
  ggtitle("Mitocondrial gene expression per cell")

# Visualize the correlation between genes detected and number of UMIs and determine 
#whether strong presence of cells with low numbers of genes/UMIs
metadata %>% 
  ggplot(aes(x=nUMI, y=nGene, color=percent.mt)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  theme(
    strip.text = element_text(size = 18, face = "bold"),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(), 
    axis.text.y  = element_text(size = 18),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    plot.title   = element_text(hjust = 0.5, face = "bold", size = 18),
    legend.text  = element_text(size = 18),   # legend item labels
    legend.title = element_text(size = 18, face = "bold"),   # legend title
    legend.key.size = unit(1.5, "cm")         # box size for legend keys
  ) +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~sample) + 
  ggtitle("Correlation between genes\ndetected and # of UMIs per cell")
##Step 3 - normalisation
merged_seurat <- NormalizeData(merged_seurat, 
                               normalization.method = "LogNormalize", 
                               scale.factor = 10000)

## Step 4 - identify highly variable features
#only want to select a few features that exhibit high cell to cell variation
merged_seurat <- FindVariableFeatures(merged_seurat, selection.method = "vst", nfeatures = 2000)

#identify the 10 most highly variable genes
top10 <- head(VariableFeatures(merged_seurat), 10)

#plot variable features with and without labels
vfeat_plot1 <- VariableFeaturePlot(merged_seurat)
LabelPoints(plot = vfeat_plot1, points = top10, repel = TRUE)

#Step 5 - now we scale the data.
merged_seurat <- ScaleData(merged_seurat, features = rownames(merged_seurat))

#take a look at slots in seurat object, above is stored in scale.data slot of seurat object 
str(merged_seurat)

#Step 6 - Perform linear dimensionality reduction
merged_seurat <- RunPCA(merged_seurat, features = VariableFeatures(object = merged_seurat))

#visualise PCA results
print(merged_seurat[["pca"]], dims = 1:5, nfeatures = 5)
DimHeatmap(merged_seurat, dims = 1:10, cells = 500, balanced = TRUE)

#determine dimensionality of the data, to choose which PC to include in downstream analysis
ElbowPlot(merged_seurat)

#Step 7 - clustering, selected 15 dimensions based on elbow plot, arguably could be around 12...
#start with less strict
merged_seurat <- FindNeighbors(merged_seurat, dims = 1:16)

#understanding resolution - higher the number the more clusters, putting multiple checks multiple
#resolutions at once to find the best
merged_seurat <- FindClusters(merged_seurat, resolution = c(0.1, 0.2, 0.3, 0.5, 0.6, 0.7))
View(merged_seurat@meta.data)

#plotting the resolutions - have a look at these...I chose 0.2
DimPlot(merged_seurat, group.by = "RNA_snn_res.0.2", label = TRUE)

#setting identity of clusters
Idents(merged_seurat) <- "RNA_snn_res.0.2"

#nonlinear dimenstionality reduction
merged_seurat <- RunUMAP(merged_seurat, dims = 1:16)
#merged_seurat <- RunTSNE(merged_seurat, dims = 1:20)

#export QCd file to LoupeBrowser, cannot have two files the same name, doesn't overwrite
create_loupe_from_seurat(
  merged_seurat,
  output_dir = "/mnt/data/home/sarahsczelecki/single-cell/seurat/outputs",
  output_name = "2w-final-mito10-sept-doublets", 
  metadata_cols = c("RNA_snn_res.0.1", "RNA_snn_res.0.2", 
                    "RNA_snn_res.0.3", "Age", "sample"))

DimPlot(merged_seurat, reduction = "umap", group.by = "sample")
DimPlot(merged_seurat, reduction = "umap", group.by = "RNA_snn_res.0.2")

#pN - number of generated artificial doublets - default is 25%, performance is largely pN-invariant
### pK value opt ---------------------------------------------------------------
#pK value optimisation (no ground truth) - from github page
sweep_res  <- paramSweep(merged_seurat, PCs = 1:16, sct = FALSE) #this number of dimensions has to match above
sweep_stats <- summarizeSweep(sweep_res, GT = FALSE)
find_pk <- find.pK(sweep_stats)

ggplot(find_pk, aes(pK, BCmetric, group = 1)) +
  geom_point() + 
  geom_line()

#pK value that corresponds to the maxmium BCMetric is the optimal pK value (based on graph)
#in our case it is 0.005
#stored to pK value

pK <- find_pk %>%
  filter(BCmetric == max(BCmetric)) %>%
  select(pK)
pK <- as.numeric(as.character(pK[[1]]))

#homotypic double proportion estimate -----------------------------------------
annotations <- merged_seurat@meta.data$RNA_snn_res.0.2 #annotations are cell clusters
homotypic_prop <- modelHomotypic(annotations)
nExp_poi <- round(0.08*nrow(merged_seurat@meta.data)) #set the % doublets from 10x sheet, exp dependent
nExp_poi.adj <- round(nExp_poi*(1-homotypic_prop))


#nExp - defines the pANN threshold used to make final doublet/singlet predictions
#This value can best be estimated from cell loading densities into the 10X/Drop-Seq device,
#and adjusted according to the estimated proportion of homotypic doublets.
#our 2024 experiment targeted 10000 cells see benchling for number of cells/sample post filtering

##run DoubletFinder---------------------------------------------------------
merged_seurat <- doubletFinder(merged_seurat,
                               PCs = 1:16,
                               pN = 0.25, 
                               pK = pK, 
                               nExp = nExp_poi.adj,
                               reuse.pANN = NULL, 
                               sct = FALSE)

#visualise doublets


p1 <- DimPlot(merged_seurat, reduction = "umap", group.by = "DF.classifications_0.25_0.005_2031") + 
  ggtitle("DoubletFinder") + 
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 18),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 18))
print(p1)
ggsave("doublet-umap.png", plot = p1, width = 7, height = 7, dpi = 800)


#number of singlets and doublets
table(merged_seurat@meta.data$DF.classifications_0.25_0.005_2031)

#doublets needs to be filtered out
merged_seurat <- subset(
  merged_seurat, 
  subset = DF.classifications_0.25_0.005_2031 == "Singlet"
)

create_loupe_from_seurat(
  merged_seurat,
  output_dir = "/mnt/data/home/sarahsczelecki/single-cell/seurat/outputs",
  output_name = "2w-final-mito10-sept-singles", 
  metadata_cols = c("RNA_snn_res.0.1", "RNA_snn_res.0.2", 
                    "RNA_snn_res.0.3", "Age", "sample"))

p1 <- DimPlot(merged_seurat, reduction = "umap", group.by = "DF.classifications_0.25_0.005_2031") + 
  ggtitle("DoubletFinder Filtered") + 
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 18),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 18))
print(p1)
ggsave("singlet-umap.png", plot = p1, width = 7, height = 7, dpi = 800)


#continue with workflow!

#DimPlot(filtered_seurat, reduction = "tsne")

#plot with labels - apparently you can tell if you have batch effects from this
#I don't really understand that so I'm going to have to revisit that. for now I will assume I need to 
#correct for batch effects. 
#now I understand this - if there are any cluster(s) that are distinct for the sample type that means
#that you probably have batch effects. It looks like we dont, so can continue without batch correction.
p1 <- DimPlot(merged_seurat, reduction = "umap", group.by = "RNA_snn_res.0.2", label = TRUE, label.size = 8) + 
  ggtitle("UMAP of Clusters") + 
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 18),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 18))
print(p1)

#Check for Batch effects graphs
#levels_age <- c("2week", "3week", "1year")
#merged_seurat@meta.data$Age <- factor(merged_seurat@meta.data$Age, levels = levels_age)

# Now plot
VlnPlot(merged_seurat, features = c("nUMI", "nGene", "percent.mt"), group.by = "sample")


p2 <- DimPlot(merged_seurat, reduction = "umap", group.by = "sample") + 
  ggtitle("Clusters by Sample") + 
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 18),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 18))
print(p2)

#data annotation
#manual pipeline based on cloupe file
#genes used are found below in heatmap section
library("ggplot2")
# First, make a copy of the cluster column
merged_seurat$cluster_id <- merged_seurat$RNA_snn_res.0.2

# Now create the cell_type column by direct assignment
merged_seurat$cell_type <- NA  # Initialize with NA

# Assign cell types for each cluster
merged_seurat$cell_type[merged_seurat$cluster_id == 0] <- "Fibroblastic Stroma 1"
merged_seurat$cell_type[merged_seurat$cluster_id == 1] <- "Fibroblastic Stroma 2"
merged_seurat$cell_type[merged_seurat$cluster_id == 2] <- "Fibroblastic Stroma 3"
merged_seurat$cell_type[merged_seurat$cluster_id == 3] <- "GCs"
merged_seurat$cell_type[merged_seurat$cluster_id == 4] <- "Theca Cells"
merged_seurat$cell_type[merged_seurat$cluster_id == 5] <- "Epithelial"
merged_seurat$cell_type[merged_seurat$cluster_id == 6] <- "Endothelial 1"
merged_seurat$cell_type[merged_seurat$cluster_id == 7] <- "Fibroblastic Stroma 4"
merged_seurat$cell_type[merged_seurat$cluster_id == 8] <- "T Lymphocytes"
merged_seurat$cell_type[merged_seurat$cluster_id == 9] <- "Fibroblastic Stroma 5"
merged_seurat$cell_type[merged_seurat$cluster_id == 10] <- "Phagocytes"
merged_seurat$cell_type[merged_seurat$cluster_id == 11] <- "Endothelial 2"
merged_seurat$cell_type[merged_seurat$cluster_id == 12] <- "B Lymphocytes"
merged_seurat$cell_type[merged_seurat$cluster_id == 13] <- "Pericytes"

p1 <- DimPlot(merged_seurat, reduction = "umap", group.by = "cell_type", label = FALSE) + 
  ggtitle("Cell Type") + 
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 18),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 18))
print(p1)
ggsave("cell-type-umap-small.png", plot = p1, width = 7, height = 5, dpi = 800)

p2 <- DimPlot(merged_seurat, reduction = "umap", group.by = "cell_type", label = FALSE) + 
  ggtitle("Cell Type") + 
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 18),
    legend.position = "none"
  )
print(p2)
ggsave("cell-type-umap-small-xleg.png", plot = p2, width = 7, height = 7, dpi = 800)

p1 <- DimPlot(merged_seurat, reduction = "umap", group.by = "sample", label = FALSE) + 
  ggtitle("Sample") + 
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 18),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 18))
print(p1)
ggsave("sample-umap-small.png", plot = p1, width = 7, height = 5, dpi = 800)

p2 <- DimPlot(merged_seurat, reduction = "umap", group.by = "sample", label = FALSE) + 
  ggtitle("Sample") + 
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 18),
    legend.position = "none"
  )
print(p2)
ggsave("age-umap-small-xleg.png", plot = p2, width = 7, height = 7, dpi = 800)

#make an Robject to the input into other parts of the pipeline
# Assuming your Seurat object is called "seurat_obj"
save(merged_seurat, file = "2w-WT-merged_seurat.RData")

# Make a table of number of cells in each cell type by age
#table(merged_seurat$cell_type, merged_seurat$Age)
table(merged_seurat$cell_type)
# First, create the count table with both cell type and age
cell_type_age_counts <- as.data.frame(table(merged_seurat$cell_type))
colnames(cell_type_age_counts) <- c("CellType", "Count")

# View the result
print(cell_type_age_counts)

# Compute percentage within each cell type (so each cell type bar sums to 100%)
#library(dplyr)
cell_type_age_counts <- cell_type_age_counts %>%
  #dplyr::group_by(CellType) %>%
  mutate(Percentage = (Count / sum(Count)) * 100) %>%
  ungroup()

cell_type_age_counts <- cell_type_age_counts %>%
  dplyr::group_by(CellType)


# Order cell types alphabetically
cell_type_age_counts$CellType <- factor(cell_type_age_counts$CellType, 
                                        levels = sort(unique(cell_type_age_counts$CellType)))

# Create stacked bar plot
p <- ggplot(cell_type_age_counts, aes(x = reorder(CellType, -Percentage), y = Percentage, fill = CellType)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = ifelse(Percentage > 0, sprintf("%.1f%%", Percentage), "")), 
            position = position_stack(vjust = 0.5), 
            size = 3.2, color = "black") +
  theme_minimal() +
  labs(title = "Cell Type Composition", 
       x = "Cell Type", y = "Percentage (%)") +
  theme(
    axis.text.x = element_text(size = 15, angle = 45, hjust = 1),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 18),
    legend.position = "none",
    plot.margin = margin(t = 10, r = 20, b = 10, l = 45)  # top, right, bottom, left
  ) + 
  ggtitle("Cell Type Composition")
print(p)

ggsave("cell-type-comp.png", plot = p, width = 8, height = 7, dpi = 800)

#find optimal clusters for the dataset
#if (!require("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("clustree")
library(clustree)
# Find markers for each cluster - only.pos TRUE is only upregulated, change to false to have both
#min.pct = min pecent, filters out genes expressed in fewer than 25% of cells in a cluster
#logfc.threshold - filters out genes with a log2foldchange of less than 0.25

#this takes a very long time (even on the HPC), so pour a coffee..apparently presto makes it faster...

#install.packages('devtools')
#devtools::install_github('immunogenomics/presto')
library(presto)
#BiocManager::install("MAST")
library(MAST)

#cluster_markers <- FindAllMarkers(merged_seurat, test.use = "MAST", only.pos = TRUE)
#cluster_markers_mast <- cluster_markers
#cluster_markers_med <- FindAllMarkers(merged_seurat, only.pos = TRUE, min.pct = 0.50, logfc.threshold = 0.50)

cluster_markers_75high <- FindAllMarkers(merged_seurat, test.use = "MAST", 
                                         only.pos = TRUE, min.pct = 0.75, logfc.threshold = 2.0)

View(cluster_markers_75high)
write.csv(cluster_markers_75high, "2w_cluster_markers_75high.csv", row.names = FALSE)


#I wonder if you can annotate this heatmap, so that you only include your markers of interest.

cluster_markers_75high %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 2) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5
DoHeatmap(merged_seurat, features = top5$gene) + NoLegend()

#define genes by functional categories to annotate with your own selected markers
bcells <- "Cd79a"
tcells <- "Cd3g"
endothelial <- c("Cd93", "Cdh5", "Pecam1")
epithelial <- c("Krt7", "Krt18", "Lgals7", "Upk3b", "Muc16")
pericytes <- c("Rgs5", "Notch3", "Rgs4")
fibroblast <- c("Col1a1", "Dcn", "Mfap4", "Mgp", "Col1a2")
granulosa <- c("Amh", "Cyp19a1", "Inha", "Hsd17b1", "Inhbb", "Fshr")
phagocyte <- "C1qa"
theca <- c("Aldh1a1", "Cyp11a1", "Cyp17a1")
gene_markers <- c(bcells, endothelial, epithelial, fibroblast, granulosa, pericytes,
                  phagocyte, theca, tcells)  # Fixed: changed "phagocytes" to "phagocyte"

# Create plot object with cell type labels
p <- DoHeatmap(merged_seurat, 
               features = gene_markers, 
               size = 5, 
               angle = 45,
               group.by = "cell_type",
               group.bar = TRUE,
               group.bar.height = 0.05, 
               label = FALSE) +
  theme(
    text = element_text(size = 16),              # Overall base font size
    axis.text.x = element_blank(),  # X-axis labels
    axis.text.y = element_text(size = 18),       # Y-axis labels
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),# Axis titles
    plot.title = element_blank(),        # Plot title
    legend.text = element_text(size = 18),       # Legend labels
    legend.title = element_text(size = 18, face = "bold"),
    legend.box.margin = margin(0, 0, 0, -10), 
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5, unit = "pt"), 
    legend.position = "right", 
    legend.justification = "top") +
  guides(fill = guide_colorbar()) +
  guides(color = "none")
print(p)

# Save the plot
ggsave("heatmap-markers-2w.png", p, width = 10, height = 8, dpi = 1200)

#---------------------------------subset for stromal cells of interest------------
#subset the seurat object
Idents(merged_seurat) <- "RNA_snn_res.0.2"
subset_seurat <- subset(merged_seurat, idents = c("0","1","2","4","7","9"))

my_colors <- c("Fibroblastic Stroma 1" = "#53B400", 
               "Fibroblastic Stroma 2" = "#00BC56",
               "Fibroblastic Stroma 3" = "#00C094",
               "Theca Cells" = "#FF66A8", 
               "Fibroblastic Stroma 4" = "#00BFC4",
               "Fibroblastic Stroma 5" = "#00B6EB")
DimPlot(subset_seurat, group.by = "cell_type", cols = my_colors)

# Re-run the analysis pipeline on the subset
subset_seurat <- NormalizeData(subset_seurat, 
                               normalization.method = "LogNormalize", 
                               scale.factor = 10000)
subset_seurat <- FindVariableFeatures(subset_seurat, selection.method = "vst", nfeatures = 2000)
subset_seurat <- ScaleData(subset_seurat, features = rownames(merged_seurat))
subset_seurat <- RunPCA(subset_seurat, features = VariableFeatures(object = merged_seurat))
ElbowPlot(subset_seurat)

subset_seurat <- FindNeighbors(subset_seurat, dims = 1:16) #16 based in elbow plot
subset_seurat <- FindClusters(subset_seurat, resolution = c(0.1, 0.3, 0.5, 0.6, 0.8, 0.7, 1))
DimPlot(subset_seurat, group.by = "RNA_snn_res.0.6", label = TRUE)
#setting identity of clusters
Idents(subset_seurat) <- "RNA_snn_res.0.6"

#nonlinear dimensionality reduction
subset_seurat <- RunUMAP(subset_seurat, dims = 1:16) #not sure if I need to also set this to the same as above check
#merged_seurat <- RunTSNE(merged_seurat, dims = 1:20)

p1 <- DimPlot(subset_seurat, reduction = "umap", group.by = "RNA_snn_res.0.6", label = FALSE) + 
  ggtitle("Subset Clusters") + 
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 18),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 18))
print(p1)
ggsave("subset-umap-small.png", plot = p1, width = 7, height = 7, dpi = 800)

p2 <- DimPlot(subset_seurat, reduction = "umap", group.by = "RNA_snn_res.0.6", label = FALSE) + 
  ggtitle("Cell Type") + 
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 18),
    legend.position = "none"
  )
print(p2)
ggsave("subset-umap-small-xleg.png", plot = p2, width = 7, height = 7, dpi = 800)

p1 <- DimPlot(subset_seurat, reduction = "umap", group.by = "Age", label = FALSE) + 
  ggtitle("Age") + 
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 18),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 18))
print(p1)
ggsave("subset-age-small.png", plot = p1, width = 7, height = 7, dpi = 800)

p2 <- DimPlot(subset_seurat, reduction = "umap", group.by = "Age", label = FALSE) + 
  ggtitle("Age") + 
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 18),
    legend.position = "none"
  )
print(p2)
ggsave("subset-umap-small-xleg.png", plot = p2, width = 7, height = 7, dpi = 800)
#export QCd file to LoupeBrowser
create_loupe_from_seurat(
  subset_seurat,
  output_dir = "/mnt/data/home/sarahsczelecki/single-cell/seurat/outputs",
  output_name = "2w-stroma-subset-singlets")

#make a seurat subset object to load into different aspects of the pipeline
save(subset_seurat, file = "2w-stroma-subset_seurat.RData")

#markers of subset seurat
cluster_markers_subset <- FindAllMarkers(subset_seurat, test.use = "MAST", 
                                         only.pos = TRUE, min.pct = 0.75, logfc.threshold = 1.0)

View(cluster_markers_subset)
write.csv(cluster_markers_subset, "2w_cluster_markers_subset.csv", row.names = FALSE)


#I wonder if you can annotate this heatmap, so that you only include your markers of interest.

cluster_markers_subset %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5

# Create plot object with cell type labels
p <- DoHeatmap(subset_seurat, 
               features = top5$gene, 
               size = 5, 
               angle = 45,
               group.by = "RNA_snn_res.0.6",
               group.bar = TRUE,
               group.bar.height = 0.05, 
               label = FALSE) +
  theme(
    text = element_text(size = 16),              # Overall base font size
    axis.text.x = element_blank(),  # X-axis labels
    axis.text.y = element_text(size = 18),       # Y-axis labels
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),# Axis titles
    plot.title = element_blank(),        # Plot title
    legend.text = element_text(size = 18),       # Legend labels
    legend.title = element_text(size = 18, face = "bold"),
    legend.box.margin = margin(0, 0, 0, -10), 
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5, unit = "pt"), 
    legend.position = "right", 
    legend.justification = "top") +
  guides(fill = guide_colorbar()) +
  guides(color = "none")
print(p)

# Save the plot
ggsave("heatmap-markers-2w-subset.png", p, width = 10, height = 11, dpi = 1200)

##########################GO/GSEA/KEGG##################################


Idents(subset_seurat) <- "RNA_snn_res.0.6"

subset_cluster_markers <- FindAllMarkers(subset_seurat, 
                                         only.pos = TRUE, 
                                         test.use = "MAST",
                                         min.pct = 0.5, logfc.threshold = 2.0)

Idents(subset_seurat) <- "Age"

subset_age_markers <- FindAllMarkers(subset_seurat, 
                                     only.pos = TRUE, 
                                     test.use = "MAST",
                                     min.pct = 0.7, logfc.threshold = 0.25)

markers_adeno_vs_1year <- FindMarkers(subset_seurat, 
                                      ident.1 = "1year_adeno", 
                                      ident.2 = "1year", 
                                      logfc.threshold = 0,   # Only return genes with logFC > 0.25
                                      min.pct = 0.1,            # Gene must be expressed in at least 10% of cells
                                      test.use = "MAST")


markers_adeno_vs_2week <- FindMarkers(subset_seurat, 
                                      ident.1 = "1year_adeno", 
                                      ident.2 = "2week", 
                                      logfc.threshold = 0,   # Only return genes with logFC > 0.25
                                      min.pct = 0.1,            # Gene must be expressed in at least 10% of cells
                                      test.use = "MAST")

markers_1year_vs_2week <- FindMarkers(subset_seurat, 
                                      ident.1 = "1year", 
                                      ident.2 = "2week", 
                                      logfc.threshold = 0,   # Only return genes with logFC > 0.25
                                      min.pct = 0.1,            # Gene must be expressed in at least 10% of cells
                                      test.use = "MAST")


View(subset_cluster_markers)
View(subset_age_markers)
write.csv(subset_cluster_markers, "subset_2week1year_cluster_markers-025logFC.csv", row.names = FALSE)
write.csv(subset_age_markers, "subset_2week1year_age_markers-025logFC.csv", row.names = FALSE)

##cell cycle check on subset data
#use built in gene markers
library(stringr)
s.genes <- Seurat::cc.genes$s.genes
g2m.genes <- Seurat::cc.genes$g2m.genes

#convert to mouse gene case
s.genes <- str_to_title(tolower(s.genes))
g2m.genes <- str_to_title(tolower(g2m.genes))

subset_seurat_cc <- CellCycleScoring(subset_seurat, s.features = s.genes, 
                                     g2m.features = g2m.genes, set.ident = TRUE)

DimPlot(subset_seurat_cc, group.by = "Phase", reduction = "umap") + 
  ggtitle("Cell Cycle Phase")
DimPlot(subset_seurat_cc, group.by = "RNA_snn_res.0.6", reduction = "umap", label = TRUE)
#####################Differential gene expressionv viz between ages################
install.packages("BiocManager")
BiocManager::install("EnhancedVolcano")
library("EnhancedVolcano")


markers_adeno_vs_1year <- markers_adeno_vs_1year %>%
  mutate(
    negLog10P = -log10(p_val_adj),
    negLog10P = ifelse(
      is.infinite(negLog10P),
      max(negLog10P[is.finite(negLog10P)], na.rm = TRUE) + 1,
      negLog10P
    ),
    regulation = case_when(
      avg_log2FC > 0 & p_val_adj < 0.05 ~ "Upregulated",
      avg_log2FC < 0 & p_val_adj < 0.05 ~ "Downregulated",
      TRUE ~ "Not significant"
    )
  )

vol1 <- EnhancedVolcano(markers_adeno_vs_1year,
                        lab = rownames(markers_adeno_vs_1year),
                        x = 'avg_log2FC',
                        y = 'p_val_adj',
                        xlab = "log2 Fold Change",
                        ylab = "-log10 Adjusted p-value",
                        title = "Volcano Plot: Adenoma vs 1-year WT",
                        pCutoff = 0.05,
                        FCcutoff = 2,
                        pointSize = 2,
                        labSize = 4,
                        colAlpha = 0.4,
                        col = c("grey10", "grey10", "grey", "purple"),  # down, up, p-adj only, both
                        legendLabels = c('NS', '', 'p < 0.05, log2FC < 2', 'p < 0.05, log2FC > 2'),
                        legendPosition = 'bottom',
                        drawConnectors = FALSE, 
                        border = 'full')
print(vol1)
ggsave("volcano_plot_adeno1year.png", plot = vol1, width = 6, height = 6, units = "in", dpi = 1200)

markers_adeno_vs_2week <- markers_adeno_vs_2week %>%
  mutate(
    negLog10P = -log10(p_val_adj),
    negLog10P = ifelse(
      is.infinite(negLog10P),
      max(negLog10P[is.finite(negLog10P)], na.rm = TRUE) + 1,
      negLog10P
    ),
    regulation = case_when(
      avg_log2FC > 0 & p_val_adj < 0.05 ~ "Upregulated",
      avg_log2FC < 0 & p_val_adj < 0.05 ~ "Downregulated",
      TRUE ~ "Not significant"
    )
  )

vol2 <- EnhancedVolcano(markers_adeno_vs_2week,
                        lab = rownames(markers_adeno_vs_2week),
                        x = 'avg_log2FC',
                        y = 'p_val_adj',
                        xlab = "log2 Fold Change",
                        ylab = "-log10 Adjusted p-value",
                        title = "Volcano Plot: Adenoma vs 2-week WT",
                        pCutoff = 0.05,
                        FCcutoff = 2,
                        pointSize = 2,
                        labSize = 4,
                        colAlpha = 0.4,
                        col = c("grey10", "grey10", "grey", "purple"),  # down, up, p-adj only, both
                        legendLabels = c('NS', '', 'p < 0.05, log2FC < 2', 'p < 0.05, log2FC > 2'),
                        legendPosition = 'bottom',
                        drawConnectors = FALSE, 
                        border = 'full')
print(vol2)
ggsave("volcano_plot_adeno2week.png", plot = vol2, width = 6, height = 6, units = "in", dpi = 1200)

markers_1year_vs_2week <- markers_1year_vs_2week %>%
  mutate(
    negLog10P = -log10(p_val_adj),
    negLog10P = ifelse(
      is.infinite(negLog10P),
      max(negLog10P[is.finite(negLog10P)], na.rm = TRUE) + 1,
      negLog10P
    ),
    regulation = case_when(
      avg_log2FC > 0 & p_val_adj < 0.05 ~ "Upregulated",
      avg_log2FC < 0 & p_val_adj < 0.05 ~ "Downregulated",
      TRUE ~ "Not significant"
    )
  )

vol3 <- EnhancedVolcano(markers_1year_vs_2week,
                        lab = rownames(markers_1year_vs_2week),
                        x = 'avg_log2FC',
                        y = 'p_val_adj',
                        xlab = "log2 Fold Change",
                        ylab = "-log10 Adjusted p-value",
                        title = "Volcano Plot: 1-year WT vs 2-week WT",
                        pCutoff = 0.05,
                        FCcutoff = 2,
                        pointSize = 2,
                        labSize = 4,
                        colAlpha = 0.4,
                        col = c("grey10", "grey10", "grey", "purple"),  # down, up, p-adj only, both
                        legendLabels = c('NS', '', 'p < 0.05, log2FC < 2', 'p < 0.05, log2FC > 2'),
                        legendPosition = 'bottom',
                        drawConnectors = FALSE, 
                        border = 'full')
print(vol3)
ggsave("volcano_plot_1year2week.png", plot = vol3, width = 6, height = 6, units = "in", dpi = 1200)


##################GENE ONTOLOGY ANALYSES OF STROMA############################

# Install BiocManager if you don't have it
#if (!require("BiocManager", quietly = TRUE))
#install.packages("BiocManager")

# Install the packages
#BiocManager::install("clusterProfiler")
#BiocManager::install("org.Mm.eg.db")

# Load the libraries to test
library(clusterProfiler)
library(org.Mm.eg.db)

stroma_age_degs <- read.csv("/mnt/rosalindData/sarah/single-cell/seurat/outputs/age-degs.csv")

# Create a list of genes for each age group
age_gene_lists <- list(
  "2week" = subset_age_markers[subset_age_markers$cluster == "2week", "gene"],
  "1year" = subset_age_markers[subset_age_markers$cluster == "1year", "gene"],
  "1year_adeno" = subset_age_markers[subset_age_markers$cluster == "1year_adeno", "gene"]
)

# Convert all to Entrez IDs using lapply
age_entrez_lists <- lapply(age_gene_lists, function(genes) {
  bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
})

# Access individual results
gene_ids_2week <- age_entrez_lists$`2week`
gene_ids_1year <- age_entrez_lists$`1year`
gene_ids_1year_adeno <- age_entrez_lists$`1year_adeno`

# Run GO analysis for each age group BP
BP_2week <- enrichGO(gene = age_entrez_lists$`2week`$ENTREZID,
                     OrgDb = org.Mm.eg.db,
                     ont = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05)

BP_1year <- enrichGO(gene = age_entrez_lists$`1year`$ENTREZID,
                     OrgDb = org.Mm.eg.db,
                     ont = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05)

BP_1year_adeno <- enrichGO(gene = age_entrez_lists$`1year_adeno`$ENTREZID,
                           OrgDb = org.Mm.eg.db,
                           ont = "BP",
                           pAdjustMethod = "BH",
                           pvalueCutoff = 0.05)

# Extract and combine top 10 from each group
extract_top_terms <- function(enrichGO_result, group_name, n = 10) {
  if(nrow(enrichGO_result@result) == 0) {
    return(data.frame())
  }
  
  result_df <- enrichGO_result@result %>%
    arrange(p.adjust) %>%
    slice_head(n = n) %>%
    mutate(Group = group_name)
  
  return(result_df)
}

# Extract top 10 from each group
top_2week_bp <- extract_top_terms(BP_2week, "2week", 5)
top_1year_bp <- extract_top_terms(BP_1year, "1year", 5)
top_1year_adeno_bp <- extract_top_terms(BP_1year_adeno, "1year_adeno", 5)

# Combine all results
combined_results_bp <- bind_rows(top_2week_bp, top_1year_bp, top_1year_adeno_bp)

combined_results_bp <- combined_results_bp %>%
  mutate(Ontology = "BP")

# Run GO analysis for each age group CC
CC_2week <- enrichGO(gene = age_entrez_lists$`2week`$ENTREZID,
                     OrgDb = org.Mm.eg.db,
                     ont = "CC",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05)

CC_1year <- enrichGO(gene = age_entrez_lists$`1year`$ENTREZID,
                     OrgDb = org.Mm.eg.db,
                     ont = "CC",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05)

CC_1year_adeno <- enrichGO(gene = age_entrez_lists$`1year_adeno`$ENTREZID,
                           OrgDb = org.Mm.eg.db,
                           ont = "CC",
                           pAdjustMethod = "BH",
                           pvalueCutoff = 0.05)



# Extract top 10 from each group
top_2week_cc <- extract_top_terms(CC_2week, "2week", 5)
top_1year_cc <- extract_top_terms(CC_1year, "1year", 5)
top_1year_adeno_cc <- extract_top_terms(CC_1year_adeno, "1year_adeno", 5)

# Combine all results
combined_results_cc <- bind_rows(top_2week_cc, top_1year_cc, top_1year_adeno_cc)

combined_results_cc <- combined_results_cc %>%
  mutate(Ontology = "CC")

############MF
# Run GO analysis for each age group MF
MF_2week <- enrichGO(gene = age_entrez_lists$`2week`$ENTREZID,
                     OrgDb = org.Mm.eg.db,
                     ont = "MF",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05)

MF_1year <- enrichGO(gene = age_entrez_lists$`1year`$ENTREZID,
                     OrgDb = org.Mm.eg.db,
                     ont = "MF",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05)

MF_1year_adeno <- enrichGO(gene = age_entrez_lists$`1year_adeno`$ENTREZID,
                           OrgDb = org.Mm.eg.db,
                           ont = "MF",
                           pAdjustMethod = "BH",
                           pvalueCutoff = 0.05)


# Extract top 10 from each group
top_2week_mf <- extract_top_terms(MF_2week, "2week", 5)
top_1year_mf <- extract_top_terms(MF_1year, "1year", 5)
top_1year_adeno_mf <- extract_top_terms(MF_1year_adeno, "1year_adeno", 5)

# Combine all results
combined_results_mf <- bind_rows(top_2week_mf, top_1year_mf, top_1year_adeno_mf)
combined_results_mf <- combined_results_mf %>%
  mutate(Ontology = "MF")

#####Combine all terms
combined_results_all <- bind_rows(combined_results_bp, combined_results_cc, combined_results_mf)

levels_age <- c("2week", "1year", "1year_adeno")
levels_GO <- c("BP", "CC", "MF")
combined_results_all$Group <- factor(combined_results_all$Group, levels = levels_age)
combined_results_all$Ontology <- factor(combined_results_all$Ontology, levels = levels_GO)

combined_results_all <- combined_results_all %>%
  arrange(Ontology, p.adjust) %>%
  mutate(Description = factor(Description, levels = rev(unique(Description))))

#plot
go <- ggplot(combined_results_all, aes(x = Ontology, y = Description)) +
  geom_point(aes(size = Count, color = -log10(p.adjust))) +
  scale_color_gradient(low = "blue", high = "red", name = "-log10(p.adj)") +
  scale_size_continuous(name = "Count", range = c(2, 8)) +
  facet_wrap(~Group, scales = "fixed", ncol = 3) +
  labs(
    x = "Ontology",
    y = "GO Terms",
    title = "Top 5 GO Terms"
  ) +
  theme_grey() +
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
    strip.text = element_text(size = 10, face = "bold")
  )
print(go)
ggsave("go-stroma-age-top5-reverse.png", plot = go, width = 9, height = 10, dpi = 800)

##################KEGG PATHWAYS###############################

# Run GO analysis for each age group KEGG Pathway Analysis
KEGG_2week <- enrichKEGG(gene = age_entrez_lists$`2week`$ENTREZID,
                         organism = "mmu",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.05)

KEGG_1year <- enrichKEGG(gene = age_entrez_lists$`1year`$ENTREZID,
                         organism = "mmu",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.05)

KEGG_1year_adeno <- enrichKEGG(gene = age_entrez_lists$`1year_adeno`$ENTREZID,
                               organism = "mmu",
                               pAdjustMethod = "BH",
                               pvalueCutoff = 0.05)

# Extract top 10 from each group
top_2week_KEGG <- extract_top_terms(KEGG_2week, "2week", 10)
top_1year_KEGG <- extract_top_terms(KEGG_1year, "1year", 10)
top_1year_adeno_KEGG <- extract_top_terms(KEGG_1year_adeno, "1year_adeno", 10)

# Combine all results
combined_results_KEGG <- bind_rows(top_2week_KEGG, top_1year_KEGG, top_1year_adeno_KEGG)

levels_age <- c("2week", "1year", "1year_adeno")
combined_results_KEGG$Group <- factor(combined_results_KEGG$Group, levels = levels_age)


KEGG <- ggplot(combined_results_KEGG, aes(x = Group, y = Description)) +
  geom_point(aes(size = Count, color = -log10(p.adjust))) +
  scale_color_gradient(low = "blue", high = "red", name = "-log10(p.adj)") +
  scale_size_continuous(name = "Count", range = c(2, 8)) +
  #facet_wrap(~subcategory, scales = "free_y", nrow = 1) +
  labs(
    x = "Age",
    y = "KEGG Pathway",
    title = "Top 10 KEGG Terms"
  ) +
  theme_grey() +
  theme(
    axis.text.y = element_text(size = 12),  # Smaller text to fit more pathways
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
    strip.text = element_text(size = 12, face = "bold"),
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.spacing = unit(0.5, "lines"),  # Adjust spacing between facets
    legend.position = "right"
  )

print(KEGG)

ggsave("stroma-KEGG-top5.png", plot = KEGG, width = 8, height = 8, dpi = 800)
##################GENE ONTOLOGY ANALYSES OF STROMA - FROM 10X DEGS############################

# Install BiocManager if you don't have it
#if (!require("BiocManager", quietly = TRUE))
#install.packages("BiocManager")

# Install the packages
#BiocManager::install("clusterProfiler")
#BiocManager::install("org.Mm.eg.db")

# Load the libraries to test
library(clusterProfiler)
library(org.Mm.eg.db)
library(dplyr)
library(tidyverse)



#manually generate DEGS from loupe browser and bring them back into R
degs_2week <- read.csv("/mnt/rosalindData/sarah/single-cell/seurat/outputs/2week-degs.csv")
degs_1year <- read.csv("/mnt/rosalindData/sarah/single-cell/seurat/outputs/1year-degs.csv")
degs_1year_adeno <- read.csv("/mnt/rosalindData/sarah/single-cell/seurat/outputs/1year_adeno-degs.csv")

#subset/filter for each age - GO
genes_2week <- degs_2week %>%
  filter(X2week.log2FC >= 0.25, X2week.p.val < 0.05) %>%
  pull(FeatureName)

genes_1year <- degs_1year %>%
  filter(X1year.log2FC >= 0.25, X1year.p.val < 0.05) %>%
  pull(FeatureName)

genes_1year_adeno <- degs_1year_adeno %>%
  filter(X1year_adeno.log2FC >= 0.25, X1year_adeno.p.val < 0.05) %>%
  pull(FeatureName)


# Convert all to Entrez IDs for GO
genes_2week <- bitr(genes_2week, 
                    fromType = "SYMBOL", 
                    toType = "ENTREZID", 
                    OrgDb = org.Mm.eg.db)

genes_1year <- bitr(genes_1year, 
                    fromType = "SYMBOL", 
                    toType = "ENTREZID", 
                    OrgDb = org.Mm.eg.db)

genes_1year_adeno <- bitr(genes_1year_adeno, 
                          fromType = "SYMBOL", 
                          toType = "ENTREZID", 
                          OrgDb = org.Mm.eg.db)



# Run GO analysis for each age group BP
BP_2week <- enrichGO(gene = genes_2week$ENTREZID,
                     OrgDb = org.Mm.eg.db,
                     ont = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05)

BP_1year <- enrichGO(gene = genes_1year$ENTREZID,
                     OrgDb = org.Mm.eg.db,
                     ont = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05)

BP_1year_adeno <- enrichGO(gene = genes_1year_adeno$ENTREZID,
                           OrgDb = org.Mm.eg.db,
                           ont = "BP",
                           pAdjustMethod = "BH",
                           pvalueCutoff = 0.05)

# Extract and combine top 10 from each group
extract_top_terms <- function(enrichGO_result, group_name, n = 10) {
  if(nrow(enrichGO_result@result) == 0) {
    return(data.frame())
  }
  
  result_df <- enrichGO_result@result %>%
    arrange(p.adjust) %>%
    slice_head(n = n) %>%
    mutate(Group = group_name)
  
  return(result_df)
}

# Extract top 10 from each group
top_2week_bp <- extract_top_terms(BP_2week, "2week", 5)
top_1year_bp <- extract_top_terms(BP_1year, "1year", 5)
top_1year_adeno_bp <- extract_top_terms(BP_1year_adeno, "1year_adeno", 5)

# Combine all results
combined_results_bp <- bind_rows(top_2week_bp, top_1year_bp, top_1year_adeno_bp)

combined_results_bp <- combined_results_bp %>%
  mutate(Ontology = "BP")

# Run GO analysis for each age group CC
CC_2week <- enrichGO(gene = genes_2week$ENTREZID,
                     OrgDb = org.Mm.eg.db,
                     ont = "CC",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05)

CC_1year <- enrichGO(gene = genes_1year$ENTREZID,
                     OrgDb = org.Mm.eg.db,
                     ont = "CC",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05)

CC_1year_adeno <- enrichGO(gene = genes_1year_adeno$ENTREZID,
                           OrgDb = org.Mm.eg.db,
                           ont = "CC",
                           pAdjustMethod = "BH",
                           pvalueCutoff = 0.05)



# Extract top 10 from each group
top_2week_cc <- extract_top_terms(CC_2week, "2week", 5)
top_1year_cc <- extract_top_terms(CC_1year, "1year", 5)
top_1year_adeno_cc <- extract_top_terms(CC_1year_adeno, "1year_adeno", 5)

# Combine all results
combined_results_cc <- bind_rows(top_2week_cc, top_1year_cc, top_1year_adeno_cc)

combined_results_cc <- combined_results_cc %>%
  mutate(Ontology = "CC")

############MF
# Run GO analysis for each age group MF
MF_2week <- enrichGO(gene = genes_2week$ENTREZID,
                     OrgDb = org.Mm.eg.db,
                     ont = "MF",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05)

MF_1year <- enrichGO(gene = genes_1year$ENTREZID,
                     OrgDb = org.Mm.eg.db,
                     ont = "MF",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05)

MF_1year_adeno <- enrichGO(gene = genes_1year_adeno$ENTREZID,
                           OrgDb = org.Mm.eg.db,
                           ont = "MF",
                           pAdjustMethod = "BH",
                           pvalueCutoff = 0.05)


# Extract top 10 from each group
top_2week_mf <- extract_top_terms(MF_2week, "2week", 5)
top_1year_mf <- extract_top_terms(MF_1year, "1year", 5)
top_1year_adeno_mf <- extract_top_terms(MF_1year_adeno, "1year_adeno", 5)

# Combine all results
combined_results_mf <- bind_rows(top_2week_mf, top_1year_mf, top_1year_adeno_mf)
combined_results_mf <- combined_results_mf %>%
  mutate(Ontology = "MF")

#####Combine all terms
combined_results_all <- bind_rows(combined_results_bp, combined_results_cc, combined_results_mf)

levels_age <- c("2week", "1year", "1year_adeno")
levels_GO <- c("BP", "CC", "MF")
combined_results_all$Group <- factor(combined_results_all$Group, levels = levels_age)
combined_results_all$Ontology <- factor(combined_results_all$Ontology, levels = levels_GO)

combined_results_all <- combined_results_all %>%
  arrange(Ontology, p.adjust) %>%
  mutate(Description = factor(Description, levels = rev(unique(Description))))

#plot
go <- ggplot(combined_results_all, aes(x = Ontology, y = Description)) +
  geom_point(aes(size = Count, color = -log10(p.adjust))) +
  scale_color_gradient(low = "blue", high = "red", name = "-log10(p.adj)") +
  scale_size_continuous(name = "Count", range = c(2, 8)) +
  facet_wrap(~Group, scales = "fixed", ncol = 3) +
  labs(
    x = "Ontology",
    y = "GO Terms",
    title = "Top 5 GO Terms"
  ) +
  theme_grey() +
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
    strip.text = element_text(size = 10, face = "bold")
  )
print(go)
ggsave("go-stroma-age-top5-degs-rev.png", plot = go, width = 9, height = 10, dpi = 800)

##################KEGG PATHWAYS###############################

# Run GO analysis for each age group KEGG Pathway Analysis
KEGG_2week <- enrichKEGG(gene = genes_2week$ENTREZID,
                         organism = "mmu",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.05)

KEGG_1year <- enrichKEGG(gene = genes_1year$ENTREZID,
                         organism = "mmu",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.05)

KEGG_1year_adeno <- enrichKEGG(gene = genes_1year_adeno$ENTREZID,
                               organism = "mmu",
                               pAdjustMethod = "BH",
                               pvalueCutoff = 0.05)

# Extract top 10 from each group
top_2week_KEGG <- extract_top_terms(KEGG_2week, "2week", 10)
top_1year_KEGG <- extract_top_terms(KEGG_1year, "1year", 10)
top_1year_adeno_KEGG <- extract_top_terms(KEGG_1year_adeno, "1year_adeno", 10)

# Combine all results
combined_results_KEGG <- bind_rows(top_2week_KEGG, top_1year_KEGG, top_1year_adeno_KEGG)

levels_age <- c("2week", "1year", "1year_adeno")
combined_results_KEGG$Group <- factor(combined_results_KEGG$Group, levels = levels_age)


KEGG <- ggplot(combined_results_KEGG, aes(x = Group, y = Description)) +
  geom_point(aes(size = Count, color = -log10(p.adjust))) +
  scale_color_gradient(low = "blue", high = "red", name = "-log10(p.adj)") +
  scale_size_continuous(name = "Count", range = c(2, 8)) +
  #facet_wrap(~subcategory, scales = "free_y", nrow = 1) +
  labs(
    x = "Age",
    y = "KEGG Pathway",
    title = "Top 10 KEGG Terms"
  ) +
  theme_grey() +
  theme(
    axis.text.y = element_text(size = 12),  # Smaller text to fit more pathways
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
    strip.text = element_text(size = 12, face = "bold"),
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.spacing = unit(0.5, "lines"),  # Adjust spacing between facets
    legend.position = "right"
  )

print(KEGG)

ggsave("stroma-KEGG-top5-degs.png", plot = KEGG, width = 8, height = 8, dpi = 800)

##################GENE ONTOLOGY ANALYSES OF STROMA - FROM FINDMARKERS PAIRWISE###########################

# Install BiocManager if you don't have it
#if (!require("BiocManager", quietly = TRUE))
#install.packages("BiocManager")

# Install the packages
#BiocManager::install("clusterProfiler")
#BiocManager::install("org.Mm.eg.db")

# Load the libraries to test
library(clusterProfiler)
library(org.Mm.eg.db)
library(dplyr)
library(tidyverse)

#subset/filter for each age - GO
genes_adeno_v_1year <- markers_adeno_vs_1year %>%
  tibble::rownames_to_column(var = "gene") %>%
  filter(avg_log2FC >= 0.25, p_val_adj < 0.05) %>%
  pull(gene)

genes_adeno_v_2week <- markers_adeno_vs_2week %>%
  tibble::rownames_to_column(var = "gene") %>%
  filter(avg_log2FC >= 0.25, p_val_adj < 0.05) %>%
  pull(gene)

genes_1year_v_2week <- markers_1year_vs_2week %>%
  tibble::rownames_to_column(var = "gene") %>%
  filter(avg_log2FC >= 0.25, p_val_adj < 0.05) %>%
  pull(gene)


# Convert all to Entrez IDs for GO
genes_adeno_v_1year <- bitr(genes_adeno_v_1year, 
                            fromType = "SYMBOL", 
                            toType = "ENTREZID", 
                            OrgDb = org.Mm.eg.db)

genes_adeno_v_2week <- bitr(genes_adeno_v_2week, 
                            fromType = "SYMBOL", 
                            toType = "ENTREZID", 
                            OrgDb = org.Mm.eg.db)

genes_1year_v_2week <- bitr(genes_1year_v_2week, 
                            fromType = "SYMBOL", 
                            toType = "ENTREZID", 
                            OrgDb = org.Mm.eg.db)



# Run GO analysis for each age group BP
BP_adeno_1year <- enrichGO(gene = genes_adeno_v_1year$ENTREZID,
                           OrgDb = org.Mm.eg.db,
                           ont = "BP",
                           pAdjustMethod = "BH",
                           pvalueCutoff = 0.05)

BP_adeno_2week <- enrichGO(gene = genes_adeno_v_2week$ENTREZID,
                           OrgDb = org.Mm.eg.db,
                           ont = "BP",
                           pAdjustMethod = "BH",
                           pvalueCutoff = 0.05)

BP_1year_2week <- enrichGO(gene = genes_1year_v_2week$ENTREZID,
                           OrgDb = org.Mm.eg.db,
                           ont = "BP",
                           pAdjustMethod = "BH",
                           pvalueCutoff = 0.05)

# Extract and combine top 10 from each group
extract_top_terms <- function(enrichGO_result, group_name, n = 10) {
  if(nrow(enrichGO_result@result) == 0) {
    return(data.frame())
  }
  
  result_df <- enrichGO_result@result %>%
    arrange(p.adjust) %>%
    slice_head(n = n) %>%
    mutate(Group = group_name)
  
  return(result_df)
}

# Extract top 10 from each group
top_adeno_1year_bp <- extract_top_terms(BP_adeno_1year, "Adeno/1Year", 5)
top_adeno_2week_bp <- extract_top_terms(BP_adeno_2week, "Adeno/2week", 5)
top_1year_2week_bp <- extract_top_terms(BP_1year_2week, "1year/2week", 5)

# Combine all results
combined_results_bp <- bind_rows(top_adeno_1year_bp, top_adeno_2week_bp, top_1year_2week_bp)

combined_results_bp <- combined_results_bp %>%
  mutate(Ontology = "BP")

##########CC##########
CC_adeno_1year <- enrichGO(gene = genes_adeno_v_1year$ENTREZID,
                           OrgDb = org.Mm.eg.db,
                           ont = "CC",
                           pAdjustMethod = "BH",
                           pvalueCutoff = 0.05)

CC_adeno_2week <- enrichGO(gene = genes_adeno_v_2week$ENTREZID,
                           OrgDb = org.Mm.eg.db,
                           ont = "CC",
                           pAdjustMethod = "BH",
                           pvalueCutoff = 0.05)

CC_1year_2week <- enrichGO(gene = genes_1year_v_2week$ENTREZID,
                           OrgDb = org.Mm.eg.db,
                           ont = "CC",
                           pAdjustMethod = "BH",
                           pvalueCutoff = 0.05)

# Extract top 10 from each group
top_adeno_1year_cc <- extract_top_terms(CC_adeno_1year, "Adeno/1Year", 5)
top_adeno_2week_cc <- extract_top_terms(CC_adeno_2week, "Adeno/2week", 5)
top_1year_2week_cc <- extract_top_terms(CC_1year_2week, "1year/2week", 5)

# Combine all results
combined_results_cc <- bind_rows(top_adeno_1year_cc, top_adeno_2week_cc, top_1year_2week_cc)

combined_results_cc <- combined_results_cc %>%
  mutate(Ontology = "CC")

############MF
# Run GO analysis for each age group MF
MF_adeno_1year <- enrichGO(gene = genes_adeno_v_1year$ENTREZID,
                           OrgDb = org.Mm.eg.db,
                           ont = "MF",
                           pAdjustMethod = "BH",
                           pvalueCutoff = 0.05)

MF_adeno_2week <- enrichGO(gene = genes_adeno_v_2week$ENTREZID,
                           OrgDb = org.Mm.eg.db,
                           ont = "MF",
                           pAdjustMethod = "BH",
                           pvalueCutoff = 0.05)

MF_1year_2week <- enrichGO(gene = genes_1year_v_2week$ENTREZID,
                           OrgDb = org.Mm.eg.db,
                           ont = "MF",
                           pAdjustMethod = "BH",
                           pvalueCutoff = 0.05)

# Extract top 10 from each group
top_adeno_1year_mf <- extract_top_terms(MF_adeno_1year, "Adeno/1Year", 5)
top_adeno_2week_mf <- extract_top_terms(MF_adeno_2week, "Adeno/2week", 5)
top_1year_2week_mf <- extract_top_terms(MF_1year_2week, "1year/2week", 5)

# Combine all results
combined_results_mf <- bind_rows(top_adeno_1year_mf, top_adeno_2week_mf, top_1year_2week_mf)

combined_results_mf <- combined_results_mf %>%
  mutate(Ontology = "MF")


#####Combine all terms
combined_results_all <- bind_rows(combined_results_bp, combined_results_cc, combined_results_mf)
levels_age <- c("Adeno/1Year", "Adeno/2week", "1year/2week")
levels_GO <- c("BP", "CC", "MF")
combined_results_all$Group <- factor(combined_results_all$Group, levels = levels_age)
combined_results_all$Ontology <- factor(combined_results_all$Ontology, levels = levels_GO)

combined_results_all <- combined_results_all %>%
  arrange(Ontology, p.adjust) %>%
  mutate(Description = factor(Description, levels = rev(unique(Description))))

#plot
go <- ggplot(combined_results_all, aes(x = Ontology, y = Description)) +
  geom_point(aes(size = Count, color = -log10(p.adjust))) +
  scale_color_gradient(low = "blue", high = "red", name = "-log10(p.adj)") +
  scale_size_continuous(name = "Count", range = c(2, 8)) +
  facet_wrap(~Group, scales = "fixed", ncol = 3) +
  labs(
    x = "Ontology",
    y = "GO Terms",
    title = "Top 5 GO Terms"
  ) +
  theme_grey() +
  theme(
    axis.text.y = element_text(size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 16),
    strip.text = element_text(size = 10, face = "bold")
  )
print(go)
ggsave("go-stroma-age-top5-degs-rev-findmarkers.png", plot = go, width = 12, height = 10, dpi = 800)

##################KEGG PATHWAYS FINDMARKERS PAIRWISE##############################
# Run GO analysis for each age group KEGG Pathway Analysis
KEGG_adeno_v_1year <- enrichKEGG(gene = genes_adeno_v_1year$ENTREZID,
                                 organism = "mmu",
                                 pAdjustMethod = "BH",
                                 pvalueCutoff = 0.05)

KEGG_adeno_v_2week <- enrichKEGG(gene = genes_adeno_v_2week$ENTREZID,
                                 organism = "mmu",
                                 pAdjustMethod = "BH",
                                 pvalueCutoff = 0.05)

KEGG_1year_v_2week <- enrichKEGG(gene = genes_1year_v_2week$ENTREZID,
                                 organism = "mmu",
                                 pAdjustMethod = "BH",
                                 pvalueCutoff = 0.05)

# Extract top 10 from each group
top_adeno_1year_KEGG <- extract_top_terms(KEGG_adeno_v_1year, "Adeno/1Year", 10)
top_adeno_2week_KEGG <- extract_top_terms(KEGG_adeno_v_2week, "Adeno/2week", 10)
top_1year_2week_KEGG <- extract_top_terms(KEGG_1year_v_2week, "1year/2week", 10)

# Combine all results
combined_results_KEGG <- bind_rows(top_adeno_1year_KEGG, top_adeno_2week_KEGG, top_1year_2week_KEGG)

levels_age <- c("Adeno/1Year", "Adeno/2week", "1year/2week")
combined_results_KEGG$Group <- factor(combined_results_KEGG$Group, levels = levels_age)


KEGG <- ggplot(combined_results_KEGG, aes(x = Group, y = Description)) +
  geom_point(aes(size = Count, color = -log10(p.adjust))) +
  scale_color_gradient(low = "blue", high = "red", name = "-log10(p.adj)") +
  scale_size_continuous(name = "Count", range = c(2, 8)) +
  #facet_wrap(~subcategory, scales = "free_y", nrow = 1) +
  labs(
    x = "Age",
    y = "KEGG Pathway",
    title = "Top 10 KEGG Terms"
  ) +
  theme_grey() +
  theme(
    axis.text.y = element_text(size = 12),  # Smaller text to fit more pathways
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
    strip.text = element_text(size = 12, face = "bold"),
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.spacing = unit(0.5, "lines"),  # Adjust spacing between facets
    legend.position = "right"
  )

print(KEGG)

ggsave("stroma-KEGG-top5-degs-pairwise.png", plot = KEGG, width = 8, height = 8, dpi = 1200)

##################GENE ONTOLOGY ANALYSES OF STROMA - FROM FINDMARKERS PAIRWISE###########################
#################################DOWNREGULATED##########################################################
#subset/filter for each age - GO
genes_adeno_v_1year <- markers_adeno_vs_1year %>%
  tibble::rownames_to_column(var = "gene") %>%
  filter(avg_log2FC <= -0.25, p_val_adj < 0.05) %>%
  pull(gene)

genes_adeno_v_2week <- markers_adeno_vs_2week %>%
  tibble::rownames_to_column(var = "gene") %>%
  filter(avg_log2FC <= -0.25, p_val_adj < 0.05) %>%
  pull(gene)

genes_1year_v_2week <- markers_1year_vs_2week %>%
  tibble::rownames_to_column(var = "gene") %>%
  filter(avg_log2FC <= -0.25, p_val_adj < 0.05) %>%
  pull(gene)


# Convert all to Entrez IDs for GO
genes_adeno_v_1year <- bitr(genes_adeno_v_1year, 
                            fromType = "SYMBOL", 
                            toType = "ENTREZID", 
                            OrgDb = org.Mm.eg.db)

genes_adeno_v_2week <- bitr(genes_adeno_v_2week, 
                            fromType = "SYMBOL", 
                            toType = "ENTREZID", 
                            OrgDb = org.Mm.eg.db)

genes_1year_v_2week <- bitr(genes_1year_v_2week, 
                            fromType = "SYMBOL", 
                            toType = "ENTREZID", 
                            OrgDb = org.Mm.eg.db)



# Run GO analysis for each age group BP
BP_adeno_1year <- enrichGO(gene = genes_adeno_v_1year$ENTREZID,
                           OrgDb = org.Mm.eg.db,
                           ont = "BP",
                           pAdjustMethod = "BH",
                           pvalueCutoff = 0.05)

BP_adeno_2week <- enrichGO(gene = genes_adeno_v_2week$ENTREZID,
                           OrgDb = org.Mm.eg.db,
                           ont = "BP",
                           pAdjustMethod = "BH",
                           pvalueCutoff = 0.05)

BP_1year_2week <- enrichGO(gene = genes_1year_v_2week$ENTREZID,
                           OrgDb = org.Mm.eg.db,
                           ont = "BP",
                           pAdjustMethod = "BH",
                           pvalueCutoff = 0.05)

# Extract and combine top 10 from each group
extract_top_terms <- function(enrichGO_result, group_name, n = 10) {
  if(nrow(enrichGO_result@result) == 0) {
    return(data.frame())
  }
  
  result_df <- enrichGO_result@result %>%
    arrange(p.adjust) %>%
    slice_head(n = n) %>%
    mutate(Group = group_name)
  
  return(result_df)
}

# Extract top 10 from each group
top_adeno_1year_bp <- extract_top_terms(BP_adeno_1year, "Adeno/1Year", 5)
top_adeno_2week_bp <- extract_top_terms(BP_adeno_2week, "Adeno/2week", 5)
top_1year_2week_bp <- extract_top_terms(BP_1year_2week, "1year/2week", 5)

# Combine all results
combined_results_bp <- bind_rows(top_adeno_1year_bp, top_adeno_2week_bp, top_1year_2week_bp)

combined_results_bp <- combined_results_bp %>%
  mutate(Ontology = "BP")

##########CC##########
CC_adeno_1year <- enrichGO(gene = genes_adeno_v_1year$ENTREZID,
                           OrgDb = org.Mm.eg.db,
                           ont = "CC",
                           pAdjustMethod = "BH",
                           pvalueCutoff = 0.05)

CC_adeno_2week <- enrichGO(gene = genes_adeno_v_2week$ENTREZID,
                           OrgDb = org.Mm.eg.db,
                           ont = "CC",
                           pAdjustMethod = "BH",
                           pvalueCutoff = 0.05)

CC_1year_2week <- enrichGO(gene = genes_1year_v_2week$ENTREZID,
                           OrgDb = org.Mm.eg.db,
                           ont = "CC",
                           pAdjustMethod = "BH",
                           pvalueCutoff = 0.05)

# Extract top 10 from each group
top_adeno_1year_cc <- extract_top_terms(CC_adeno_1year, "Adeno/1Year", 5)
top_adeno_2week_cc <- extract_top_terms(CC_adeno_2week, "Adeno/2week", 5)
top_1year_2week_cc <- extract_top_terms(CC_1year_2week, "1year/2week", 5)

# Combine all results
combined_results_cc <- bind_rows(top_adeno_1year_cc, top_adeno_2week_cc, top_1year_2week_cc)

combined_results_cc <- combined_results_cc %>%
  mutate(Ontology = "CC")

############MF
# Run GO analysis for each age group MF
MF_adeno_1year <- enrichGO(gene = genes_adeno_v_1year$ENTREZID,
                           OrgDb = org.Mm.eg.db,
                           ont = "MF",
                           pAdjustMethod = "BH",
                           pvalueCutoff = 0.05)

MF_adeno_2week <- enrichGO(gene = genes_adeno_v_2week$ENTREZID,
                           OrgDb = org.Mm.eg.db,
                           ont = "MF",
                           pAdjustMethod = "BH",
                           pvalueCutoff = 0.05)

MF_1year_2week <- enrichGO(gene = genes_1year_v_2week$ENTREZID,
                           OrgDb = org.Mm.eg.db,
                           ont = "MF",
                           pAdjustMethod = "BH",
                           pvalueCutoff = 0.05)

# Extract top 10 from each group
top_adeno_1year_mf <- extract_top_terms(MF_adeno_1year, "Adeno/1Year", 5)
top_adeno_2week_mf <- extract_top_terms(MF_adeno_2week, "Adeno/2week", 5)
top_1year_2week_mf <- extract_top_terms(MF_1year_2week, "1year/2week", 5)

# Combine all results
combined_results_mf <- bind_rows(top_adeno_1year_mf, top_adeno_2week_mf, top_1year_2week_mf)

combined_results_mf <- combined_results_mf %>%
  mutate(Ontology = "MF")


#####Combine all terms
combined_results_all <- bind_rows(combined_results_bp, combined_results_cc, combined_results_mf)
levels_age <- c("Adeno/1Year", "Adeno/2week", "1year/2week")
levels_GO <- c("BP", "CC", "MF")
combined_results_all$Group <- factor(combined_results_all$Group, levels = levels_age)
combined_results_all$Ontology <- factor(combined_results_all$Ontology, levels = levels_GO)

combined_results_all <- combined_results_all %>%
  arrange(Ontology, p.adjust) %>%
  mutate(Description = factor(Description, levels = rev(unique(Description))))

#plot
go <- ggplot(combined_results_all, aes(x = Ontology, y = Description)) +
  geom_point(aes(size = Count, color = -log10(p.adjust))) +
  scale_color_gradient(low = "blue", high = "red", name = "-log10(p.adj)") +
  scale_size_continuous(name = "Count", range = c(2, 8)) +
  facet_wrap(~Group, scales = "fixed", ncol = 3) +
  labs(
    x = "Ontology",
    y = "GO Terms",
    title = "Top 5 GO Terms"
  ) +
  theme_grey() +
  theme(
    axis.text.y = element_text(size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 16),
    strip.text = element_text(size = 10, face = "bold")
  )
print(go)
ggsave("go-stroma-age-top5-degs-rev-findmarkers-downregulated.png", plot = go, width = 14, height = 10, dpi = 800)

##################KEGG PATHWAYS FINDMARKERS PAIRWISE##############################
# Run GO analysis for each age group KEGG Pathway Analysis
KEGG_adeno_v_1year <- enrichKEGG(gene = genes_adeno_v_1year$ENTREZID,
                                 organism = "mmu",
                                 pAdjustMethod = "BH",
                                 pvalueCutoff = 0.05)

KEGG_adeno_v_2week <- enrichKEGG(gene = genes_adeno_v_2week$ENTREZID,
                                 organism = "mmu",
                                 pAdjustMethod = "BH",
                                 pvalueCutoff = 0.05)

KEGG_1year_v_2week <- enrichKEGG(gene = genes_1year_v_2week$ENTREZID,
                                 organism = "mmu",
                                 pAdjustMethod = "BH",
                                 pvalueCutoff = 0.05)

# Extract top 10 from each group
top_adeno_1year_KEGG <- extract_top_terms(KEGG_adeno_v_1year, "Adeno/1Year", 10)
top_adeno_2week_KEGG <- extract_top_terms(KEGG_adeno_v_2week, "Adeno/2week", 10)
top_1year_2week_KEGG <- extract_top_terms(KEGG_1year_v_2week, "1year/2week", 10)

# Combine all results
combined_results_KEGG <- bind_rows(top_adeno_1year_KEGG, top_adeno_2week_KEGG, top_1year_2week_KEGG)

levels_age <- c("Adeno/1Year", "Adeno/2week", "1year/2week")
combined_results_KEGG$Group <- factor(combined_results_KEGG$Group, levels = levels_age)


KEGG <- ggplot(combined_results_KEGG, aes(x = Group, y = Description)) +
  geom_point(aes(size = Count, color = -log10(p.adjust))) +
  scale_color_gradient(low = "blue", high = "red", name = "-log10(p.adj)") +
  scale_size_continuous(name = "Count", range = c(2, 8)) +
  #facet_wrap(~subcategory, scales = "free_y", nrow = 1) +
  labs(
    x = "Age",
    y = "KEGG Pathway",
    title = "Top 10 KEGG Terms"
  ) +
  theme_grey() +
  theme(
    axis.text.y = element_text(size = 12),  # Smaller text to fit more pathways
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
    strip.text = element_text(size = 12, face = "bold"),
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.spacing = unit(0.5, "lines"),  # Adjust spacing between facets
    legend.position = "right"
  )

print(KEGG)

ggsave("stroma-KEGG-top5-degs-pairwise-downregulated.png", plot = KEGG, width = 8, height = 8, dpi = 1200)


##################KEGG PATHWAYS###############################

# Run GO analysis for each age group KEGG Pathway Analysis
KEGG_2week <- enrichKEGG(gene = genes_2week$ENTREZID,
                         organism = "mmu",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.05)

KEGG_1year <- enrichKEGG(gene = genes_1year$ENTREZID,
                         organism = "mmu",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.05)

KEGG_1year_adeno <- enrichKEGG(gene = genes_1year_adeno$ENTREZID,
                               organism = "mmu",
                               pAdjustMethod = "BH",
                               pvalueCutoff = 0.05)

# Extract top 10 from each group
top_2week_KEGG <- extract_top_terms(KEGG_2week, "2week", 10)
top_1year_KEGG <- extract_top_terms(KEGG_1year, "1year", 10)
top_1year_adeno_KEGG <- extract_top_terms(KEGG_1year_adeno, "1year_adeno", 10)

# Combine all results
combined_results_KEGG <- bind_rows(top_2week_KEGG, top_1year_KEGG, top_1year_adeno_KEGG)

levels_age <- c("2week", "1year", "1year_adeno")
combined_results_KEGG$Group <- factor(combined_results_KEGG$Group, levels = levels_age)


KEGG <- ggplot(combined_results_KEGG, aes(x = Group, y = Description)) +
  geom_point(aes(size = Count, color = -log10(p.adjust))) +
  scale_color_gradient(low = "blue", high = "red", name = "-log10(p.adj)") +
  scale_size_continuous(name = "Count", range = c(2, 8)) +
  #facet_wrap(~subcategory, scales = "free_y", nrow = 1) +
  labs(
    x = "Age",
    y = "KEGG Pathway",
    title = "Top 10 KEGG Terms"
  ) +
  theme_grey() +
  theme(
    axis.text.y = element_text(size = 12),  # Smaller text to fit more pathways
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
    strip.text = element_text(size = 12, face = "bold"),
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.spacing = unit(0.5, "lines"),  # Adjust spacing between facets
    legend.position = "right"
  )

print(KEGG)

ggsave("stroma-KEGG-top5-degs.png", plot = KEGG, width = 8, height = 8, dpi = 800)

#############################GSEA of Stroma Subset#########################
#subset/filter for each age - GSEA
gene_ranks_2week <- degs_2week %>%
  dplyr::filter(
    !is.na(X2week.log2FC),
    X2week.log2FC >= 0.25,
    X2week.p.val < 0.05
  ) %>%
  dplyr::select(FeatureName, X2week.log2FC) %>%
  tibble::deframe()

gene_ranks_2week <- sort(gene_ranks_2week, decreasing = TRUE)

gene_ranks_1year <- degs_1year %>%
  dplyr::filter(
    !is.na(X1year.log2FC),
    X1year.log2FC >= 0.25,
    X1year.p.val < 0.05
  ) %>%
  dplyr::select(FeatureName, X1year.log2FC) %>%
  tibble::deframe()

gene_ranks_1year <- sort(gene_ranks_1year, decreasing = TRUE)

gene_ranks_1year_adeno <- degs_1year_adeno %>%
  dplyr::filter(
    !is.na(X1year_adeno.log2FC),
    X1year_adeno.log2FC >= 0.25,
    X1year_adeno.p.val < 0.05
  ) %>%
  dplyr::select(FeatureName, X1year_adeno.log2FC) %>%
  tibble::deframe()

gene_ranks_1year_adeno <- sort(gene_ranks_1year_adeno, decreasing = TRUE)

##function to convert named numeric vector into entrez ID while maintaining ranking
convert_symbols_to_entrez <- function(gene_rank_vector, OrgDb = org.Mm.eg.db) {
  # Use bitr to convert gene symbols to Entrez IDs
  conversion <- bitr(names(gene_rank_vector),
                     fromType = "SYMBOL",
                     toType = "ENTREZID",
                     OrgDb = OrgDb)
  
  # Merge logFC (values) with converted IDs
  df <- data.frame(SYMBOL = names(gene_rank_vector),
                   logFC = gene_rank_vector)
  
  merged <- merge(conversion, df, by = "SYMBOL")
  
  # Build final vector: values are logFC, names are ENTREZID
  gene_list <- merged$logFC
  names(gene_list) <- merged$ENTREZID
  
  # Sort for GSEA
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  return(gene_list)
}

#this stores ranked entrezids in a numeric vector

gene_ranks_2week <- convert_symbols_to_entrez(gene_ranks_2week)
gene_ranks_1year <- convert_symbols_to_entrez(gene_ranks_1year)
gene_ranks_1year_adeno <- convert_symbols_to_entrez(gene_ranks_1year_adeno)

#download gmt file - I need to load this into the HPC, because it needs permissions
#to download something from the internet - 
gmt_url <- "https://www.gsea-msigdb.org/gsea/msigdb/download_file.jsp?filePath=/msigdb/release/2025.1.Mm/mh.all.v2025.1.Mm.entrez.gmt"
gmt <- "Mm.all.v2025.1.Mm.symbols.gmt"

download.file(gmt_url, gmt, mode = "wb")  # mode = "wb" for binary-safe download

gmt_file <- read.gmt(gmt)

gsea_2week <- GSEA(gene_ranks_2week, 
                   pvalueCutoff = 0.05, 
                   pAdjustMethod = "BH")


