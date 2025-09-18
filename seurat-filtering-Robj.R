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

setwd("/mnt/data/home/sarahsczelecki/single-cell/seurat/outputs")

##get data into R, directory should be the one that includes matrix.mtx, genes.tsv, barcodes.tsv
g1_2week <- Read10X(data.dir = 
                      "/mnt/data/home/sarahsczelecki/single-cell/2week_g1/sample_filtered_feature_bc_matrix")
g2_2week <- Read10X(data.dir = 
                      "/mnt/data/home/sarahsczelecki/single-cell/2week_g2/sample_filtered_feature_bc_matrix")
g3_2week <- Read10X(data.dir = 
                      "/mnt/data/home/sarahsczelecki/single-cell/2week_g3/sample_filtered_feature_bc_matrix")
g1_1year <- Read10X(data.dir = 
                      "/mnt/data/home/sarahsczelecki/single-cell/1year_g1/sample_filtered_feature_bc_matrix")
g2_1year <- Read10X(data.dir = 
                      "/mnt/data/home/sarahsczelecki/single-cell/1year_g2/sample_filtered_feature_bc_matrix")
g3_1year <- Read10X(data.dir = 
                      "/mnt/data/home/sarahsczelecki/single-cell/1year_g3/sample_filtered_feature_bc_matrix")
g1_3week <- Read10X(data.dir = 
                      "/mnt/data/home/sarahsczelecki/single-cell/3week_g1/sample_filtered_feature_bc_matrix")

# Turn count matrix into a Seurat object (output is a Seurat object) - 7 single cell datasets, removed RT data
seur1 <- CreateSeuratObject(counts = g1_2week,
                            project = "2week", 
                            min.features = 100, 
                            min.cells = 5)

seur2 <- CreateSeuratObject(counts = g2_2week,
                            project = "2week", 
                            min.features = 100, 
                            min.cells = 5)

seur3 <- CreateSeuratObject(counts = g3_2week,
                            project = "2week", 
                            min.features = 100, 
                            min.cells = 5)

seur4 <- CreateSeuratObject(counts = g1_3week,
                            project = "3week", 
                            min.features = 100, 
                            min.cells = 5)

seur5 <- CreateSeuratObject(counts = g1_1year,
                            project = "1year", 
                            min.features = 100, 
                            min.cells = 5)

seur6 <- CreateSeuratObject(counts = g2_1year,
                            project = "1year_adeno", 
                            min.features = 100, 
                            min.cells = 5)

seur7 <- CreateSeuratObject(counts = g3_1year,
                            project = "1year", 
                            min.features = 100, 
                            min.cells = 5)

#explore metadata
head(seur1@meta.data)
head(seur2@meta.data)
head(seur3@meta.data)
head(seur4@meta.data)
head(seur5@meta.data)
head(seur6@meta.data)
head(seur7@meta.data)

#there are three categories in the metadata orig.ident, nCount_RNA, nFeature_RNA

#merge data
merged_seurat <- merge(x = seur1, 
                       y = list(seur2, seur3), 
                       add.cell.id = c("2week_g1","2week_g2", "2week_g3"))

merged_seurat <- merge(x = seur1, 
                        y = list(seur2, seur3, seur4, seur5, seur6, seur7), 
                        add.cell.id = c("2week_g1","2week_g2", "2week_g3", "3week_g1", 
                                        "1year_g1", "1year_g2", "1year_g3"))

merged_seurat <- merge(x = seur2, 
                        y = list(seur3, seur5, seur7), 
                        add.cell.id = c("2week_g2", "2week_g3", 
                                        "1year_g1", "1year_g3"))

merged_seurat <- merge(x = seur1, 
                       y = list(seur6, seur7), 
                       add.cell.id = c("2week_g1","1year_g2", "1year_g3"))

merged_seurat <- merge(x = seur3, 
                       y = list(seur6, seur7), 
                       add.cell.id = c("2week_g3","1year_g2", "1year_g3"))

merged_seurat <- merge(x = seur2, 
                       y = list(seur4, seur5), 
                       add.cell.id = c("2week_g2","3week_g1", "1year_g1"))

merged_seurat <- merge(x = seur2, 
                       y = list(seur4, seur5, seur6), 
                       add.cell.id = c("2week_g2","3week_g1", "1year_g1", "1year_g2"))



# Concatenate the count matrices of both samples together
merged_seurat <- JoinLayers(merged_seurat)

# Check that the merged object has the appropriate sample-specific prefixes
head(merged_seurat@meta.data)
tail(merged_seurat@meta.data)
View(merged_seurat@meta.data)

#In order to create the appropriate plots for the quality control analysis, 
#we need to calculate some additional metrics. These include:
#number of genes (Features) detected per UMI: 
#this metric with give us an idea of the complexity of our dataset (more genes detected per UMI, more complex our data)
#mitochondrial ratio: 
#this metric will give us a percentage of cell reads originating from the mitochondrial genes

# Calculate mitochondrial gene percentage
merged_seurat[["percent.mt"]] <- PercentageFeatureSet(merged_seurat, pattern = "^mt-")

# Visualize QC metrics
VlnPlot(merged_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Filter out low-quality cells
merged_seurat <- subset(merged_seurat, subset = (nFeature_RNA >= 200) & 
                          (nFeature_RNA <= 6000) & 
                          (percent.mt <= 45) & 
                          (nCount_RNA >= 500) &
                          (nCount_RNA <= 50000))

#visualise QC metrics after filters
VlnPlot(merged_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Add number of genes per UMI for each cell to metadata
merged_seurat$log10GenesPerUMI <- 
  log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)

#We are a now all set with quality metrics required for assessing our data. 
#However, we would like to include some additional information 
#that would be useful to have in our metadata including: 
#cell IDs and condition information.

# Create metadata dataframe
metadata <- merged_seurat@meta.data

# Add cell IDs to metadata
metadata$cells <- rownames(metadata)
View(metadata)

# Create sample column
#install.packages("stringr")
library(stringr)
metadata$sample <- NA
metadata$sample[which(str_detect(metadata$cells, "^2week_g1_"))] <- "2week_g1"
metadata$sample[which(str_detect(metadata$cells, "^2week_g2_"))] <- "2week_g2"
metadata$sample[which(str_detect(metadata$cells, "^2week_g3_"))] <- "2week_g3"
#metadata$sample[which(str_detect(metadata$cells, "^3week_g1_"))] <- "3week_g1"
#metadata$sample[which(str_detect(metadata$cells, "^1year_g1_"))] <- "1year_g1"
#metadata$sample[which(str_detect(metadata$cells, "^1year_g2_"))] <- "1year_g2"
#metadata$sample[which(str_detect(metadata$cells, "^1year_g3_"))] <- "1year_g3"

# rename some of the existing columns in our metadata dataframe to be more intuitive
metadata <- metadata %>%
  dplyr::rename(Age = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)

# Add metadata back to Seurat object
merged_seurat@meta.data <- metadata

# Create .RData object to load at any time
save(merged_seurat, 
     file="2w-final-mito45.RData")
#########################END RDATA HERE#########################################

########################Dimensonality/Normalisation############################
levels <- c("2week_g1", "1year_g2", "1year_g3")

levels_age <- c("2week", "1year", "1year_adeno")

metadata$sample <- factor(metadata$sample, levels = levels)
metadata$Age <- factor(metadata$Age, levels = levels_age)

#check number of cells per sample
metadata %>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")

#check number of cells per age group
metadata %>% 
  ggplot(aes(x=Age, fill=Age)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")


# Visualize the number UMIs/transcripts per cell
metadata %>% 
  ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500) + 
  ggtitle("UMI counts (transcripts) per cell")

# Visualize the distribution of genes detected per cell via histogram
metadata %>% 
  ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300) + 
  ggtitle("Genes detected per cell")

# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI (novelty score)
metadata %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8) + 
  ggtitle("Genes detected per UMI (novelty score)")

# Visualize the distribution of mitochondrial gene expression detected per cell
metadata %>% 
  ggplot(aes(color=sample, x=percent.mt, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2) + 
  ggtitle("Mitocondrial gene expression detected per cell")

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
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~sample) + 
  ggtitle("Correlation between genes detected and number of UMIs per cell")

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

#plotting the resolutions - have a look at these...I chose 0.1
DimPlot(merged_seurat, group.by = "RNA_snn_res.0.2", label = TRUE)

#setting identity of clusters
Idents(merged_seurat) <- "RNA_snn_res.0.2"

#nonlinear dimenstionality reduction
merged_seurat <- RunUMAP(merged_seurat, dims = 1:16)
#merged_seurat <- RunTSNE(merged_seurat, dims = 1:20)

#export QCd file to LoupeBrowser, cannot have two files the same name, doesn't overwrite
metadata_cols <- c("RNA_snn_res.0.2", "Age")

create_loupe_from_seurat(
  merged_seurat,
  output_dir = "/mnt/data/home/sarahsczelecki/single-cell/seurat/outputs",
  output_name = "2w1y-final-mito10-2", 
  metadata_cols = c("RNA_snn_res.0.2", "Age", "sample"))
 
DimPlot(merged_seurat, reduction = "umap")
#DimPlot(filtered_seurat, reduction = "tsne")

#plot with labels - apparently you can tell if you have batch effects from this
#I don't really understand that so I'm going to have to revisit that. for now I will assume I need to 
#correct for batch effects. 
#now I understand this - if there are any cluster(s) that are distinct for the sample type that means
#that you probably have batch effects. It looks like we dont, so can continue without batch correction.
p1 <- DimPlot(merged_seurat, reduction = "umap", group.by = "RNA_snn_res.0.1", label = TRUE)
print(p1)

#Check for Batch effects graphs
levels_age <- c("2week", "1year", "1year_adeno")
merged_seurat@meta.data$Age <- factor(merged_seurat@meta.data$Age, levels = levels_age)

# Now plot
VlnPlot(merged_seurat, features = c("nUMI", "nGene", "percent.mt"), group.by = "Age")
p2 <- DimPlot(merged_seurat, reduction = "umap", group.by = "Age")
print(p2)

#correct for batch effects - will have to look up scripts if I need that...
#I think then in that case if there are batch effects, then you run the integration
#anchors part and then rerun the pipeline. Since my data doesn't have batch effects, I'm going to skip this for this
#analysis

#Add clusters - manual list generation, automated assignment
markers <- read.csv("data/home/sarahsczelecki/single-cell/seurat/test-list-broad.csv")
markers_by_cluster <- AverageExpression(merged_seurat, 
                                        features = markers$gene, 
                                        group.by = "RNA_snn_res.0.2")

head(rownames(merged_seurat))

# The result is typically a list with one element per assay
# Let's extract the RNA assay (most common)
avg_expr <- markers_by_cluster$RNA

# Convert to a more manageable format
avg_expr_df <- as.data.frame(avg_expr)

# Visualize the expression matrix as a heatmap
library(pheatmap)
pheatmap(avg_expr_df, scale = "row")

# First, make a copy of the cluster column
merged_seurat$cluster_id <- merged_seurat$RNA_snn_res.0.2

# Now create the cell_type column by direct assignment
merged_seurat$cell_type <- NA  # Initialize with NA

# Assign cell types for each cluster
merged_seurat$cell_type[merged_seurat$cluster_id == 0] <- "Fibroblastic Stroma 1"
merged_seurat$cell_type[merged_seurat$cluster_id == 1] <- "Fibroblastic Stroma 2"
merged_seurat$cell_type[merged_seurat$cluster_id == 2] <- "Endothelial 1"
merged_seurat$cell_type[merged_seurat$cluster_id == 3] <- "GCs"
merged_seurat$cell_type[merged_seurat$cluster_id == 4] <- "Fibroblastic Stroma 3"
merged_seurat$cell_type[merged_seurat$cluster_id == 5] <- "Steroidogenic Theca 1"
merged_seurat$cell_type[merged_seurat$cluster_id == 6] <- "Phagocytes 1"
merged_seurat$cell_type[merged_seurat$cluster_id == 7] <- "Epithelial 1"
merged_seurat$cell_type[merged_seurat$cluster_id == 8] <- "T Lymphocytes 1"
merged_seurat$cell_type[merged_seurat$cluster_id == 9] <- "Steroidogenic Theca 2"
merged_seurat$cell_type[merged_seurat$cluster_id == 10] <- "Epithelial 2"
merged_seurat$cell_type[merged_seurat$cluster_id == 11] <- "Pericytes"
merged_seurat$cell_type[merged_seurat$cluster_id == 12] <- "Endothelial 2"
merged_seurat$cell_type[merged_seurat$cluster_id == 13] <- "Phagocytes 2"
merged_seurat$cell_type[merged_seurat$cluster_id == 14] <- "B Lymphocytes"
merged_seurat$cell_type[merged_seurat$cluster_id == 15] <- "T Lymphocytes 2"

p <- DimPlot(merged_seurat, reduction = "umap", group.by = "cell_type", label = FALSE) + 
  theme(
    text = element_text(size = 16),           # Overall text size
    axis.text = element_text(size = 14),      # Axis tick labels
    axis.title = element_text(size = 14),     # Axis titles (umap_1, umap_2)
    legend.text = element_text(size = 14),    # Legend labels
    legend.title = element_text(size = 14)) +  # Legend title
  ggtitle("Cell Type")# Plot title
print(p)
ggsave("cell-type-umap-small.png", plot = p, width = 7, height = 5, dpi = 800)


p2 <- DimPlot(merged_seurat, reduction = "umap", group.by = "Age", label = FALSE) + 
  theme(
    text = element_text(size = 16),           # Overall text size
    axis.text = element_text(size = 14),      # Axis tick labels
    axis.title = element_text(size = 14),     # Axis titles (umap_1, umap_2)
    legend.text = element_text(size = 14),    # Legend labels
    legend.title = element_text(size = 14)) +  # Legend title
  ggtitle("Age")# Plot title
print(p2)
ggsave("age-umap-small.png", plot = p2, width = 7, height = 5, dpi = 800)


# Make a table of number of cells in each cell type by age
table(merged_seurat$cell_type, merged_seurat$Age)

# Convert to dataframe with better formatting
cell_type_age_counts <- as.data.frame(table(merged_seurat$cell_type, merged_seurat$Age))
colnames(cell_type_age_counts) <- c("CellType", "Age", "Count")

# View the result
print(cell_type_age_counts)

# First, create the count table with both cell type and age
cell_type_age_counts <- as.data.frame(table(merged_seurat$cell_type, merged_seurat$Age))
colnames(cell_type_age_counts) <- c("CellType", "Age", "Count")

# Compute percentage within each cell type (so each cell type bar sums to 100%)
library(dplyr)
cell_type_age_counts <- cell_type_age_counts %>%
  group_by(CellType) %>%
  mutate(Percentage = (Count / sum(Count)) * 100) %>%
  ungroup()

# Order cell types alphabetically
cell_type_age_counts$CellType <- factor(cell_type_age_counts$CellType, 
                                        levels = sort(unique(cell_type_age_counts$CellType)))

# Create stacked bar plot
ggplot(cell_type_age_counts, aes(x = CellType, y = Percentage, fill = Age)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = ifelse(Percentage > 5, sprintf("%.1f%%", Percentage), "")), 
            position = position_stack(vjust = 0.5), 
            size = 3, color = "white") +
  theme_minimal() +
  labs(title = "Cell Type Composition by Age", 
       x = "Cell Type", y = "Percentage (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("2week" = "#F8766D", 
                               "1year" = "#00BA38", 
                               "1year_adeno" = "#619CFF")) +
  ylim(0, 100)

#############separate fibroblast cell types
# First, create the count table with both cell type and age
cell_type_age_counts <- as.data.frame(table(merged_seurat$cell_type, merged_seurat$Age))
colnames(cell_type_age_counts) <- c("CellType", "Age", "Count")

# Filter to only include Fibroblastic Stroma 1, 2, and 3
fibroblast_types <- c("Fibroblastic Stroma 1", "Fibroblastic Stroma 2", "Fibroblastic Stroma 3")
cell_type_age_counts <- cell_type_age_counts[cell_type_age_counts$CellType %in% fibroblast_types, ]

# Compute percentage within each cell type (so each cell type bar sums to 100%)
library(dplyr)
cell_type_age_counts <- cell_type_age_counts %>%
  group_by(CellType) %>%
  mutate(Percentage = (Count / sum(Count)) * 100) %>%
  ungroup()

# Order cell types alphabetically (or numerically: 1, 2, 3)
cell_type_age_counts$CellType <- factor(cell_type_age_counts$CellType, 
                                        levels = sort(fibroblast_types))

# Create stacked bar plot
ggplot(cell_type_age_counts, aes(x = CellType, y = Percentage, fill = Age)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = ifelse(Percentage > 5, sprintf("%.1f%%", Percentage), "")), 
            position = position_stack(vjust = 0.5), 
            size = 5, color = "black") +
  theme_minimal() +
  labs(title = "Fibroblastic Stroma Composition by Age", 
       x = "Cell Type", y = "Percentage (%)") +
  theme(
    text = element_text(size = 16),              # Overall base font size
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),  # X-axis labels
    axis.text.y = element_text(size = 14),       # Y-axis labels
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 16),# Axis titles
    plot.title = element_text(size = 18),        # Plot title
    legend.text = element_text(size = 14),       # Legend labels
    legend.title = element_text(size = 16)       # Legend title
  ) +
  scale_fill_manual(values = c("2week" = "#F8766D", 
                               "1year" = "#00BA38", 
                               "1year_adeno" = "#619CFF")) +
  ylim(0, 100)


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

cluster_markers_80high <- FindAllMarkers(merged_seurat, test.use = "MAST", 
                                         only.pos = TRUE, min.pct = 0.75, logfc.threshold = 2.0)

View(cluster_markers_80high)
write.csv(cluster_markers_80high, "2w1y_cluster_markers_75high.csv", row.names = FALSE)

library(patchwork)

# Create individual plots
enpep_plot <- FeaturePlot(merged_seurat, features = "Enpep")
hhip_plot <- FeaturePlot(merged_seurat, features = "Hhip")
col1a1_plot <- FeaturePlot(merged_seurat, features = "Col1a1")
age_plot <- DimPlot(merged_seurat, group.by = "Age")
cluster_plot <- DimPlot(merged_seurat, group.by = "cell_type") + 
  theme(legend.text = element_text(size = 7),
        legend.key.size = unit(0.3, "cm"))

plot1 <- (enpep_plot | plot_spacer()) / (hhip_plot | cluster_plot) / (col1a1_plot | age_plot) 
print(plot1)
ggsave("plot-enpep-age.png", plot1, width = 9, height = 12, dpi = 800)


DimPlot(merged_seurat, reduction = "umap", group.by = "cluster_id", label = TRUE) +
  ggtitle("By Cluster ID")

#plots for CAF markers
vim_plot <- FeaturePlot(merged_seurat, features = "Vim")
ccn1_plot <- FeaturePlot(merged_seurat, features = "Ccn1")
col1a1_plot <- FeaturePlot(merged_seurat, features = "Col1a1")
col1a2_plot <- FeaturePlot(merged_seurat, features = "Col1a2")
ccl19_plot <- FeaturePlot(merged_seurat, features = "Ccl19")
s100_plot <- FeaturePlot(merged_seurat, feature = "S100a4")
acta2_plot <- FeaturePlot(merged_seurat, feature = "Acta2")
fap_plot <- FeaturePlot(merged_seurat, feature = "Fap")
pdgfra_plot <- FeaturePlot(merged_seurat, feature = "Pdgfra")
pdgfrb_plot <- FeaturePlot(merged_seurat, feature = "Pdgfrb")

age_plot <- DimPlot(merged_seurat, group.by = "Age")
cluster_plot <- DimPlot(merged_seurat, group.by = "cell_type") + 
  theme(legend.text = element_text(size = 7),
        legend.key.size = unit(0.3, "cm"))

plot1 <- ((vim_plot | ccn1_plot) / (col1a1_plot | col1a2_plot) / (ccl19_plot | s100_plot) / (age_plot | cluster_plot))
print(plot1)
ggsave("caf-stroma.png", plot1, width = 12, height = 12, dpi = 800)

plot2 <- ((acta2_plot | fap_plot) / (pdgfra_plot | pdgfrb_plot) / (age_plot | cluster_plot))
print(plot2)
ggsave("caf-stroma2.png", plot2, width = 12, height = 12, dpi = 800)

#I wonder if you can annotate this heatmap, so that you only include your markers of interest.

cluster_markers_80high %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 2) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5
DoHeatmap(merged_seurat, features = top5$gene) + NoLegend()

#define genes by functional categories to annotate with your own selected markers
#define genes by functional categories to annotate with your own selected markers
#define genes by functional categories to annotate with your own selected markers
bcells <- "Cd79a"
tcells <- "Cd3g"
endothelial <- c("Cd93", "Cdh5", "Pecam1")
epithelial <- c("Krt7", "Krt18", "Lgals7", "Upk3b", "Muc16")
pericytes <- c("Rgs5", "Notch3", "Rgs4")
fibroblast <- c("Col1a1", "Dcn", "Mfap4")
granulosa <- c("Amh", "Cyp19a1", "Inha", "Hsd17b1")
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
               group.bar.height = 0.03, 
               label = FALSE  # Changed from FALSE to TRUE
) +
  theme(axis.text.y = element_text(size = 16),
        plot.margin = margin(1, 1, 1, 1),
        legend.text = element_text(size = 14), 
        legend.title = element_text(size = 16)) +
  guides(fill = guide_colorbar())
print(p)

# Save the plot
ggsave("heatmap.png", p, width = 12, height = 8, dpi = 1200)

#---------------------------------subset for stromal cells of interest------------
#subset the seurat object
Idents(merged_seurat) <- "cluster_id"
subset_seurat <- subset(merged_seurat, idents = c("0","1","4"))


my_colors <- c("Fibroblastic Stroma 1" = "#0BB702", 
               "Fibroblastic Stroma 2" = "#00BE67", 
               "Fibroblastic Stroma 3" = "#00C19A")
DimPlot(subset_seurat, group.by = "cell_type", cols = my_colors)

DimPlot(subset_seurat, group.by = "Age")

# Re-run the analysis pipeline on the subset
subset_seurat <- NormalizeData(subset_seurat, 
                               normalization.method = "LogNormalize", 
                               scale.factor = 10000)
subset_seurat <- FindVariableFeatures(subset_seurat, selection.method = "vst", nfeatures = 2000)
subset_seurat <- ScaleData(subset_seurat, features = rownames(merged_seurat))
subset_seurat <- RunPCA(subset_seurat, features = VariableFeatures(object = merged_seurat))
ElbowPlot(subset_seurat)

subset_seurat <- FindNeighbors(subset_seurat, dims = 1:17) #17 based in elbow plot
subset_seurat <- FindClusters(subset_seurat, resolution = c(0.1, 0.3, 0.5, 0.6, 0.8, 0.7, 1))
DimPlot(subset_seurat, group.by = "RNA_snn_res.0.6", label = TRUE)
DimPlot(subset_seurat, group.by = "Age", label = TRUE)
#setting identity of clusters
Idents(subset_seurat) <- "RNA_snn_res.0.6"

#nonlinear dimensionality reduction
subset_seurat <- RunUMAP(subset_seurat, dims = 1:15) #not sure if I need to also set this to 16? check
#merged_seurat <- RunTSNE(merged_seurat, dims = 1:20)

p <- DimPlot(subset_seurat, group.by = "RNA_snn_res.0.6", label = TRUE) + 
  theme(
    text = element_text(size = 18),           # Overall text size
    axis.text = element_text(size = 16),      # Axis tick labels
    axis.title = element_text(size = 16),     # Axis titles (umap_1, umap_2)
    legend.text = element_text(size = 16),    # Legend labels
    legend.title = element_text(size = 16)) +  # Legend title
  ggtitle("Fibroblastic Stroma Subset Clusters")# Plot title
print(p)
ggsave("subclusters-umap-small.png", plot = p, width = 7, height = 5, dpi = 800)


p2 <- DimPlot(subset_seurat, reduction = "umap", group.by = "Age", label = FALSE) + 
  theme(
    text = element_text(size = 18),           # Overall text size
    axis.text = element_text(size = 16),      # Axis tick labels
    axis.title = element_text(size = 16),     # Axis titles (umap_1, umap_2)
    legend.text = element_text(size = 16),    # Legend labels
    legend.title = element_text(size = 16)) +  # Legend title
  ggtitle("Age")# Plot title
print(p2)
ggsave("subsetage-umap-small.png", plot = p2, width = 7, height = 5, dpi = 800)


library(patchwork)

# Create individual plots
enpep_plot2 <- FeaturePlot(subset_seurat, features = "Enpep")
hhip_plot2 <- FeaturePlot(subset_seurat, features = "Hhip")
col1a1_plot2 <- FeaturePlot(subset_seurat, features = "Col1a1")
age_plot2 <- DimPlot(subset_seurat, group.by = "Age")
cluster_plot2 <- DimPlot(subset_seurat, group.by = "RNA_snn_res.0.6") + 
  theme(legend.text = element_text(size = 7),
        legend.key.size = unit(0.3, "cm"))

plot1 <- (enpep_plot2 | plot_spacer()) / (hhip_plot2 | cluster_plot2) / (col1a1_plot2 | age_plot2) 
print(plot1)
ggsave("plot-enpep-age-subset.png", plot1, width = 9, height = 12, dpi = 800)


#export QCd file to LoupeBrowser
create_loupe_from_seurat(
  subset_seurat,
  output_dir = "/mnt/RosalindData/sarah/single-cell/seurat/outputs",
  output_name = "2week-1year-stroma-seurat-subset")

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

#############POTENTCY OF STROMA#########################################
# Load required libraries
library(Seurat)
#devtools::install_github("digitalcytometry/cytotrace2", subdir = "cytotrace2_r")
library(CytoTRACE2)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(RColorBrewer)

# Step 2: Prepare expression matrix for CytoTRACE
# CytoTRACE works with raw counts (not normalized data)
expression_matrix <- GetAssayData(subset_seurat, assay = "RNA", slot = "counts")

# Convert sparse matrix to regular matrix (required by CytoTRACE)
expression_matrix <- as.matrix(expression_matrix)

# Check matrix dimensions
print(paste("Expression matrix dimensions:", 
            nrow(expression_matrix), "genes x", 
            ncol(expression_matrix), "cells"))


# Step 3: Extract metadata
cell_metadata <- subset_seurat@meta.data
cluster_info <- cell_metadata$RNA_snn_res.0.6 #need to change this to the correct clusters
names(cluster_info) <- rownames(cell_metadata)

# Step 4: Run CytoTRACE analysis
print("Running CytoTRACE analysis...")
print("This may take several minutes depending on dataset size...")

# Run CytoTRACE with default parameters
cytotrace_results <- cytotrace2(expression_matrix)

# generate prediction and phenotype association plots with plotData function
plots <- plotData(cytotrace2_result = cytotrace_results, 
                  annotation = annotation,
                  expression_data = expression_matrix
)

plots$CytoTRACE2_Potency_UMAP

# Step 5: Extract CytoTRACE scores
cytotrace_scores <- cytotrace_results$CytoTRACE2_Score
names(cytotrace_scores) <- rownames(cytotrace_results)  # These must match cell names in Seurat
subset_seurat$CytoTRACE_score <- cytotrace_scores[Cells(subset_seurat)]

subset_seurat$CytoTRACE_potency <- cytotrace_results$CytoTRACE2_Potency

FeaturePlot(subset_seurat, features = "CytoTRACE_score", 
            cols = c("lightgrey", "blue"), pt.size = 0.5) +
  labs(title = "CytoTRACE2 Score on UMAP")

age_potency <- DimPlot(subset_seurat, group.by = "Age", 
        label = FALSE, repel = TRUE) +
  ggtitle("Age") 
print(age_potency)

potency <- DimPlot(subset_seurat, group.by = "CytoTRACE_potency", 
               label = FALSE, repel = TRUE) +
  ggtitle("Potency") 
print(potency)

potency_combo <- potency / age_potency
print(potency_combo)
ggsave("potency.png", plot = potency, width = 7, height = 5, dpi = 800)



###########################slingshot - single trajectory#########################
library(Seurat)
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
                 start.clus = '0')  # Replace '1' with your starting cluster

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
library(Seurat)
library(slingshot)
library(SingleCellExperiment)
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
                 start.clus = '0',
                 stretch = 0.8,
                 smoother = 'smooth.spline',
                 shrink = 1.0)

pt <- slingPseudotime(sds)
n_lineages <- ncol(pt)

# Calculate layout dimensions
n_cols <- ceiling(sqrt(n_lineages + 1))
n_rows <- ceiling((n_lineages + 1) / n_cols)

# Setup the plotting device with your adjusted margins
par(mfrow = c(n_rows, n_cols), 
    mar = c(3, 2, 1, 4),    # You can adjust these values as needed
    oma = c(2, 2, 2, 0))  

colors_pt <- colorRampPalette(c("grey", "lightblue", "navyblue"))(100)

# Plot each lineage
for(i in 1:n_lineages) {
  current_pt <- pt[,i]
  pt_range <- range(current_pt, na.rm = TRUE)
  
  plot(reducedDims(sce)$UMAP, 
       col = colors_pt[cut(current_pt, breaks = 100)],
       pch = 16, 
       cex = 0.5,
       xlab = "UMAP 1",     # Changed to explicit UMAP labels
       ylab = "UMAP 2",
       main = paste("Lineage", i),
       xaxt = "s",          # Ensure x-axis is drawn
       yaxt = "s")          # Ensure y-axis is drawn
  
  # Add curves
  lines(SlingshotDataSet(sds), 
        lwd = 2,
        col = 'black',
        type = 'curves')
  
  lines(SlingshotDataSet(sds), 
        lwd = 3,
        col = 'red',
        type = 'curves',
        linInd = i)
  
  # Add pseudotime legend
  par(xpd = TRUE)
  legend_image <- as.raster(matrix(rev(colors_pt), ncol = 1))
  rasterImage(legend_image, 
              xleft = par("usr")[2] + 0.1 * diff(par("usr")[1:2]), 
              ybottom = par("usr")[3], 
              xright = par("usr")[2] + 0.2 * diff(par("usr")[1:2]), 
              ytop = par("usr")[4])
  
  # Add numerical labels
  text(x = par("usr")[2] + 0.25 * diff(par("usr")[1:2]),
       y = par("usr")[4],
       labels = round(pt_range[2], 1),
       cex = 0.7)
  text(x = par("usr")[2] + 0.25 * diff(par("usr")[1:2]),
       y = par("usr")[3],
       labels = round(pt_range[1], 1),
       cex = 0.7)
  text(x = par("usr")[2] + 0.25 * diff(par("usr")[1:2]),
       y = mean(par("usr")[3:4]),
       labels = "Pseudotime",
       srt = 90,
       cex = 0.8)
  par(xpd = FALSE)
}

# Summary plot
if (n_rows * n_cols > n_lineages) {
  plot(reducedDims(sce)$UMAP, 
       col = rainbow(length(unique(clusters)))[as.numeric(colData(sce)$GMM)],
       pch = 16, 
       cex = 0.5,
       xlab = "UMAP 1",     # Changed to explicit UMAP labels
       ylab = "UMAP 2",
       main = "All Trajectories",
       xaxt = "s",          # Ensure x-axis is drawn
       yaxt = "s")          # Ensure y-axis is drawn
  
  lines(SlingshotDataSet(sds), 
        lwd = 2,
        col = 'black',
        type = 'curves')
  
  # Add cluster legend
  par(xpd = TRUE)
  legend(x = par("usr")[2] + 0.1 * diff(par("usr")[1:2]),
         y = par("usr")[4],
         legend = levels(clusters), 
         col = rainbow(length(unique(clusters))),
         pch = 16,
         cex = 0.7,
         title = "Clusters")
  par(xpd = FALSE)
}

# Reset par settings
par(mfrow = c(1,1), mar = c(5,4,4,2) + 0.1, oma = c(0,0,0,0))

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
root_cells <- colnames(subset_cds)[clusters(subset_cds) == 0]

subset_cds <- order_cells(subset_cds, reduction_method = "UMAP", 
                          root_cells = root_cells)

plot_cells(subset_cds, 
           color_cells_by = "pseudotime", 
           label_groups_by_cluster = FALSE,
           show_trajectory_graph = TRUE,
           trajectory_graph_color = "grey80", 
           trajectory_graph_segment_size = 1.5,
           label_branch_points = FALSE, 
           label_roots = FALSE, 
           label_leaves = FALSE, 
           cell_size = 1.0, 
           graph_label_size = 4)
