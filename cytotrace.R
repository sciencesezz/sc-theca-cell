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
library(CytoTRACE2)
library(pheatmap)
library(RColorBrewer)


#import R object
#load("2w3w1y-WT-merged_seurat.Robj")
load("2w-stroma-subset_seurat.RData")
View(subset_seurat)
head(subset_seurat)
tail(subset_seurat)

#############################POTENTCY#########################################
# Step 1: Prepare expression matrix for CytoTRACE
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

annotation <- data.frame(
  cell_id = colnames(expression_matrix),     # must match your expression data
  phenotype = cluster_info[ colnames(expression_matrix) ]             # or Age / cell type labels
)
head(annotation)

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

scoreplot <- FeaturePlot(subset_seurat, features = "CytoTRACE_score", 
            cols = c("lightgrey", "blue"), pt.size = 0.1) +
  labs(title = "CytoTRACE2 Score on UMAP")
print(scoreplot)
ggsave("scoreplot-by-clus-age.png", plot = scoreplot, width = 7, height = 7, dpi = 800)

potency <- DimPlot(subset_seurat, group.by = "CytoTRACE_potency", 
                   label = FALSE, repel = TRUE) +
  ggtitle("Potency") 
print(potency)
ggsave("potency.png", plot = potency, width = 7, height = 7, dpi = 800)

score <- ggplot(subset_seurat@meta.data, aes(x = RNA_snn_res.0.6, y = CytoTRACE_score, fill = RNA_snn_res.0.6)) +
  geom_violin(trim = TRUE, width = 1.8) +
  geom_boxplot(width = 0.2, outlier.size = 0.6, fill = "white") +
  facet_wrap(~ Age) +
  theme_classic() +
  labs(title = "CytoTRACE Score by Cluster and Age",
       x = "Cluster ID",
       y = "CytoTRACE Score") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        strip.text = element_text(size = 12, face = "bold"),
        legend.position = "none")
print(score)

ggsave("score-by-clus-age.png", plot = score, width = 8, height = 6, dpi = 800)
