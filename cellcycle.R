###script for cell cyle check#####
getwd()
#install.packages("Seurat")
library(Seurat)
#install.packages("ggplot2")
library(ggplot2)
library(sctransform)
library(tidyr)
library(dplyr)


setwd("/mnt/data/home/sarahsczelecki/single-cell/seurat/outputs")

#import R object
#load("2w1y1yadeno-final-mito45.RData")
load("2w-stroma-subset_seurat.RData")
View(subset_seurat)

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

p  <- DimPlot(subset_seurat_cc, group.by = "Phase", reduction = "umap") + 
  ggtitle("Cell Cycle Phase") +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 18),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 18))
ggsave("cellcycle-subset.png", plot = p, width = 7, height = 7, dpi = 800)

DimPlot(subset_seurat_cc, group.by = "RNA_snn_res.0.6", reduction = "umap", label = TRUE)