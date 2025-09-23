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
library(patchwork)
library(pheatmap)
library(RColorBrewer)
#remotes::install_github('chris-mcginnis-ucsf/DoubletFinder', force = TRUE)
library(DoubletFinder)


#import R object
#load("2w3w1y-WT-merged_seurat.Robj")
load("2w-final-mito10.RData") #actually I think this is wrong I need to do it to the individual seurat objects
#then I need to merge those together

#split seurat object by sample ID, then run doublet finder on each of those datasets
#let's run 2 week as an example first
#more sensitive to heterotypic doublets - 2 diff cell types, rather than 2 cells that are the same

#need seurat object
#above

#standard workflow
#merged_seurat <- NormalizeData(object = merged_seurat)
#merged_seurat <- FindVariableFeatures(object = merged_seurat)
#merged_seurat <- ScaleData(object = merged_seurat)
#merged_seurat <- RunPCA(object = merged_seurat)
#ElbowPlot(merged_seurat) #select 16 dimensions


#workflow I have used
##Step 3 - normalisation
merged_seurat <- NormalizeData(merged_seurat, 
                               normalization.method = "LogNormalize", 
                               scale.factor = 10000)

## Step 4 - identify highly variable features
#only want to select a few features that exhibit high cell to cell variation
merged_seurat <- FindVariableFeatures(merged_seurat, selection.method = "vst", nfeatures = 2000)

#Step 5 - now we scale the data.
merged_seurat <- ScaleData(merged_seurat, features = rownames(merged_seurat))

#Step 6 - Perform linear dimensionality reduction
merged_seurat <- RunPCA(merged_seurat, features = VariableFeatures(object = merged_seurat))
ElbowPlot(merged_seurat) #select 16 dimensions

#Step 7 - clustering, select dimensions based on elbowplot above
merged_seurat <- FindNeighbors(merged_seurat, dims = 1:16)

#understanding resolution - higher the number the more clusters, putting multiple checks multiple
#resolutions at once to find the best
merged_seurat <- FindClusters(merged_seurat, resolution = c(0.1, 0.2, 0.3, 0.5, 0.6, 0.7))
#plotting the resolutions - have a look at these...I chose 0.2
DimPlot(merged_seurat, group.by = "RNA_snn_res.0.2", label = TRUE)
#setting identity of clusters
Idents(merged_seurat) <- "RNA_snn_res.0.2"

#nonlinear dimenstionality reduction
merged_seurat <- RunUMAP(merged_seurat, dims = 1:16)
DimPlot(merged_seurat, group.by = "RNA_snn_res.0.2", label = TRUE)

#pK value optimisation (no ground truth) - from github page
sweep_res  <- paramSweep(merged_seurat, PCs = 1:16, sct = FALSE)
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

#pN





