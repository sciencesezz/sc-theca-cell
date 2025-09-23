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


#import R object - filtered, I think I can also add this to the seurat pipeline.
#load("2w3w1y-WT-merged_seurat.Robj")
load("2w-final-mito10.RData") #actually I think this is wrong I need to do it to the individual seurat objects
#then I need to merge those together

#split seurat object by sample ID, then run doublet finder on each of those datasets
#let's run 2 week as an example first
#more sensitive to heterotypic doublets - 2 diff cell types, rather than 2 cells that are the same

#need seurat object
#above

#standard workflow - basic, might not match my previous data
#merged_seurat <- NormalizeData(object = merged_seurat)
#merged_seurat <- FindVariableFeatures(object = merged_seurat)
#merged_seurat <- ScaleData(object = merged_seurat)
#merged_seurat <- RunPCA(object = merged_seurat)
#ElbowPlot(merged_seurat) #select 16 dimensions


#workflow I have used to match other scripts
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

#continue with workflow!
