library(Seurat)
library(HGNChelper)
library(openxlsx)
library(ggplot2)
library(ggfortify)
library(jcolors)
require(gridExtra)
library(grid)
library(matrixStats)
library(tidyverse)
library(RColorBrewer)
library(parallel)
library(openxlsx)
library(plotly)


source('./util_funcs.R')

## Prepare Marker list
markers <- read_tsv('../../data/input/fromJosh/MammaryMarkers.tsv', col_names = F)
header.ind <- grep('>', markers$X1)
markers.list <- list()
for(i in 1:(length(header.ind) - 1)){
  markers.list <- c(markers.list, markers$X1[(header.ind[i]+1):(header.ind[i+1]-1)])
}
markers.list <- c(markers.list, markers$X1[length(header.ind) + 1])
marker.groups <- gsub('>', '', markers$X1[header.ind])
names(markers.list) <- marker.groups
markers.list <- lapply(markers.list, function(ll){
  tmp <- strsplit(ll, split = ' ')[[1]]
  return(tmp[tmp !=''])
})

saveRDS(markers.list, '../../data/input/rds/markers.list')
markers.list <- readRDS('../../data/input/rds/markers.list')

num.cores <- detectCores(all.tests = FALSE, logical = TRUE)


T0.data <- Read10X(data.dir = "../../data/input/scRNA-seq/processedCounts/T0/filtered_feature_bc_matrix/")
T1.data <- Read10X(data.dir = "../../data/input/scRNA-seq/processedCounts/T1/filtered_feature_bc_matrix/")
T2.data <- Read10X(data.dir = "../../data/input/scRNA-seq/processedCounts/T2/filtered_feature_bc_matrix/")
T3.data <- Read10X(data.dir = "../../data/input/scRNA-seq/processedCounts/T3/filtered_feature_bc_matrix/")
T4.data <- Read10X(data.dir = "../../data/input/scRNA-seq/processedCounts/T4/filtered_feature_bc_matrix/")


S.O.T0 <- CreateSeuratObject(counts = T0.data,  min.cells = 3, min.features = 200)
S.O.T0[["percent.mt"]] <- PercentageFeatureSet(S.O.T0, pattern = "^mt-")
plot1 <- FeatureScatter(S.O.T0, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(S.O.T0, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1
plot2
S.O.T0 <- subset(S.O.T0, subset = nFeature_RNA > 150 & nFeature_RNA < 6100 & percent.mt < 18)
S.O.T0 <- prep_S.O(S.O.T0)
S.O.T0 <- addCellId(S.O.T0, markers.list)
S.O.T0@meta.data$data.id <- 'T0'
DimPlot(S.O.T0, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'customclassif', dims = c(1,3))        

S.O.T1 <- CreateSeuratObject(counts = T1.data,  min.cells = 3, min.features = 200)
S.O.T1[["percent.mt"]] <- PercentageFeatureSet(S.O.T1, pattern = "^mt-")
plot1 <- FeatureScatter(S.O.T1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(S.O.T1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1
plot2
S.O.T1 <- subset(S.O.T1, subset = nFeature_RNA > 150 & nFeature_RNA < 6100 & percent.mt < 15)
S.O.T1 <- prep_S.O(S.O.T1)
S.O.T1 <- addCellId(S.O.T1, markers.list)
S.O.T1@meta.data$data.id <- 'T1'
DimPlot(S.O.T1, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'customclassif')        

S.O.T2 <- CreateSeuratObject(counts = T2.data,  min.cells = 3, min.features = 200)
S.O.T2[["percent.mt"]] <- PercentageFeatureSet(S.O.T2, pattern = "^mt-")
plot1 <- FeatureScatter(S.O.T2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(S.O.T2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#plot1
#plot2
S.O.T2 <- subset(S.O.T2, subset = nFeature_RNA > 150 & nFeature_RNA < 6100 & percent.mt < 40)
S.O.T2 <- prep_S.O(S.O.T2)
S.O.T2 <- addCellId(S.O.T2, markers.list)
S.O.T2@meta.data$data.id <- 'T2'
DimPlot(S.O.T2, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'customclassif', dims = c(1,3))        

S.O.T3 <- CreateSeuratObject(counts = T3.data,  min.cells = 3, min.features = 200)
S.O.T3[["percent.mt"]] <- PercentageFeatureSet(S.O.T3, pattern = "^mt-")
plot1 <- FeatureScatter(S.O.T3, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(S.O.T3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1
plot2
S.O.T3 <- subset(S.O.T3, subset = nFeature_RNA > 150 & nFeature_RNA < 10000 & percent.mt < 40)
S.O.T3 <- prep_S.O(S.O.T3)
S.O.T3 <- addCellId(S.O.T3, markers.list)
S.O.T3@meta.data$data.id <- 'T3'
DimPlot(S.O.T3, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'customclassif')        


S.O.T4 <- CreateSeuratObject(counts = T4.data,  min.cells = 3, min.features = 200)
S.O.T4[["percent.mt"]] <- PercentageFeatureSet(S.O.T4, pattern = "^mt-")
plot1 <- FeatureScatter(S.O.T4, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(S.O.T4, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1
plot2
S.O.T4 <- subset(S.O.T4, subset = nFeature_RNA > 150 & nFeature_RNA < 8000 & percent.mt < 30)
S.O.T4 <- prep_S.O(S.O.T4)
S.O.T4 <- addCellId(S.O.T4, markers.list)
S.O.T4@meta.data$data.id <- 'T4'
DimPlot(S.O.T4, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'customclassif')        


S.O.list <- list(T0 = S.O.T0, T1 = S.O.T1, T2 = S.O.T2, T3 = S.O.T3, T4 = S.O.T4)

# ## Downsample to 6000 cells
# set.seed(100)
# S.O.list <- mclapply(S.O.list, function(S.O){
#   S.O <- subset(x = S.O, downsample = 8000)
#   
# }, mc.cores = num.cores)



### Merging data
alldata.lab <- merge(S.O.list[[1]], S.O.list[2:5], add.cell.ids=c("T0","T1","T2","T3", "T4"))


### Anchoring data
## shared PCA/UMAP
S.Os <- SplitObject(alldata.lab, split.by = 'data.id')
S.O.list <- lapply(X = S.Os, FUN = function(x) {
  ## Extract the count data
  
  ## extract the count data from each as.matrix(S.O.list[[1]][["RNA"]]@data)
  ## Replace genes with Bdiv orthologous when needed
  ## recreate the new Seurat object.
  DefaultAssay(x) <- 'RNA'
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 6000)
})

features <- SelectIntegrationFeatures(object.list = S.O.list, nfeatures = 6000)
anchors <- FindIntegrationAnchors(object.list = S.O.list, 
                                  anchor.features = features)
S.O.integrated <- IntegrateData(anchorset = anchors)
# switch to integrated assay. Make sure to set to RNA for Differential Expression
DefaultAssay(S.O.integrated) <- "integrated"
S.O.integrated <- ScaleData(object = S.O.integrated, verbose = FALSE)
S.O.integrated <- RunPCA(S.O.integrated, features = VariableFeatures(object = S.O.integrated))
S.O.integrated <- FindNeighbors(S.O.integrated, dims = 1:10, reduction = 'pca')
S.O.integrated <- FindClusters(S.O.integrated, resolution = 0.4)
S.O.integrated <- RunUMAP(S.O.integrated, dims = 1:13)

## Test plots
Idents(S.O.integrated) <- 'phase'
p <- DimPlot(S.O.integrated, reduction = "umap", 
             #group.by = "cell", 
             #split.by = 'spp',
             pt.size = 1,
             #shape.by='spp',
             label = F, label.size = 4) + #NoLegend() + 
  theme(panel.spacing = unit(0.5, "lines")) + 
  theme(axis.text.x = element_text(face="bold", size=12, angle=0)) +
  theme(axis.text.y = element_text(face="bold", size=12, angle=0)) +
  theme(
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  )


plot(p)

saveRDS(S.O.integrated, '../Input/toxo_cdc/rds/S.O.intra_extra_crk2_ark3_lables_anchored_integrated.rds')
