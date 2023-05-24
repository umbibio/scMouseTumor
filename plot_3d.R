library(plotly)
library(Seurat)
library(tidyverse)
library(openxlsx)


source('./util_funcs.R')

getPcaUmapMetaData <- function(S.O){
  pc <- S.O@reductions$pca@cell.embeddings
  pc <- data.frame(pc) %>% dplyr::mutate(Sample = rownames(pc)) %>% 
    transmute(Sample = Sample, PC_1 = PC_1, PC_2 = PC_2, PC_3 = PC_3)
  umap <- S.O@reductions$umap@cell.embeddings
  umap <- data.frame(umap) %>% dplyr::mutate(Sample = rownames(umap)) %>% 
    transmute(Sample = Sample, UMAP_1 = UMAP_1, UMAP_2 = UMAP_2, UMAP_3 = UMAP_3)
  
  # meta.data <- data.frame(Sample = rownames(S.O@meta.data), customclassif = S.O@meta.data$customclassif, 
  #                        clusters = S.O@meta.data$seurat_clusters)
  meta.data <- data.frame(Sample = rownames(S.O@meta.data), celltype = S.O@meta.data$cell.types, 
                          clusters = S.O@meta.data$seurat_clusters)
  
  meta.data <- left_join(meta.data,
                         pc, by = 'Sample')
  meta.data <- left_join(meta.data, umap, by = 'Sample')
  return(meta.data)  
}


S.O <- readRDS('../../data/input/fromJosh/s2_fvb50subset_new.rds')
S.O@meta.data$cell.types <- S.O@active.ident

pca.tg <- getPcaUmapMetaData(S.O)




fig <- plot_ly(pca.tg, x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, color = ~celltype)
               #colors = c('#f8766d', '#a2a400', '#00bf7d', '#04b0f6', '#e76bf3')) 
fig <- fig %>% add_markers(size = 2)
fig <- fig %>% layout(scene = list(xaxis = list(title = 'UMAP_1'),
                                   yaxis = list(title = 'UMAP_2'),
                                   zaxis = list(title = 'UMAP_3')))

fig

Idents(S.O) <- 'celltype'
FeaturePlot(S.O, features = c("Pdgfra", 'Pdgfrb'), reduction = 'umap', label = T)


PCAs <- pca.tg %>% select(contains('PC'))
PCAs$type <- 'cells'
PCAs$color <- pca.tg$phase

traces <- sds.data %>% dplyr::select(contains('PC'))
traces$type <- 'curve'
traces$color <- 'trace'

df <- bind_rows(PCAs, traces)
df$color <- factor(df$color, levels = c('G1.a', 'G1.b', 'S', 'M', 'C', 'trace'))
df$type <- factor(df$type, c('cells', 'curve'))

fig <- plot_ly(df, x = ~PC_1, y = ~PC_2, z = ~PC_3, color = ~color, 
               colors = c('#f8766d', '#a2a400', '#00bf7d', '#04b0f6', '#e76bf3', 'black')) 
fig <- fig %>% add_markers(size = c(4,4)[df$type])
fig <- fig %>% layout(scene = list(xaxis = list(title = 'PC_1'),
                                   yaxis = list(title = 'PC_2'),
                                   zaxis = list(title = 'PC_3')))

fig

plot(sds.data$PC_1[sds.data$cell.ord], sds.data$PC_2[sds.data$cell.ord], type = 'l')


#### isoMap
library(RDRToolbox)
iso <- Isomap(as.matrix(S.O@assays$RNA@data), dims=3, k=10, plotResiduals=TRUE)
plot(iso$dim2)

