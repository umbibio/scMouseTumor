library(RColorBrewer)
library(Seurat)
library(tidyverse)
library(openxlsx)



#library(sctransform)

source('./util_funcs.R')

num.cores <- detectCores(all.tests = FALSE, logical = TRUE)

S.O <- readRDS('../../data/input/fromJosh/s2_fvb50subset_new.rds')
S.O@meta.data$cell.types <- S.O@active.ident

## Important contrasts to make:
## 1. RS-LP vs LP
## 2. Tumor vs RS_LP 
## 3. RS_LP to B
## 4. Tumor vs LP

# ## Differential gene expression
Idents(S.O) <- 'cell.types'
RS_LP.vs.LP <- FindMarkers(object = S.O, ident.1 = 'RS-LP', ident.2 = 'LP')

RS_LP.vs.LP$GeneID <- rownames(RS_LP.vs.LP)
RS_LP.vs.LP <- RS_LP.vs.LP %>% dplyr::filter(abs(avg_log2FC) > log2(1.5) & p_val_adj < 0.05)
RS_LP.vs.LP <- RS_LP.vs.LP %>% arrange(desc(abs(avg_log2FC)))
FeaturePlot(object = S.O, 
            features = c('Krt14', 'Trf'), 
            cols = c("grey", "blue"), reduction = "umap")


saveRDS(RS_LP.vs.LP, '../../data/input/rds/RS_LP.vs.LP_DEGs.rds')



Tumor.vs.RS_LP <- FindMarkers(object = S.O, ident.1 = 'tumor', ident.2 = 'RS-LP')

Tumor.vs.RS_LP$GeneID <- rownames(Tumor.vs.RS_LP)
Tumor.vs.RS_LP <- Tumor.vs.RS_LP %>% dplyr::filter(abs(avg_log2FC) > log2(1.5) & p_val_adj < 0.05)
Tumor.vs.RS_LP <- Tumor.vs.RS_LP %>% arrange(desc(abs(avg_log2FC)))
FeaturePlot(object = S.O, 
            features = c('Krt6a', 'Ptn'), 
            cols = c("grey", "blue"), reduction = "umap")


saveRDS(Tumor.vs.RS_LP, '../../data/input/rds/Tumor.vs.RS_LP_DEGs.rds')


Tumor.vs.LP <- FindMarkers(object = S.O, ident.1 = 'tumor', ident.2 = 'LP')

Tumor.vs.LP$GeneID <- rownames(Tumor.vs.LP)
Tumor.vs.LP <- Tumor.vs.LP %>% dplyr::filter(abs(avg_log2FC) > log2(1.5) & p_val_adj < 0.05)
Tumor.vs.LP <- Tumor.vs.LP %>% arrange(desc(abs(avg_log2FC)))
FeaturePlot(object = S.O, 
            features = c('Krt6a', 'Ptn'), 
            cols = c("grey", "blue"), reduction = "umap")


saveRDS(Tumor.vs.LP, '../../data/input/rds/Tumor.vs.LP_DEGs.rds')


RS_LP.vs.B <- FindMarkers(object = S.O, ident.1 = 'RS-LP', ident.2 = 'B')

RS_LP.vs.B$GeneID <- rownames(RS_LP.vs.B)
RS_LP.vs.B <- RS_LP.vs.B %>% dplyr::filter(abs(avg_log2FC) > log2(1.5) & p_val_adj < 0.05)
RS_LP.vs.B <- RS_LP.vs.B %>% arrange(desc(abs(avg_log2FC)))
FeaturePlot(object = S.O, 
            features = c('Tagln', '2810417H13Rik'), 
            cols = c("grey", "blue"), reduction = "umap")


saveRDS(RS_LP.vs.B, '../../data/input/rds/RS_LP.vs.B_DEGs.rds')
