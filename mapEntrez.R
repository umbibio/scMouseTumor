library(tidyverse)
library(openxlsx)


## Mouse-Human orthologs
## Following files were downloaded from BioMart using 
## Homologous and Features
mouse_human <- read.csv('../../data/input/mouseGenomics/mouse_human_orth.csv')
human_entrz <- read.csv('../../data/input/mouseGenomics/human_ensemble_entrez.csv')

colnames(mouse_human) <- c('mouseGeneID', 'mouseGeneName', 'humanGeneID', 'humanGeneName')
mouse_human$entrez <- human_entrz$NCBI.gene..formerly.Entrezgene..ID[match(mouse_human$humanGeneID, human_entrz$Gene.stable.ID)]

write.xlsx(mouse_human, '../../data/input/mouseGenomics/mouse_human_orth_entrez.xlsx')

three.tissue <- read_tsv('~/Desktop/three_tissue_network_subset.tsv')

FOXM1.net <- three.tissue %>% dplyr::filter(src_symbol == 'FOXM1')

FOXM1.net <- left_join(FOXM1.net , mouse_human, by = c('trg_symbol' = 'Human.gene.name'))
DEGs <- read.xlsx('../../data/input/fromJosh/upRS_LPvLP.xlsx')

FOXM1.net <- inner_join(FOXM1.net , DEGs, by = c('Gene.name' = 'GeneID'))

FOXM1.net <- FOXM1.net %>% dplyr::filter(abs(avg_log2FC) > log2(1.5))


write.xlsx(FOXM1.net, '~/Desktop/FOXM1_out.xlsx')
