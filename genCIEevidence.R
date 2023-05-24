library(tidyverse)
library(openxlsx)


## Entrez IDS
mouse_human <- read.xlsx('../../data/input/mouseGenomics/mouse_human_orth_entrez.xlsx')

## DEGs
RS_LP.vs.LP <- readRDS('../../data/input/rds/RS_LP.vs.LP_DEGs.rds')
RS_LP.vs.LP <- inner_join(RS_LP.vs.LP, mouse_human, by = c('GeneID' = 'mouseGeneName'))
RS_LP.vs.LP.evidence <- RS_LP.vs.LP %>% transmute(entrez, FC = avg_log2FC, pvalue = p_val)
write.csv(RS_LP.vs.LP.evidence, '../../data/input/CIE/RS_LP.vs.LP.evidence.csv')


Tumor.vs.RS_LP <- readRDS('../../data/input/rds/Tumor.vs.RS_LP_DEGs.rds')
Tumor.vs.RS_LP <- inner_join(Tumor.vs.RS_LP, mouse_human, by = c('GeneID' = 'mouseGeneName'))
Tumor.vs.RS_LP.evidence <- Tumor.vs.RS_LP %>% transmute(entrez, FC = avg_log2FC, pvalue = p_val)
write.csv(Tumor.vs.RS_LP.evidence, '../../data/input/CIE/Tumor.vs.RS_LP.evidence.csv')


Tumor.vs.LP <- readRDS('../../data/input/rds/Tumor.vs.LP_DEGs.rds')
Tumor.vs.LP <- inner_join(Tumor.vs.LP, mouse_human, by = c('GeneID' = 'mouseGeneName'))
Tumor.vs.LP.evidence <- Tumor.vs.LP %>% transmute(entrez, FC = avg_log2FC, pvalue = p_val)
write.csv(Tumor.vs.LP.evidence, '../../data/input/CIE/Tumor.vs.LP.evidence.csv')


RS_LP.vs.B <- readRDS('../../data/input/rds/RS_LP.vs.B_DEGs.rds')
RS_LP.vs.B <- inner_join(RS_LP.vs.B, mouse_human, by = c('GeneID' = 'mouseGeneName'))
RS_LP.vs.B.evidence <- RS_LP.vs.B %>% transmute(entrez, FC = avg_log2FC, pvalue = p_val)
write.csv(RS_LP.vs.B.evidence, '../../data/input/CIE/RS_LP.vs.B.evidence.csv')



