library(tidyverse)
library(openxlsx)


## Entrez IDS
mouse_human <- read.xlsx('../../data/input/mouseGenomics/mouse_human_orth_entrez.xlsx')

DEGs <- read.xlsx('../../data/input/fromJosh/TumorvLP.xlsx')

DEGs <- inner_join(DEGs, mouse_human, by = c('GeneNam' = 'mouseGeneName'))

evidence <- DEGs %>% transmute(entrez, log2FC = avg_log2FC, qval = p_val) %>%
  dplyr::filter(abs(log2FC ) > log2(1.5), qval <= 0.05)

write.csv(evidence, '../../data/input/CIE/Tumor_vs_LP.csv')
