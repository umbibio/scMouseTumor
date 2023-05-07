library(openxlsx)

lp2tumor <- read.xlsx('../../data/input/fromJosh/lp2tumor_Ident.xlsx')
s2 <- read.xlsx('../../data/input/fromJosh/s2_fvb50_epiIdent.xlsx')

cell.idents <- bind_rows(lp2tumor, s2)
cell.idents <- cell.idents %>% distinct()

saveRDS(cell.idents, '../../data/input/rds/cellIdents.rds')
