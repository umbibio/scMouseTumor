library(Seurat)
library(ggplot2)
library(ggfortify)
library(jcolors)
require(gridExtra)
library(grid)
library(gam)
library(princurve)
library(parallel)
library(tidyverse)
library(MyEllipsefit)
library(sctransform)
library(openxlsx)
library(doParallel)
library(tidytext)
library(ggrepel)
library(dtwclust)
library(geomtextpath)


source('./util_funcs.R')

## Plot specific splines
spline.fits <- readRDS('../../data/input/rds/spline_fits.rds')
trans.time <- readRDS( '../../data/input/rds/trans_time.rds')


## Turn the data into wide format (time by gene) and center & scale each gene
spline.fits.wide <- spline.fits %>% 
  pivot_wider(names_from = 'GeneID', values_from = 'y') %>% 
  mutate_at(vars(-x), ~scale(., center = T, scale = T)) %>%
  as.data.frame()


spline.fits.long <- spline.fits.wide  %>% 
  pivot_longer(-x, names_to = 'GeneID', values_to = 'expr')


plot_trends <- function(my.GeneID, trans.time){
  
  my.sc.rna <- spline.fits.long %>% dplyr::filter(GeneID == my.GeneID)
  trans.time <- trans.time %>% arrange(cellType)
  orders <- levels(trans.time$cellType)
  cellType <- as.character(trans.time$cellType)
  my.sc.rna <- my.sc.rna %>% 
    mutate(cellType = ifelse(x >= trans.time$strt[1] & x < trans.time$stp[1], (cellType[1]),
                             ifelse(x >= trans.time$strt[2] & x < trans.time$stp[2], cellType[2],
                                    ifelse(x >= trans.time$strt[3] & x < trans.time$stp[3], cellType[3],
                                           cellType[4]))))
  
  my.sc.rna$cellType <- factor(my.sc.rna$cellType, levels = orders)
  p  <- ggplot(my.sc.rna, aes(x= x,y=expr)) +
    geom_line(color = 'blue',alpha = 0.8, linewidth = 0.8)+ ylim(c(min(my.sc.rna$expr) - 0.1, max(my.sc.rna$expr) + 0.1)) +
    geom_line(aes(x= x,y=min(expr), color = cellType), linewidth = 2) + 
    theme_bw(base_size = 14) + #ylim(c(min(my.sc.rna$expr) - 0.1, max(my.sc.rna$expr) + 0.1)) +
    theme(legend.position = "right") +
    ylab('normExpr') + xlab('Time') +
    theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
    theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
    theme(strip.background = element_rect(colour="black", fill="white",
                                          size=0.5, linetype="solid")) +
    ggtitle(paste('expression curve', my.GeneID)) +
    theme(strip.text = element_text(size = 14, face="bold", angle = 0)) +
    theme(
      plot.title = element_text(size=14, face = "bold.italic", color = 'black'),
      axis.title.x = element_text(size=14, face="bold", hjust = 1),
      axis.title.y = element_text(size=14, face="bold")
    )
  
  
  return(p)
}


map.id <- function(gene.name){
  my.GeneID <- unique(spline.fits.long$GeneID[grep(gene.name, spline.fits.long$GeneID, ignore.case = T)])
}

my.GeneID <- map.id('foxm1')
my.GeneID <- map.id('EZH2')
my.GeneID <- map.id('FOSL1')
my.GeneID <- map.id('LMNB1')
my.GeneID <- map.id('TEAD4')
my.GeneID <- map.id('BRD4')
my.GeneID <- map.id('PSMB5')
my.GeneID <- map.id('HIF1A')[1]
my.GeneID <- map.id('Trp63')[1]

p <- plot_trends(my.GeneID, trans.time)
plot(p)

my.genes <- c(map.id('foxm1'), map.id('FOSL1'), map.id('RUNX1'), map.id('prmt5'))

for(i in 1:length(my.genes)){
  p <- plot_trends(my.genes[i], trans.time)
  f.n <- paste("../../data/output/figs/", my.genes[i], '.pdf', sep = '')
  ggsave(filename=f.n,
         plot=p,
         width = 6, height = 6,
         units = "in", # other options are "in", "cm", "mm"
         dpi = 300
  )
  
  
}

