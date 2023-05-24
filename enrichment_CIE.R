library(openxlsx)
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(fdrtool)



## ATAC
in.dir <- '../../data/output/CIE/'
all.files <- list.files(in.dir)

all.clust.items <- list()
for(f in all.files){
  contrast <- gsub('TernarytcChIP', '', gsub('\\.evidence\\.csv\\.tsv', '', f))
  tmp <- read_tsv(paste(in.dir, f, sep = ''))
  tmp <- tmp %>% dplyr::filter(adj.pval.down < 0.05 | adj.pval.up < 0.05)
  tmp$contrast <- contrast
  tmp$dir <- ifelse(tmp$adj.pval.up <  tmp$adj.pval.down, 'up', 'down')
  tmp$qval <- ifelse(tmp$dir == 'up', tmp$adj.pval.up, tmp$adj.pval.down)
  tmp <- tmp %>% transmute(name, dir, qval, contrast)
  all.clust.items <- c(all.clust.items, list(tmp))
}


all.clust.items <- bind_rows(all.clust.items)


all.clust.items$name <- factor(all.clust.items$name, 
                               levels = unique(all.clust.items$name[order(all.clust.items$qval)]))

write.xlsx(all.clust.items, '../../data/output/CIE/all_active_TFs.xlsx')

## Category of contrasts
p <- ggplot(all.clust.items, aes(x = name, y = contrast)) + 
  geom_point(aes(colour = dir, size = -log(qval))) +
  theme_bw(base_size = 14) +
  #scale_colour_gradient(limits=c(0, 0.01), low="red") +
  ylab(NULL) + xlab(NULL) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6, face="bold")) + 
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10, face="bold")) +
  theme(legend.position="none") +
  theme(strip.background = element_rect(colour="black", fill="white", 
                                        size=0.5, linetype="solid")) + guides(color = none)+
  ggtitle(strsplit(in.dir,split = "/")[[1]][5]) +
  theme(plot.title = element_text(size = 10))

plot(p)


ggsave(filename="../../data/output/figs/Ternary_active_TFs.pdf", 
       plot=p,
       width = 18, height = 4, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)

