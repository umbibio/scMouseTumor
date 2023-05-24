library(ggplot2)
library(openxlsx)
library(tidytext)
library(ggrepel)
library(dtwclust)
library(geomtextpath)

## After Running the script on Watson
## nohup Rscript Josh_lp2tumorTF_Rscript.R ./lognormCounts.rds ./pseudotimeData.rds ./out/spline_fits.rds


##
spline.fits <- readRDS('../../data/input/rds/spline_fits.rds')
trans.time <- readRDS( '../../data/input/rds/trans_time.rds')

length(unique(spline.fits$GeneID))


## Combining all
RS_LP.vs.LP$contrast <- 'RS-LP.vs.LP'
Tumor.vs.RS_LP$contrast <- 'Tumor.vs.RS-LP'
RS_LP.vs.B$contrast <- 'RS-LP.vs.B'

marker.genes <- bind_rows(RS_LP.vs.LP, Tumor.vs.RS_LP, RS_LP.vs.B)


## Turn the data into wide format (time by gene) and center & scale each gene
sc.rna.dtw.wide <- spline.fits %>% 
  pivot_wider(names_from = 'GeneID', values_from = 'y') %>% 
  mutate_at(vars(-1), scale) %>%
  as.data.frame()


sc.rna.mu.scale <- sc.rna.dtw.wide %>% 
  pivot_longer(-x, names_to = 'GeneID', values_to = 'expr')

sc.rna.mu.scale <- sc.rna.mu.scale[!is.na(sc.rna.mu.scale$expr), ]

## Look at Marker genes only (no down genes)
#sc.rna.mu.scale <- sc.rna.mu.scale %>% dplyr::filter(GeneID %in% marker.genes$GeneID)
## GGplot cluster graphs
my.TFs <- c('Mybl2',
            'Foxm1',
            'Lmnb1',
            'E2f1',
            'E2f7',
            'Ezh2',
            'Fosl1',
            'Tp63',
            'Znf750')


## Clustering IMCs
active.TFs <- read.xlsx('../../data/output/CIE/all_active_TFs.xlsx')
mourse.TFs <- tolower(active.TFs$name)
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
mourse.TFs <- firstup(mourse.TFs)

k <- 3
clust.mat <- sc.rna.dtw.wide[,colnames(sc.rna.dtw.wide) %in%  unique(sc.rna.mu.scale$GeneID)]
clust.mat.Tfs <- sc.rna.dtw.wide[,colnames(sc.rna.dtw.wide) %in%  unique(mourse.TFs)]
clust.mat.Tfs <- sc.rna.dtw.wide[,colnames(sc.rna.dtw.wide) %in%  my.TFs]

hc_dtw <- dtwClustCurves(clust.mat, nclust = k)
hc_dtw.TFs <- dtwClustCurves(clust.mat.Tfs, nclust = k)

plot(hc_dtw, type = 'sc')
plot(hc_dtw.TFs, type = 'sc')




## GGplot cluster graphs
dat.long <- sc.rna.mu.scale %>% dplyr::filter(GeneID %in% my.TFs)
dat.long <- dat.long %>% 
  transmute(time = x, GeneID = GeneID, normExpr = expr)

clust.info <- data.frame(GeneID = colnames(clust.mat.Tfs), cluster = cutree(hc_dtw.TFs, k = k))

long.clust <- left_join(dat.long, clust.info, by = 'GeneID')
long.clust.tab <- long.clust %>% dplyr::select(GeneID, cluster) %>% distinct()
#write.xlsx(long.clust.tab, '../../data/output/tables/cluster_sig_genes.xlsx')

long.clust$label <- NA
long.clust$label[which(long.clust$time == 0.5)] <- long.clust$GeneID[which(long.clust$time == 0.5)]


long.clust$cluster <- paste('C', long.clust$cluster, sep = '')
long.clust$cluster <- factor(long.clust$cluster, levels = unique(sort(long.clust$cluster)))

p  <- ggplot(long.clust, aes(x= time,y=normExpr)) +
  geom_path(aes(color = GeneID),alpha = 0.8, size = 0.8)+ 
  theme_bw(base_size = 16) +
  theme(legend.position = "right") +
  ylab('normExpr') + xlab('Time') + #+ ylim(c(0, 5)) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 18, face="bold")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 18, face="bold")) +
  theme(strip.background = element_rect(colour="black", fill="gray90",
                                        size=0.5, linetype="solid")) +
  theme(strip.text = element_text(size = 18, face="bold", angle = 0)) + 
  
  coord_cartesian(xlim = c(0,1.1)) + 
  geom_text_repel(aes(label = label), size = 4, fontface = "bold",
                  box.padding = unit(0.6, "lines"),
                  max.overlaps = 300,
                  #segment.angle = 180,
                  nudge_x = 0.25, 
                  nudge_y = 0.25,
                  hjust=0.25,
                  #nudge_x=0.25, 
                  segment.size = 0.1,
                  na.rm = TRUE)+ 
  facet_wrap(.~cluster, scales = 'free') +
  
  
  #ggtitle(titles[i]) +
  theme(
    plot.title = element_text(size=18, face = "bold.italic", color = 'red'),
    axis.title.x = element_text(size=18, face="bold", hjust = 1),
    axis.title.y = element_text(size=18, face="bold")
  ) + 
  theme(#legend.position = c(0.15, 0.85),
    legend.position = 'none',
    legend.title = element_text(colour="black", size=12, 
                                face="bold"),
    legend.text = element_text(colour="black", size=12, 
                               face="bold"))






plot(p)

ggsave(filename="../Output/AP2sReview/figs/ap2_clusteres.pdf",
       plot=p,
       width = 14, height = 10,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)


