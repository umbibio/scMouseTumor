library(ggplot2)
library(openxlsx)
library(tidytext)
library(ggrepel)
library(dtwclust)
library(geomtextpath)

## After Running the script on Watson
## nohup Rscript Josh_lp2tumorTF_Rscript.R ./lognormCounts.rds ./pseudotimeData.rds ./out/spline_fits.rds


##
spline.fits <- readRDS('../Input/Ap2Review/Josh/out/spline_fits_lambda_01.rds')
length(unique(spline.fits$GeneID))

gene.set1 <- read.xlsx('../Input/Ap2Review/Josh/lp2tumor_files/upTumorvLP.xlsx')
gene.set2 <- read.xlsx('../Input/Ap2Review/Josh/lp2tumor_files/upRS_LPvLP.xlsx')

gene.set1 <- gene.set1 %>% dplyr::filter(abs(avg_log2FC) > log2(1.5) & p_val_adj < 0.01)
gene.set2 <- gene.set2 %>% dplyr::filter(abs(avg_log2FC) > log2(1.5) & p_val_adj < 0.01)

sig.genes <- unique(c(gene.set1$GeneID, gene.set2$GeneID))

spline.fits.filt <- spline.fits %>% dplyr::filter(GeneID %in% sig.genes)
dtw.wide <- spline.fits.filt %>% 
  pivot_wider(names_from = 'GeneID', values_from = 'y') %>% 
  mutate_at(vars(!matches("^x$")), ~scale(., center = T, scale = T)) %>%
  as.data.frame()


mu.scale <- dtw.wide %>% 
  pivot_longer(-x, names_to = 'GeneID', values_to = 'expr')


dtw.wide.no.scale <- spline.fits.filt %>% 
  pivot_wider(names_from = 'GeneID', values_from = 'y') %>% 
  mutate_at(vars(!matches("^x$")), ~scale(., center = F, scale = F)) %>%
  as.data.frame()


mu.no.scale <- dtw.wide.no.scale %>% 
  pivot_longer(-x, names_to = 'GeneID', values_to = 'expr')

## Clustering IMCs

k <- 6
clust.mat <- dtw.wide[,colnames(dtw.wide) %in%  sig.genes]
hc_dtw <- dtwClustCurves(clust.mat, nclust = k)

plot(hc_dtw, type = 'sc')


## GGplot cluster graphs
dat.long <- mu.no.scale
dat.long <- dat.long %>% 
  transmute(time = x, GeneID = GeneID, normExpr = expr)

clust.info <- data.frame(GeneID = colnames(clust.mat), cluster = cutree(hc_dtw, k = k))

long.clust <- left_join(dat.long, clust.info, by = 'GeneID')
long.clust.tab <- long.clust %>% dplyr::select(GeneID, cluster) %>% distinct()
write.xlsx(long.clust.tab, '../Input/Ap2Review/Josh/out/cluster_sig_genes.xlsx')

long.clust$label <- NA
long.clust$label[which(long.clust$time == 0.5)] <- long.clust$GeneID[which(long.clust$time == 0.5)]

# sc.rna.sc.atac.joint.long$label <- NA
# AP2s.clusts <- sc.rna.sc.atac.joint.long %>% group_by(Name) %>% summarise(cluster.RNA = cluster.RNA[1]) %>% 
#   ungroup() %>% group_by(cluster.RNA) %>% mutate(lab.time = 1:n() %% 6)
# 
# sc.rna.sc.atac.joint.long$lab.time <- floor(sc.rna.sc.atac.joint.long$time) 
# sc.rna.sc.atac.joint.long <- left_join(sc.rna.sc.atac.joint.long, AP2s.clusts, by = c('Name', 'lab.time', 'cluster.RNA'))
# sc.rna.sc.atac.joint.long$label[sc.rna.sc.atac.joint.long$time == sc.rna.sc.atac.joint.long$lab.time] <-
#   sc.rna.sc.atac.joint.long$Name[sc.rna.sc.atac.joint.long$time == sc.rna.sc.atac.joint.long$lab.time]

long.clust$cluster <- paste('C', long.clust$cluster, sep = '')
long.clust$cluster <- factor(long.clust$cluster, levels = unique(sort(long.clust$cluster)))

p  <- ggplot(long.clust, aes(x= time,y=normExpr)) +
  geom_path(aes(color = GeneID),alpha = 0.8, size = 0.8)+ 
  theme_bw(base_size = 16) +
  theme(legend.position = "right") +
  ylab('normExpr') + xlab('Time') + ylim(c(0, 5)) + 
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


