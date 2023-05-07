library(tidyverse)
library(openxlsx)
library(doParallel)
library(npreg)

xx <- readRDS('../../data/input/fromJosh/lognormCounts.rds')
yy <- readRDS('../../data/input/fromJosh/pseudotimeData.rds')

yy <- data.frame(cells = names(yy), pt = c(yy))

xx <- as.data.frame(as.matrix(xx))
xx$GeneName <- rownames(xx)
xx.long <- xx %>% pivot_longer(-GeneName, names_to = 'cells', values_to = 'expr')


xx.long <- left_join(xx.long, yy, by = 'cells')
xx.long$st <- (xx.long$pt - min(xx.long$pt))/ (max(xx.long$pt) - min(xx.long$pt))

cell.idents <- readRDS('../../data/input/rds/cellIdents.rds')
xx.long <- left_join(xx.long, cell.idents, by = c('cells' = 'cellID'))
xx.long$cellType <- factor(xx.long$cellType, levels = c("LP", "RS-LP", "RS-T", "tumor"))
cell.to.time <- xx.long %>% dplyr::select(st, cellType) %>% distinct()
plot(cell.to.time$cellType, cell.to.time$st)

## Identify transition points
trans.time <- cell.to.time %>% group_by(cellType) %>% summarise(qq = quantile(st, probs = 0.75))
n.time <- length(tans.time$qq) 
trans.time$strt <- c(0,(tans.time$qq[1:(n.time - 1)] + tans.time$qq[2:(n.time)]) / 2)
trans.time$stp <- c((tans.time$qq[1:(n.time - 1)] + tans.time$qq[2:(n.time)]) / 2, 1)

## For parallel calculations

num.cores <- 16


genes <- unique(xx.long$GeneName)
i = which(genes == 'S100a1')
i = which(genes == 'Krt5')

lbx <- 2.0

## Expression graphs

spline.fits <- mclapply(1:length(genes), function(i){
  tmp <- xx.long %>% dplyr::filter(GeneName == genes[i]) %>%
    transmute(GeneID = GeneName, x = st, y = expr)
  
  
  y <- tmp$y
  t <- tmp$x
  
  rna.sp <- smooth.spline(t, y, lambda = 0.1)
  
  # w <- rep(1, length(y))
  # w[which(y == 0)] <- 1/3
  # 
  # rna.sp <- smooth.spline(t, y, cv = T, w = w)
  # sparx <- rna.sp$spar
  # 
  # if(sparx < lbx){
  #   
  #   rna.sp <- smooth.spline(t, y, spar = lbx, w = w)
  #   
  # }
  
  rna.sp <- predict(rna.sp, seq(0, 1, by = 0.01))
  
  ######
  #plot(t, y, col = 'red')
  #points(rna.sp$x, rna.sp$y, type = 'l', col = 'green', lwd = 2)
  mu <- data.frame(x = rna.sp$x, y = rna.sp$y)
  
  mu <- data.frame(GeneID = rep(tmp$GeneID[1], length(mu[,1])), x = mu[,1], y = mu[,2])
  return(mu)
}, mc.cores = num.cores)

spline.fits <- bind_rows(spline.fits)
saveRDS(spline.fits, '../../data/output/scMouseTumor/spline_fits.rds')


## Plot specific splines
spline.fits <- readRDS('../../data/output/scMouseTumor/spline_fits.rds')

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
    geom_line(color = 'blue',alpha = 0.8, size = 0.8)+ ylim(c(min(my.sc.rna$expr) - 0.1, max(my.sc.rna$expr) + 0.1)) +
    geom_line(aes(x= x,y=min(expr), color = cellType), size = 2) + 
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
my.GeneID <- map.id('PHF2')[2]
my.GeneID <- map.id('FOSL1')
my.GeneID <- map.id('RUNX1')
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

