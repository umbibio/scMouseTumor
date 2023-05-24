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
n.time <- length(trans.time$qq) 
trans.time$strt <- c(0,(trans.time$qq[1:(n.time - 1)] + trans.time$qq[2:(n.time)]) / 2)
trans.time$stp <- c((trans.time$qq[1:(n.time - 1)] + trans.time$qq[2:(n.time)]) / 2, 1)

saveRDS(trans.time, '../../data/input/rds/trans_time.rds')
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
saveRDS(spline.fits, '../../data/input/rds/spline_fits.rds')
