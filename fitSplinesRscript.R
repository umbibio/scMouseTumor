library(tidyverse)
library(openxlsx)
library(doParallel)
library(npreg)

options <- commandArgs(trailingOnly = TRUE)

in.file1 <- options[1]
in.file2 <- options[2]
out.file <- options[3]

xx <- readRDS(in.file1)
yy <- readRDS(in.file2)

#xx <- readRDS('../Input/Ap2Review/Josh/lp2tumor_files/lognormCounts.rds')
#yy <- readRDS('../Input/Ap2Review/Josh/lp2tumor_files/pseudotimeData.rds')

yy <- data.frame(cells = names(yy), pt = c(yy))

xx <- as.data.frame(xx)
xx$GeneName <- rownames(xx)
xx.long <- xx %>% pivot_longer(-GeneName, names_to = 'cells', values_to = 'expr')


xx.long <- left_join(xx.long, yy, by = 'cells')
xx.long$st <- (xx.long$pt - min(xx.long$pt))/ (max(xx.long$pt) - min(xx.long$pt))


#xx.long.traj1 <- xx.long %>% dplyr::filter(ident %in% c('LP', 'RS-LP'))

## For parallel calculations

num.cores <- detectCores(all.tests = FALSE, logical = TRUE)


genes <- unique(xx.long$GeneName)
#i = which(genes == 'S100a1')
#i = which(genes == 'Krt5')

lbx <- 1.7

## Expression graphs

spline.fits <- mclapply(1:10, function(i){
  tmp <- xx.long %>% dplyr::filter(GeneName == genes[i]) %>%
    transmute(GeneID = GeneName, x = st, y = expr)
  
  
  y <- tmp$y
  t <- tmp$x
  
  
  w <- rep(1, length(y))
  w[which(y == 0)] <- 1/2
  
  rna.sp <- smooth.spline(t, y, cv = T, w = w)
  sparx <- rna.sp$spar
  
  if(sparx < lbx){
    
    rna.sp <- smooth.spline(t, y, spar = lbx, w = w)
    
  }
  
  rna.sp <- predict(rna.sp, seq(0, 1, by = 0.01))
  
  ######
  #plot(t, y, col = 'red')
  #points(rna.sp$x, rna.sp$y, type = 'l', col = 'green', lwd = 2)
  mu <- data.frame(x = rna.sp$x, y = rna.sp$y)
  
  mu <- data.frame(GeneID = rep(tmp$GeneID[1], length(mu[,1])), x = mu[,1], y = mu[,2])
  return(mu)
}, mc.cores = num.cores)

spline.fits <- bind_rows(spline.fits)
saveRDS(spline.fits, out.file)


