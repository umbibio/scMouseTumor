library(openxlsx)
library(tidyverse)
library(seqinr)


GTF.tab <- read.table('~/work/Shailja/Josh/gencode.vM32.annotation.gtf',  sep = '\t',quote = '', row.names = NULL)

GTF.tab.trans <- GTF.tab %>% dplyr::filter(V3 == 'transcript')

GTF.tab.trans$gene_id <- gsub(" ", "", gsub("\"", "", gsub('gene_id', "", 
                                                           unlist(lapply(strsplit(GTF.tab.trans$V9, split = ';'), `[[`, 1)))))
GTF.tab.trans$transcript_id <- gsub(" ", "", gsub("\"", "", gsub('transcript_id', "", 
                                                                 unlist(lapply(strsplit(GTF.tab.trans$V9, split = ';'), `[[`, 2)))))

GTF.tab.trans$gene_name <- gsub(" ", "", gsub("\"", "", gsub('gene_name', "", 
                                                         unlist(lapply(strsplit(GTF.tab.trans$V9, split = ';'), `[[`, 4)))))
GTF.tab.trans$gene_type <- gsub(" ", "", gsub("\"", "", gsub('gene_type', "", 
                                                             unlist(lapply(strsplit(GTF.tab.trans$V9, split = ';'), `[[`, 3)))))

GTF.tab.trans <- GTF.tab.trans %>% mutate(up.strt = ifelse(V7 == '+', V4 - 1500, V5 + 1),
                                          up.end = ifelse(V7 == '+', V4 - 1, V5 + 1500))


promoter <- GTF.tab.trans %>% transmute(chr = V1, strt = up.strt, stp = up.end, 
                                        name = '.', score = '.', strand = V7, gene_id = gene_id, transcript_id = transcript_id,
                                        gene_name = gene_name, gene_type = gene_type)

promoter$stp <- as.numeric(format(promoter$stp, scientific=F))
promoter$strt <- as.numeric(format(promoter$strt, scientific=F))

promoter <- promoter %>% dplyr::filter(strt < stp & strt > 0 & stp > 0)
write.table(promoter, '~/work/Shailja/Josh/promoter.bed', quote = F, sep = '\t', row.names = F, col.names = F)


## Filter for gene clusters
long.clust.tab <- read.xlsx('../Input/Ap2Review/Josh/out/cluster_sig_genes.xlsx')

length(long.clust.tab$GeneID)
sum(long.clust.tab$GeneID %in% promoter$gene_name)

promoter.clust <- inner_join(promoter, long.clust.tab, by = c('gene_name' = 'GeneID')) 
write.table(promoter.clust, '~/work/Shailja/Josh/promoter_clust.bed', quote = F, sep = '\t', row.names = F, col.names = F)

for(i in unique(promoter.clust$cluster)){
  tmp <- promoter.clust %>% dplyr::filter(cluster == i)
  out.file <- paste0('~/work/Shailja/Josh/promoter_clust-', i, '.bed')
  write.table(tmp, out.file, quote = F, sep = '\t', row.names = F, col.names = F)
  
}

GTF.tab.trans[grep('prmt5', GTF.tab.trans$gene_name, ignore.case = T),]
# chr.names <- read.table('~/work/Shailja/Josh/names.txt', sep = ' ',quote = '', row.names = NULL, fill = T)
# chr.names <- chr.names %>% transmute(ID = V1, Chr = V4)
# chr.names$ID <- gsub('>', '', chr.names$ID)
# chr.names$Chr <- gsub(',', '', chr.names$Chr)
# chr.names <- chr.names[-1, ] %>% group_by(ID) %>% summarise(Chr.name = Chr[1])
# chr.names <- chr.names[-grep('NT', chr.names$ID),]
# chr.names <- chr.names[-grep('NW', chr.names$ID),]
# fa <- read.fasta('~/work/Shailja/Josh/GCF_000001635.27_GRCm39_genomic.fna')
# xx <- fa[which(names(fa) %in% chr.names$ID)]
# names(xx) <- paste('chr', chr.names$Chr.name, sep = '')
# write.fasta(xx, names = names(xx), file.out = '~/work/Shailja/Josh/Clean_GRCm39_genomic.fa')
