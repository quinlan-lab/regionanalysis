args = commandArgs(trailingOnly = TRUE)
library(dplyr) 
library('reader')
#ccrs <- read.table('/Users/arq5x/Google Drive/CCR/Data/gnomadbased-ccrs.bed.gz', header=FALSE)
ccrs <- read.table(args[1], header=FALSE)
colnames(ccrs) <- c("chrom","start", "end", "gene", "transcript", "exon", "ranges", "varflag", "syn_density", "cpg", "cov_score", "resid", "resid_pctile", "ccr_pct")

write.table(ccrs %>% filter(ccr_pct >= 99) %>% group_by(gene, ranges) %>% summarise() %>% group_by(gene) %>% count() %>% filter(n>1), file=args[2],row.names=FALSE,sep="\t",quote=FALSE)
