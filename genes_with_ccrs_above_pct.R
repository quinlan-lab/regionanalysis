library(dplyr) 
ccrs <- read.table('/uufs/chpc.utah.edu/common/home/u1021864/analysis/essentials/gnomadbased-ccrs.bed.gz', header=FALSE)
colnames(ccrs) <- c("chrom","start", "end", "gene", "transcript", "exon", "ranges", "varflag", "syn_density", "cpg", "cov_score", "resid", "resid_pctile", "ccr_pct")

ccr_gene_50 <- dim(ccrs %>% filter(ccr_pct >= 50) %>% select(gene) %>% unique())[1]
ccr_gene_80 <- dim(ccrs %>% filter(ccr_pct >= 80) %>% select(gene) %>% unique())[1]
ccr_gene_90 <- dim(ccrs %>% filter(ccr_pct >= 90) %>% select(gene) %>% unique())[1]
ccr_gene_95 <- dim(ccrs %>% filter(ccr_pct >= 95) %>% select(gene) %>% unique())[1]
ccr_gene_99 <- dim(ccrs %>% filter(ccr_pct >= 99) %>% select(gene) %>% unique())[1]
ccr_gene_995 <- dim(ccrs %>% filter(ccr_pct >= 99.5) %>% select(gene) %>% unique())[1]
ccr_gene_999 <- dim(ccrs %>% filter(ccr_pct >= 99.9) %>% select(gene) %>% unique())[1]

pdf(file = "/uufs/chpc.utah.edu/common/home/u1021864/public_html/randomplots/genes_with_ccrs_above_pct.pdf")

bp = barplot(c(ccr_gene_50, ccr_gene_80, ccr_gene_90, ccr_gene_95, ccr_gene_99, ccr_gene_995, ccr_gene_999), 
        names.arg = c("50", "80", "90", "95", "99", "99.5", "99.9"), 
        xlab = "Minimum CCR percentile",
        ylab = "Number of genes with >= CCR",col = "dodgerblue",
        ylim=c(0,18500))

text(x = bp, y = c(ccr_gene_50, ccr_gene_80, ccr_gene_90, ccr_gene_95, ccr_gene_99, ccr_gene_995, ccr_gene_999), 
     label = c(ccr_gene_50, ccr_gene_80, ccr_gene_90, ccr_gene_95, ccr_gene_99, ccr_gene_995, ccr_gene_999), 
     pos = 3, cex = 0.8, col = "black")

dev.off()
