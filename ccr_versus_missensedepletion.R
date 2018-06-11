args = commandArgs (trailingOnly = TRUE)
library(dplyr) 
library(hexbin)
library(readr)
ccr_v_mpcreg <- read.table(args[1], header=TRUE)

ccr_v_mpcreg_clean <- ccr_v_mpcreg %>% filter(mpcreg >= 0)
y <- ccr_v_mpcreg_clean %>% filter(ccr >= as.numeric(args[2]))
pdf(file = paste(args[3], "/ccr", args[2], "_v_mpc.pdf", sep=""))

breaks <- pretty(y$mpcreg,40)
hist(y$mpcreg, breaks=breaks, xlab="Missense depletion", ylab=paste("Number of CCRs at or above the ", args[2], "th percentile", sep=""), col=ifelse(breaks >= 0.4,'lightblue','lightgrey'), main="", las=1)
abline(v = 0.4, col="black", lwd=3, lty=2)

dev.off()

# constrained regions that are unique to CCR
z <- y %>% filter(mpcreg > 0.4)
write.table(z, paste(args[3], "/ccr", args[2], "vmpc(supp_table_4).tsv", sep=""), sep="\t", row.names=F, quote=F)
