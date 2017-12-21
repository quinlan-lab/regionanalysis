args = commandArgs (trailingOnly = TRUE)
library(dplyr) 
library(hexbin)
library(readr)
ccr_v_mpcreg <- read.table(args[1], header=TRUE)

pdf(file = args[2])
ccr_v_mpcreg_clean <- ccr_v_mpcreg %>% filter(mpcreg >= 0)
y <- ccr_v_mpcreg_clean %>% filter(ccr >= 99)

breaks <- pretty(y$mpcreg,40)
hist(y$mpcreg, breaks=breaks, xlab="Missense depletion", ylab="Number of CCRs at or above the 99th percentile", col=ifelse(breaks >= 0.4,'lightblue','lightgrey'), main="", las=1)
abline(v = 0.4, col="black", lwd=3, lty=2)

dev.off()

# constrained regions that are unique to CCR
z <- y %>% filter(mpcreg > 0.4)
write.table(z, args[3], sep="\t", row.names=F, quote=F)
