dat <- read.csv("C:/Users/rose/Dropbox/current/root research/docs 2014/DAMM_MCNiP/decompvbioc_MCNiP.csv", header=TRUE, sep=",", na.strings="NA", dec=".", strip.white=TRUE)
dat$Root.DOCN <- as.factor (dat$Root.DOCN) 
                                   
library(ggplot2)
ggplot(dat, aes(x=ex, y=DB, color=Root.DOCN)) + geom_point() + geom_line() + ylab("SOM decomposition (mg cm-3 hr-1) / Microbial biomass (mg cm-3)") + xlab("Exudate input (mg cm-3 hr-1)") #+
#theme_bw() +
# theme(panel.grid.major=element_blank(),panel.grid.minor = element_blank(), plot.background = element_rect(fill="white", colour = NA,), legend.position = c(1,1), legend.justification=c(1,1)) + theme(legend.background=element_rect(fill="white", colour="black")) 

ggplot(dat, aes(x=ex, y=SOC, color=Root.DOCN)) + geom_point() + geom_line() + ylab("SOC pool (mg cm-3)") + xlab("Exudate input (mg cm-3 hr-1)")
