library(plyr)
library(ggplot2)
library(scales)

#putting all the data together in order to merge
#common colnames: date (as.date), tree

#import all data####
#this dataset comes from 'exudates_raw_R' ##gC/groot/day
exudates <- read.csv("C:/Users/rose/Dropbox/current/root research/data 2014/Ch2/exudates_gcgr.csv", header=TRUE, sep=",", na.strings="NA", dec=".", strip.white=TRUE)
exudates <- exudates[-c(1,5:9)]
exudates <- ddply(exudates, c("species","date"), summarise, N=length(gcgr), mean = mean(gcgr), sd=sd(gcgr), se=sd/sqrt(N))
colnames(exudates) <- c("tree", "date", "N", "exudates.gcgr", "sd.exudates", "se.exudates")

# #this dataset comes from 'biomass2012' ##fine roots
 biomass.live <- read.csv("C:/Users/rose/Dropbox/current/root research/data 2014/Ch2/live_biomass_plot_means.csv", header=TRUE, sep=",", na.strings="NA", dec=".", strip.white=TRUE) 
 colnames(biomass.live) <- c("X", "date", "tree", "o.m", "N", "biomass.live.gm2", "sd.biomass.live", "se.biomass.live")
 biomass.live <- biomass.live[-1]

#calculate MCNiP seasonal model inputs####
mean.ex.gcgr <- mean(exudates$exudates.gcgr)
mean.biomass.gm2 <- mean(c(304.7527, 241.27, 190.0082))
ex.gCm2day <- mean.ex.gcgr*mean.biomass.gm2
ex.mgcm3hr <- (mean.ex.gcgr*mean.biomass.gm2)*1000/(0.15*100*100*100*24)
ex.gCm2yr <- ex.gCm2day*30.4*12
turnover <- 5.5 #years turnover averaging all years for each stand and then averaging the averages using all minirhiz/Burton method #3.134 years turnover time average of two years oak, hem, 1 year ash using Burton 2000 method from Abramoff_DOE_342014
ttime <- 1/mean(turnover)
avail <- ttime*mean.biomass.gm2*0.4 #g C (xC:N ratio)/m2/yr available root litter for microbes
av.mgcm3hr <-  avail*1000/(0.15*100*100*100*24*365)

yearex <- ex.mgcm3hr*24*365
yearturn <- av.mgcm3hr*24*365

ex.turn <- ex.mgcm3hr + av.mgcm3hr

#calculate MCNiP seasonal model inputs PER TREE, (FRAM, QURU, TSCA)####
tree.ex.gcgr <- ddply(exudates, c("tree"), summarise, N=length(exudates.gcgr), mean = mean(exudates.gcgr, na.rm=TRUE), sd=sd(exudates.gcgr, na.rm=TRUE), se=sd/sqrt(N))
tree.biomass.gm2 <- c( 190.0082, 304.7527, 241.27)

tree.ex.gCm2day <- tree.ex.gcgr$mean*tree.biomass.gm2

tree.ex.mgcm3hr <- (tree.ex.gCm2day)*1000/(0.15*100*100*100*24)
tree.ex.gCm2yr <- tree.ex.gCm2day*30.4*12

tree.turnover <- c(1.46, 2.39, 12.95) #FRAM, QURU, TSCA
tree.time <- 1/tree.turnover
tree.avail <- tree.time*tree.biomass.gm2$mean*0.4 #g C (xC:N ratio)/m2/yr available root litter for microbes

tree.av.mgcm3hr <- tree.avail*1000/(0.15*100*100*100*24*365)

tree.ex.mgcm3hr  
tree.av.mgcm3hr
tree.boz <- tree.ex.mgcm3hr + tree.av.mgcm3hr

#calc vmax####
df <-read.csv("C:/Users/rose/Dropbox/current/root research/data 2014/Ch4/mydataFlux.csv", sep=",",header=T, na.strings="NA")

#uptake and dep have identical eqns and parameters
a <- 1.0815e11
ea <- 61.77
R <- 0.008314
vmax <- a* exp(-ea/(R*(df$SoilT+273)))
vmax.w <- a* exp(-ea/(R*(df$SoilT+273+5)))

sum(vmax)
sum(vmax.w)
