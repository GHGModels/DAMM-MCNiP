#upload####
df <-read.csv("/Users/rzabramoff/Documents/records/dropbox backup/dropbox.6.1.2015/BU/Dissertation_Research/data 2014/Ch4/mydataFlux.csv", sep=",",header=T, na.strings="NA")
bulkAM.lit <-read.csv("/Users/rzabramoff/Documents/records/dropbox backup/dropbox.6.1.2015/BU/Dissertation_Research/data 2014/Ch4/bulkAMlit_032215.csv", sep=",",header=T, na.strings="NA")
rhizoAM.lit <-read.csv("/Users/rzabramoff/Documents/records/dropbox backup/dropbox.6.1.2015/BU/Dissertation_Research/data 2014/Ch4/rhizoAMlit_032215.csv", sep=",",header=T, na.strings="NA")
bulkEM.lit <-read.csv("/Users/rzabramoff/Documents/records/dropbox backup/dropbox.6.1.2015/BU/Dissertation_Research/data 2014/Ch4/bulkEMlit_032215.csv", sep=",",header=T, na.strings="NA")
rhizoEM.lit <-read.csv("/Users/rzabramoff/Documents/records/dropbox backup/dropbox.6.1.2015/BU/Dissertation_Research/data 2014/Ch4/rhizoEMlit_032215.csv", sep=",",header=T, na.strings="NA")
bulkAM.mic <-read.csv("/Users/rzabramoff/Documents/records/dropbox backup/dropbox.6.1.2015/BU/Dissertation_Research/data 2014/Ch4/bulkAMmic_032215.csv", sep=",",header=T, na.strings="NA")
rhizoAM.mic <-read.csv("/Users/rzabramoff/Documents/records/dropbox backup/dropbox.6.1.2015/BU/Dissertation_Research/data 2014/Ch4/rhizoAMmic_032215.csv", sep=",",header=T, na.strings="NA")
bulkEM.mic <-read.csv("/Users/rzabramoff/Documents/records/dropbox backup/dropbox.6.1.2015/BU/Dissertation_Research/data 2014/Ch4/bulkEMmic_032215.csv", sep=",",header=T, na.strings="NA")
rhizoEM.mic <-read.csv("/Users/rzabramoff/Documents/records/dropbox backup/dropbox.6.1.2015/BU/Dissertation_Research/data 2014/Ch4/rhizoEMmic_032215.csv", sep=",",header=T, na.strings="NA")
bulkAM <-read.csv("/Users/rzabramoff/Documents/records/dropbox backup/dropbox.6.1.2015/BU/Dissertation_Research/data 2014/Ch4/bulkAM_032215.csv", sep=",",header=T, na.strings="NA")
rhizoAM <-read.csv("/Users/rzabramoff/Documents/records/dropbox backup/dropbox.6.1.2015/BU/Dissertation_Research/data 2014/Ch4/rhizoAM_032215.csv", sep=",",header=T, na.strings="NA")
bulkEM <-read.csv("/Users/rzabramoff/Documents/records/dropbox backup/dropbox.6.1.2015/BU/Dissertation_Research/data 2014/Ch4/bulkEM_032215.csv", sep=",",header=T, na.strings="NA")
rhizoEM <-read.csv("/Users/rzabramoff/Documents/records/dropbox backup/dropbox.6.1.2015/BU/Dissertation_Research/data 2014/Ch4/rhizoEM_032215.csv", sep=",",header=T, na.strings="NA")

#bulk rhizo scale####
#convert mg/cm3/hr to mg/m2/hr
frac = 0.2048506
AM.lit <- NULL
AM.lit <- (rhizoAM.lit*frac+bulkAM.lit*(1-frac))*(0.1*100*100*100)
EM.lit <- NULL
EM.lit <- (rhizoEM.lit*frac+bulkEM.lit*(1-frac))*(0.1*100*100*100)
AM.mic <- NULL
AM.mic <- (rhizoAM.mic*frac+bulkAM.mic*(1-frac))*(0.1*100*100*100)
EM.mic <- NULL
EM.mic <- (rhizoEM.mic*frac+bulkEM.mic*(1-frac))*(0.1*100*100*100)
AM <- NULL
AM <- (rhizoAM*frac+bulkAM*(1-frac))*(0.1*100*100*100)
EM <- NULL
EM <- (rhizoEM*frac+bulkEM*(1-frac))*(0.1*100*100*100)

#plots####
df$timeNA <- ifelse(is.na(df$scale), NA, df$time)

#save plots code 6 panel####
pdf(("C:/Users/rose/Dropbox/current/root research/data 2014/Ch4/micCN_031615.pdf"))

#shrinkdis####
par(mfrow=c(2,3))
par(mar = c(5,5,5,1))
par(oma = c(0,0,0,0))
#Litter C min
plot(df$timeNA, AM.lit$CMIN, type="l", ylim = c(0,200), main = "Input C:N", col =6, cex = 0.5, cex.axis = 1.3, cex.lab = 1.5, xlab = "Day of year", ylab = expression("C flux (mg " ~ m^{-2} ~ hr^{-1} ~ ")"))
points(df$timeNA, EM.lit$CMIN, col = 4, cex = 0.5, type="l")
legend("topleft",c("Narrow C:N","Wide C:N"),
       pch = c(16,16),
       col=c(6,4),cex=1.3, bty = "n"
)
text(75,230, "a)", cex = 1.5, xpd = NA, font = 2)

#Microbial C min
plot(df$timeNA, EM.mic$CMIN, type="l", ylim = c(0,200), main = "Microbial biomass C:N", col =4, cex = 0.5, cex.axis = 1.3, cex.lab = 1.5, xlab = "Day of year", ylab = expression("C flux (mg " ~ m^{-2} ~ hr^{-1} ~ ")"))
points(df$timeNA, AM.mic$CMIN, col = 6, cex = 0.5, type="l")
text(75,230, "b)", cex = 1.5, xpd = NA, font = 2)

#C min
plot(df$timeNA, AM$CMIN, type="l", ylim = c(0,200), main = "All C:N", col =6, cex = 0.5, cex.axis = 1.3, cex.lab = 1.5, xlab = "Day of year", ylab = expression("C flux (mg " ~ m^{-2} ~ hr^{-1} ~ ")"))
points(df$timeNA, EM$CMIN, col = 4, cex = 0.5, type="l")
text(75,230, "c)", cex = 1.5, xpd = NA, font = 2)

#Litter N min
plot(df$timeNA, AM.lit$NMIN, type="l", ylim = c(0,4), col =6, cex = 0.5, cex.axis = 1.3, cex.lab = 1.5, xlab = "Day of year", ylab = expression("N flux (mg " ~ m^{-2} ~ hr^{-1} ~ ")"))
points(df$timeNA, EM.lit$NMIN, col = 4, cex = 0.5, type="l")
text(75,4.5, "d)", cex = 1.5, xpd = NA, font = 2)

#Microbial Nmin
plot(df$timeNA, AM.mic$NMIN, type="l", ylim= c(0,4), col =6, cex = 0.5, cex.axis = 1.3, cex.lab = 1.5, xlab = "Day of year", ylab = expression("N flux (mg " ~ m^{-2} ~ hr^{-1} ~ ")"))
points(df$timeNA, EM.mic$NMIN, col = 4, cex = 0.5, type="l")
text(75,4.5, "e)", cex = 1.5, xpd = NA, font = 2)

#N min
plot(df$timeNA, AM$NMIN, type="l", ylim = c(0,4), col =6, cex = 0.5, cex.axis = 1.3, cex.lab = 1.5, xlab = "Day of year", ylab = expression("N flux (mg " ~ m^{-2} ~ hr^{-1} ~ ")"))
points(df$timeNA, EM$NMIN, col = 4, cex = 0.5, type="l")
text(75,4.5, "f)", cex = 1.5, xpd = NA, font = 2)
#shrinked####

dev.off()

#scale
a<- sum(AM$CMIN)*0.001 #mg/m2/hr to g/m2/yr
b<- sum(EM$CMIN)*0.001 
c<- sum(AM$NMIN)*0.001 
d<- sum(EM$NMIN)*0.001 

(a-b)/a
(c-d)/c

a<- sum(AM.lit$CMIN)*0.001 #mg/m2/hr to g/m2/yr
b<- sum(EM.lit$CMIN)*0.001 
c<- sum(AM.lit$NMIN)*0.001 
d<- sum(EM.lit$NMIN)*0.001 

a<- sum(AM.mic$CMIN)*0.001 #mg/m2/hr to g/m2/yr
b<- sum(EM.mic$CMIN)*0.001 
c<- sum(AM.mic$NMIN)*0.001 
d<- sum(EM.mic$NMIN)*0.001 