#upload####
#warm is at 5deg, sum and win warm is at 3deg
df <-read.csv("C:/Users/rose/Dropbox/current/root research/data 2014/Ch4/mydataFlux.csv", sep=",",header=T, na.strings="NA")

bulk <-read.csv("C:/Users/rose/Dropbox/current/root research/data 2014/Ch4/bulk_032115.csv", sep=",",header=T, na.strings="NA")
rhizo <-read.csv("C:/Users/rose/Dropbox/current/root research/data 2014/Ch4/rhizo_032115.csv", sep=",",header=T, na.strings="NA")

bulkwarm5 <-read.csv("C:/Users/rose/Dropbox/current/root research/data 2014/Ch4/bulkwarm5_032215.csv", sep=",",header=T, na.strings="NA")
rhizowarm5 <-read.csv("C:/Users/rose/Dropbox/current/root research/data 2014/Ch4/rhizowarm5_032215.csv", sep=",",header=T, na.strings="NA")
bulkcold5 <-read.csv("C:/Users/rose/Dropbox/current/root research/data 2014/Ch4/bulkcold5_032215.csv", sep=",",header=T, na.strings="NA")
rhizocold5 <-read.csv("C:/Users/rose/Dropbox/current/root research/data 2014/Ch4/rhizocold5_032215.csv", sep=",",header=T, na.strings="NA")

bulkwarm1.5x <-read.csv("C:/Users/rose/Dropbox/current/root research/data 2014/Ch4/bulkwarm1.5x_032215.csv", sep=",",header=T, na.strings="NA")
rhizowarm1.5x <-read.csv("C:/Users/rose/Dropbox/current/root research/data 2014/Ch4/rhizowarm1.5x_032215.csv", sep=",",header=T, na.strings="NA")
bulkcold0.5x <-read.csv("C:/Users/rose/Dropbox/current/root research/data 2014/Ch4/bulkwarm0.5x_032215.csv", sep=",",header=T, na.strings="NA")
rhizocold0.5x <-read.csv("C:/Users/rose/Dropbox/current/root research/data 2014/Ch4/rhizowarm0.5x_032215.csv", sep=",",header=T, na.strings="NA")


bulkwet <-read.csv("C:/Users/rose/Dropbox/current/root research/data 2014/Ch4/bulkwet_032215.csv", sep=",",header=T, na.strings="NA")
rhizowet <-read.csv("C:/Users/rose/Dropbox/current/root research/data 2014/Ch4/rhizowet_032215.csv", sep=",",header=T, na.strings="NA")
bulkdry <-read.csv("C:/Users/rose/Dropbox/current/root research/data 2014/Ch4/bulkdry_032115.csv", sep=",",header=T, na.strings="NA")
rhizodry <-read.csv("C:/Users/rose/Dropbox/current/root research/data 2014/Ch4/rhizodry_032115.csv", sep=",",header=T, na.strings="NA")

#bulk rhizo scale####
#convert mg/cm3/hr to mg/m2/hr
frac = 0.2048506
normal <- NULL
normal <- (rhizo*frac+bulk*(1-frac))*(0.1*100*100*100)
warm5 <- NULL
warm5 <- (rhizowarm5*frac+bulkwarm5*(1-frac))*(0.1*100*100*100)
warmx <- NULL
warmx <- (rhizowarm1.5x*frac+bulkwarm1.5x*(1-frac))*(0.1*100*100*100)
cold5 <- NULL
cold5 <- (rhizocold5*frac+bulkcold5*(1-frac))*(0.1*100*100*100)
coldx <- NULL
coldx <- (rhizocold0.5x*frac+bulkcold0.5x*(1-frac))*(0.1*100*100*100)
wet <- NULL
wet <- (rhizowet*frac+bulkwet*(1-frac))*(0.1*100*100*100)
dry <- NULL
dry <- (rhizodry*frac+bulkdry*(1-frac))*(0.1*100*100*100)

#plots####
df$timeNA <- ifelse(is.na(df$scale), NA, df$time)

#save plots code 4 panel CMIN####
pdf(("C:/Users/rose/Dropbox/current/root research/data 2014/Ch4/warmwet_031615.pdf"))

par(mfrow=c(2,2))
par(mar = c(5,5,1,1))
par(oma = c(0,0,0,0))
#warming
plot(df$timeNA, warm5$CMIN, type="l", ylim = c(0,300), col =6, cex = 0.5, cex.axis = 1.3, cex.lab = 1.5, xlab = "Day of year", ylab = expression("C flux (mg " ~ m^{-2} ~ hr^{-1} ~ ")"))
points(df$timeNA, normal$CMIN, col = 1, cex = 0.5, type="l")
points(df$timeNA, cold5$CMIN, col = 4, cex = 0.5, type="l")
legend("topleft",c("Ambient","+5", "-5"),
       pch = c(16,16,16),
       col=c(1,6,4),cex=1.3, bty = "n"
)

#winter warming
plot(df$timeNA, winwarm$CMIN, type="l", ylim = c(0,300), col =6, cex = 0.5, cex.axis = 1.3, cex.lab = 1.5, xlab = "Day of year", ylab = expression("C flux (mg " ~ m^{-2} ~ hr^{-1} ~ ")"))
points(df$timeNA, normal$CMIN, col = 4, cex = 0.5, type="l")
legend("topleft",c("Ambient","Winter warming"),
       pch = c(16,16),
       col=c(4,6),cex=1.3, bty = "n"
)

#0.5x soilM
plot(df$timeNA, dry$CMIN, type="l", ylim = c(0,300), col =6, cex = 0.5, cex.axis = 1.3, cex.lab = 1.5, xlab = "Day of year", ylab = expression("C flux (mg " ~ m^{-2} ~ hr^{-1} ~ ")"))
points(df$timeNA, normal$CMIN, col = 4, cex = 0.5, type="l")
legend("topleft",c("Ambient","0.5x Soil Moisture"),
       pch = c(16,16),
       col=c(4,6),cex=1.3, bty = "n"
)

#2x soilM
plot(df$timeNA, wet$CMIN, type="l", ylim= c(0,300), col =6, cex = 0.5, cex.axis = 1.3, cex.lab = 1.5, xlab = "Day of year", ylab = expression("C flux (mg " ~ m^{-2} ~ hr^{-1} ~ ")"))
points(df$timeNA, normal$CMIN, col = 4, cex = 0.5, type="l")
legend("topleft",c("Ambient","2x Soil Moisture"),
       pch = c(16,16),
       col=c(4,6),cex=1.3, bty = "n"
)

dev.off()

