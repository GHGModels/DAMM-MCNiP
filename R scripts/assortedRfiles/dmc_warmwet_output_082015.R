#upload####
#warm is at 5deg, sum and win warm is at 3deg
df <-read.csv("/Users/rzabramoff/Documents/records/dropbox backup/dropbox.6.1.2015/BU/Dissertation_Research/data 2014/Ch4/mydataFlux.csv", sep=",",header=T, na.strings="NA")

bulk <-read.csv("/Users/rzabramoff/Documents/records/dropbox backup/dropbox.6.1.2015/BU/Dissertation_Research/data 2014/Ch4/bulk_032115.csv", sep=",",header=T, na.strings="NA")
rhizo <-read.csv("/Users/rzabramoff/Documents/records/dropbox backup/dropbox.6.1.2015/BU/Dissertation_Research/data 2014/Ch4/rhizo_032115.csv", sep=",",header=T, na.strings="NA")

bulkwarm5 <-read.csv("/Users/rzabramoff/Documents/records/dropbox backup/dropbox.6.1.2015/BU/Dissertation_Research/data 2014/Ch4/bulkwarm5_032215.csv", sep=",",header=T, na.strings="NA")
rhizowarm5 <-read.csv("/Users/rzabramoff/Documents/records/dropbox backup/dropbox.6.1.2015/BU/Dissertation_Research/data 2014/Ch4/rhizowarm5_032215.csv", sep=",",header=T, na.strings="NA")
bulkcold5 <-read.csv("/Users/rzabramoff/Documents/records/dropbox backup/dropbox.6.1.2015/BU/Dissertation_Research/data 2014/Ch4/bulkcold5_032215.csv", sep=",",header=T, na.strings="NA")
rhizocold5 <-read.csv("/Users/rzabramoff/Documents/records/dropbox backup/dropbox.6.1.2015/BU/Dissertation_Research/data 2014/Ch4/rhizocold5_032215.csv", sep=",",header=T, na.strings="NA")

bulkwarm1.5x <-read.csv("/Users/rzabramoff/Documents/records/dropbox backup/dropbox.6.1.2015/BU/Dissertation_Research/data 2014/Ch4/bulkwarm1.5x_032215.csv", sep=",",header=T, na.strings="NA")
rhizowarm1.5x <-read.csv("/Users/rzabramoff/Documents/records/dropbox backup/dropbox.6.1.2015/BU/Dissertation_Research/data 2014/Ch4/rhizowarm1.5x_032215.csv", sep=",",header=T, na.strings="NA")
bulkcold0.5x <-read.csv("/Users/rzabramoff/Documents/records/dropbox backup/dropbox.6.1.2015/BU/Dissertation_Research/data 2014/Ch4/bulkwarm0.5x_032215.csv", sep=",",header=T, na.strings="NA")
rhizocold0.5x <-read.csv("/Users/rzabramoff/Documents/records/dropbox backup/dropbox.6.1.2015/BU/Dissertation_Research/data 2014/Ch4/rhizowarm0.5x_032215.csv", sep=",",header=T, na.strings="NA")


bulkwet <-read.csv("/Users/rzabramoff/Documents/records/dropbox backup/dropbox.6.1.2015/BU/Dissertation_Research/data 2014/Ch4/bulkwet_032215.csv", sep=",",header=T, na.strings="NA")
rhizowet <-read.csv("/Users/rzabramoff/Documents/records/dropbox backup/dropbox.6.1.2015/BU/Dissertation_Research/data 2014/Ch4/rhizowet_032215.csv", sep=",",header=T, na.strings="NA")
bulkdry <-read.csv("/Users/rzabramoff/Documents/records/dropbox backup/dropbox.6.1.2015/BU/Dissertation_Research/data 2014/Ch4/bulkdry_032215.csv", sep=",",header=T, na.strings="NA")
rhizodry <-read.csv("/Users/rzabramoff/Documents/records/dropbox backup/dropbox.6.1.2015/BU/Dissertation_Research/data 2014/Ch4/rhizodry_032215.csv", sep=",",header=T, na.strings="NA")

bulkdryhot <-read.csv("/Users/rzabramoff/Documents/records/dropbox backup/dropbox.6.1.2015/BU/Dissertation_Research/data 2014/Ch4/bulkdryhot_032215.csv", sep=",",header=T, na.strings="NA")
rhizodryhot <-read.csv("/Users/rzabramoff/Documents/records/dropbox backup/dropbox.6.1.2015/BU/Dissertation_Research/data 2014/Ch4/rhizodryhot_032215.csv", sep=",",header=T, na.strings="NA")

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

dryhot<- NULL
dryhot <- (rhizodryhot*frac+bulkdryhot*(1-frac))*(0.1*100*100*100)

#plots####
df$timeNA <- ifelse(is.na(df$scale), NA, df$time)

#save plots code 4 panel fluxes####
pdf(("/Users/rzabramoff/Documents/records/dropbox backup/dropbox.6.1.2015/BU/Dissertation_Research/data 2014/Ch4/warmwet_032215.pdf"))

#shrinkdis####
par(mfrow=c(2,2))
par(mar = c(5,5,2,1))
par(oma = c(0,0,0,0))
#warming C
plot(df$timeNA, warm5$CMIN, type="l", main = "C mineralization", ylim = c(0,350), col =6, cex = 0.5, cex.axis = 1.3, cex.lab = 1.3, xlab = "Day of year", ylab = expression("C flux (mg " ~ m^{-2} ~ hr^{-1} ~ ")"))
points(df$timeNA, normal$CMIN, col = 1, cex = 0.5, type="l")
points(df$timeNA, cold5$CMIN, col = 4, cex = 0.5, type="l")
legend("topleft",c("Ambient",expression("+5"~degree~ "C"), expression("-5"~degree~ "C")),
       pch = c(16,16,16),
       col=c(1,6,4),cex=1, bty = "n"
)
text(75,385, "a)", cex = 1.5, xpd = NA, font = 2)

#warming N
plot(df$timeNA, warm5$NMIN, type="l", main = "N mineralization", ylim = c(0,5), col =6, cex = 0.5, cex.axis = 1.3, cex.lab = 1.3, xlab = "Day of year", ylab = expression("N flux (mg " ~ m^{-2} ~ hr^{-1} ~ ")"))
points(df$timeNA, normal$NMIN, col = 1, cex = 0.5, type="l")
points(df$timeNA, cold5$NMIN, col = 4, cex = 0.5, type="l")
legend("topleft",c("Ambient",expression("+5"~degree~ "C"), expression("-5"~degree~ "C")),
       pch = c(16,16,16),
       col=c(1,6,4),cex=1, bty = "n"
)
text(75,5.5, "b)", cex = 1.5, xpd = NA, font = 2)

#soilM C
plot(df$timeNA, wet$CMIN, type="l", ylim = c(0,350), col =6, cex = 0.5, cex.axis = 1.3, cex.lab = 1.3, xlab = "Day of year", ylab = expression("C flux (mg " ~ m^{-2} ~ hr^{-1} ~ ")"))
points(df$timeNA, normal$CMIN, col = 1, cex = 0.5, type="l")
points(df$timeNA, dry$CMIN, col = 4, cex = 0.5, type="l")

legend("topleft",c("Ambient", expression("1.5x" ~theta), expression("0.5x" ~ theta)),
       pch = c(16,16,16),
       col=c(1,6,4),cex=1, bty = "n"
)
text(75,385, "c)", cex = 1.5, xpd = NA, font = 2)

#soilM N
plot(df$timeNA, wet$NMIN, type="l", ylim = c(0,5), col =6, cex = 0.5, cex.axis = 1.3, cex.lab = 1.3, xlab = "Day of year", ylab = expression("N flux (mg " ~ m^{-2} ~ hr^{-1} ~ ")"))
points(df$timeNA, normal$NMIN, col = 1, cex = 0.5, type="l")
points(df$timeNA, dry$NMIN, col = 4, cex = 0.5, type="l")

legend("topleft",c("Ambient", expression("1.5x" ~theta), expression("0.5x" ~ theta)),
       pch = c(16,16,16),
       col=c(1,6,4),cex=1, bty = "n"
)
text(75,5.5, "d)", cex = 1.5, xpd = NA, font = 2)
#shrinked####

dev.off()

#save plots code 4 panel diffs####
pdf(("/Users/rzabramoff/Documents/records/dropbox backup/dropbox.6.1.2015/BU/Dissertation_Research/data 2014/Ch4/warmwetdiff_032215.pdf"))

#shrinkdis####
par(mfrow=c(2,2))
par(mar = c(5,5,2,1))
par(oma = c(0,0,0,0))
#warming C
plot(df$timeNA, warm5$CMIN-normal$CMIN, type="l", main = "C mineralization", ylim = c(0,200), col =6, cex = 0.5, cex.axis = 1.3, cex.lab = 1.3, xlab = "Day of year", ylab = expression("C flux (mg " ~ m^{-2} ~ hr^{-1} ~ ")"))
points(df$timeNA, normal$CMIN-cold5$CMIN, col = 4, cex = 0.5, type="l")
legend("topleft",c(expression("+5"~degree~ "C - Ambient"), expression("Ambient - 5"~degree~ "C")),
       pch = c(16,16),
       col=c(6,4),cex=1, bty = "n"
)
text(75,220, "a)", cex = 1.5, xpd = NA, font = 2)

#warming N
plot(df$timeNA, warm5$NMIN-normal$NMIN, type="l", main = "N mineralization", ylim = c(-2,3), col =6, cex = 0.5, cex.axis = 1.3, cex.lab = 1.3, xlab = "Day of year", ylab = expression("N flux (mg " ~ m^{-2} ~ hr^{-1} ~ ")"))
points(df$timeNA, normal$NMIN-cold5$NMIN, col = 4, cex = 0.5, type="l")
legend("topleft",c(expression("+5"~degree~ "C - Ambient"), expression("Ambient - 5"~degree~ "C")),
       pch = c(16,16),
       col=c(6,4),cex=1, bty = "n"
)
text(75,3.4, "b)", cex = 1.5, xpd = NA, font = 2)

#soilM C
plot(df$timeNA, wet$CMIN-normal$CMIN, type="l", main = "C mineralization", ylim = c(0,200), col =6, cex = 0.5, cex.axis = 1.3, cex.lab = 1.3, xlab = "Day of year", ylab = expression("C flux (mg " ~ m^{-2} ~ hr^{-1} ~ ")"))
points(df$timeNA, normal$CMIN-dry$CMIN, col = 4, cex = 0.5, type="l")
legend("topleft",c(expression("1.5x" ~theta ~ "- Ambient"),expression("Ambient - 0.5x" ~theta)),
       pch = c(16,16),
       col=c(6,4),cex=1, bty = "n"
)
text(75,220, "c)", cex = 1.5, xpd = NA, font = 2)

#soilM N
plot(df$timeNA, wet$NMIN-normal$NMIN, type="l", main = "N mineralization", ylim = c(-2,3), col =6, cex = 0.5, cex.axis = 1.3, cex.lab = 1.3, xlab = "Day of year", ylab = expression("N flux (mg " ~ m^{-2} ~ hr^{-1} ~ ")"))
points(df$timeNA, normal$NMIN-dry$NMIN, col = 4, cex = 0.5, type="l")
legend("topleft",c(expression("1.5x" ~theta ~ "- Ambient"),expression("Ambient - 0.5x" ~theta)),
                   pch = c(16,16),
                   col=c(6,4),cex=1, bty = "n"
)
text(75,3.4, "d)", cex = 1.5, xpd = NA, font = 2)
#shrinked####

dev.off()

#plot cumulative mins####
pdf(("/Users/rzabramoff/Documents/records/dropbox backup/dropbox.6.1.2015/BU/Dissertation_Research/data 2014/Ch4/warmwetcumsum_032215.pdf"))
#shrinkdis####
par(mfrow=c(2,2))
par(mar = c(5,5,2,1))
par(oma = c(0,0,0,0))
#warming C
plot(df$timeNA, cumsum(warm5$CMIN*0.001), type="l", main = "C mineralization", ylim = c(0,400), col =6, cex = 0.5, cex.axis = 1.3, cex.lab = 1.3, xlab = "Day of year", ylab = expression("Cumulative C flux (g  " ~ m^{-2} ~ ")"))
points(df$timeNA, cumsum(normal$CMIN*0.001), col = 1, cex = 0.5, type="l")
points(df$timeNA, cumsum(cold5$CMIN*0.001), col = 4, cex = 0.5, type="l")
legend("topleft",c("Ambient",expression("+5"~degree~ "C"), expression("-5"~degree~ "C")),
       pch = c(16,16,16),
       col=c(1,6,4),cex=1, bty = "n"
)
text(68,440, "a)", cex = 1.5, xpd = NA, font = 2)

#warming N
plot(df$timeNA, cumsum(warm5$NMIN*0.001), type="l", main = "N mineralization", ylim = c(0,4), col =6, cex = 0.5, cex.axis = 1.3, cex.lab = 1.3, xlab = "Day of year", ylab = expression("Cumulative N flux (g  " ~ m^{-2} ~ ")"))
points(df$timeNA, cumsum(normal$NMIN*0.001), col = 1, cex = 0.5, type="l")
points(df$timeNA, cumsum(cold5$NMIN*0.001), col = 4, cex = 0.5, type="l")
legend("topleft",c("Ambient",expression("+5"~degree~ "C"), expression("-5"~degree~ "C")),
       pch = c(16,16,16),
       col=c(1,6,4),cex=1, bty = "n"
)
text(68,4.38, "b)", cex = 1.5, xpd = NA, font = 2)

#soilM C
plot(df$timeNA, cumsum(wet$CMIN*0.001), type="l", ylim = c(0,500), col =6, cex = 0.5, cex.axis = 1.3, cex.lab = 1.3, xlab = "Day of year", ylab = expression("Cumulative C flux (g  " ~ m^{-2} ~ ")"))
points(df$timeNA, cumsum(normal$CMIN*0.001), col = 1, cex = 0.5, type="l")
points(df$timeNA, cumsum(dry$CMIN*0.001), col = 4, cex = 0.5, type="l")

legend("topleft",c("Ambient", "1.5x Soil moisture","0.5x Soil moisture"),
       pch = c(16,16,16),
       col=c(1,6,4),cex=1, bty = "n"
)
text(68,545, "c)", cex = 1.5, xpd = NA, font = 2)

#soilM N
plot(df$timeNA, cumsum(wet$NMIN*0.001), type="l", ylim = c(0,4), col =6, cex = 0.5, cex.axis = 1.3, cex.lab = 1.3, xlab = "Day of year", ylab = expression("Cumulative N flux (g  " ~ m^{-2} ~ ")"))
points(df$timeNA, cumsum(normal$NMIN*0.001), col = 1, cex = 0.5, type="l")
points(df$timeNA, cumsum(dry$NMIN*0.001), col = 4, cex = 0.5, type="l")

legend("topleft",c("Ambient", "1.5x Soil moisture","0.5x Soil moisture"),
       pch = c(16,16,16),
       col=c(1,6,4),cex=1, bty = "n"
)
text(68,4.38, "d)", cex = 1.5, xpd = NA, font = 2)

#shrinked####
dev.off()


a<- sum(warm5$CMIN)*0.001 #mg/m2/hr to g/m2/yr
b<- sum(cold5$CMIN)*0.001
c<- sum(wet$CMIN)*0.001
d<- sum(dry$CMIN)*0.001
e <- sum(normal$CMIN)*0.001

a<- sum(warm5$NMIN)*0.001 #mg/m2/hr to g/m2/yr
b<- sum(cold5$NMIN)*0.001
c<- sum(wet$NMIN)*0.001
d<- sum(dry$NMIN)*0.001
e <- sum(normal$NMIN)*0.001

f<- sum(dryhot$CMIN)*0.001

#cV

sd(wet$NMIN)/mean(wet$NMIN)
sd(dry$NMIN)/mean(dry$NMIN)
sd(normal$NMIN)/mean(normal$NMIN)

#dryhot fig
pdf(("/Users/rzabramoff/Documents/records/dropbox backup/dropbox.6.1.2015/BU/Dissertation_Research/data 2014/Ch4/dryhot_032215.pdf"))

par(mfrow=c(1,1))
par(mar = c(5,5,2,1))
par(oma = c(0,0,0,0))
#dry hot (ambient line, hot line, dry hot line)
plot(df$timeNA, warm5$CMIN, type="l", ylim = c(0,350), col =6, cex = 0.5, cex.axis = 1.3, cex.lab = 1.3, xlab = "Day of year", ylab = expression("C efflux (mg " ~ m^{-2} ~ hr^{-1} ~ ")"))
points(df$timeNA, normal$CMIN, col = 1, cex = 0.5, type="l")
points(df$timeNA, dryhot$CMIN, col = 3, cex = 0.5, type="l")
legend("topleft",c("Ambient",expression("+5"~degree~ "C"), expression("+5"~degree~ "C and 0.90x" ~theta)),
       pch = c(16,16,16),
       col=c(1,6,3),cex=1, bty = "n"
)
dev.off()