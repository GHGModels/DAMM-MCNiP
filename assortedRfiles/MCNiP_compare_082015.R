#setwd("/Users/rzabramoff/Documents/records/dropbox backup/dropbox.6.1.2015/BU/Dissertation_Research/data 2014/Ch4")

#upload####
df <-read.csv("/Users/rzabramoff/Documents/records/dropbox backup/dropbox.6.1.2015/BU/Dissertation_Research/data 2014/Ch4/mydataFlux.csv", sep=",",header=T, na.strings="NA")
#mcnip <-read.csv("/Users/rzabramoff/Documents/records/dropbox backup/dropbox.6.1.2015/BU/Dissertation_Research/data 2014/Ch4/mcnip_031615.csv", sep=",",header=T, na.strings="NA")
#dmc <-read.csv("/Users/rzabramoff/Documents/records/dropbox backup/dropbox.6.1.2015/BU/Dissertation_Research/data 2014/Ch4/dmc_032015.csv", sep=",",header=T, na.strings="NA") #200 year run time
#non <-read.csv("/Users/rzabramoff/Documents/records/dropbox backup/dropbox.6.1.2015/BU/Dissertation_Research/data 2014/Ch4/non_032015.csv", sep=",",header=T, na.strings="NA")  #200 year run time
#with new soilM
mcnip <-read.csv("/Users/rzabramoff/Documents/records/dropbox backup/dropbox.6.1.2015/BU/Dissertation_Research/data 2014/Ch4/mcnip_032115.csv", sep=",",header=T, na.strings="NA")
dmc <-read.csv("/Users/rzabramoff/Documents/records/dropbox backup/dropbox.6.1.2015/BU/Dissertation_Research/data 2014/Ch4/dmc_032115.csv", sep=",",header=T, na.strings="NA") #200 year run time
non <-read.csv("/Users/rzabramoff/Documents/records/dropbox backup/dropbox.6.1.2015/BU/Dissertation_Research/data 2014/Ch4/non_032115.csv", sep=",",header=T, na.strings="NA")  #200 year run time

kath <- read.csv("/Users/rzabramoff/Dropbox (Personal)/bu/dissertation_research/ch2/Kath_resp.csv", header=TRUE, sep=",", na.strings="NA", dec=".", strip.white=TRUE)
kath$Date <- as.Date(kath$Date)
#kath$root <- kath$Control - kath$CorrectedFlux
plot(kath$Date, kath$CorrectedFlux) #gC m-2 day-1
#points(kath$Date, kath$Control, col = 2)
#points(kath$Date, kath$root, col = 4)

#conversions####
df$mcnip <- mcnip$CMIN*10000*10 #convert from mg/cm3/hr to mg/m2/hr to a 10cm detph
df$dmc <- dmc$CMIN*10000*10
df$non <- non$CMIN*10000*10
kath$flux <- kath$CorrectedFlux*1000/24

#plots####
df$timeNA <- ifelse(is.na(df$scale), NA, df$time)

#save plots code####
#pdf(("/Users/rzabramoff/Documents/records/dropbox backup/dropbox.6.1.2015/BU/Dissertation_Research/data 2014/Ch4/plot_032615.pdf"))
#shrink dis####
par(mfrow=c(2,2))
par(mar = c(5,5,1,1))
par(oma = c(0,0,0,0))
#DAMM
plot(df$timeNA, df$scale, cex = 0.5, cex.axis = 1.3, cex.lab = 1.3, xlab = "Day of year", ylab = expression("C efflux (mg " ~ m^{-2} ~ hr^{-1} ~ ")"))
points(df$timeNA, df$DAMM, col = 6, cex = 0.5)
legend(110,323,c("Trenched plots","DAMM"),
       pch = c(16,16),
       col=c("black",6),cex=1.3, bty = "n"
)
text(75,330, "a)", cex = 1.7, xpd = NA, font = 2)

#MCNiP
plot(df$timeNA, df$scale, cex = 0.5, cex.axis = 1.3, cex.lab = 1.3, xlab = "Day of year", ylab = expression("C efflux (mg " ~ m^{-2} ~ hr^{-1} ~ ")"))
points(df$timeNA, df$mcnip, col = 4, cex = 0.5)
legend(110,323,c("MCNiP"),
       pch = c(16),
       col=c(4),cex=1.3, bty = "n"
)
text(75,330, "b)", cex = 1.7, xpd = NA, font = 2)

#DMC
plot(df$timeNA, df$scale, cex = 0.5, cex.axis = 1.3, cex.lab = 1.3, xlab = "Day of year", ylab = expression("C efflux (mg " ~ m^{-2} ~ hr^{-1} ~ ")"))
points(df$timeNA, df$dmc, col = 2, cex = 0.5)
legend(110,323,c("DAMM-MCNiP"),
       pch = c(16),
       col=c(2),cex=1.3, bty = "n"
)
text(75,330, "c)", cex = 1.7, xpd = NA, font = 2)

#NON
plot(df$timeNA, df$scale, cex = 0.5, cex.axis = 1.3, cex.lab = 1.3, xlab = "Day of year", ylab = expression("C efflux (mg " ~ m^{-2} ~ hr^{-1} ~ ")"))
points(df$timeNA, df$non, col = 3, cex = 0.5)
legend(110,323,c("DAMM-MCP"),
       pch = c(16),
       col=c(3),cex=1.3, bty = "n"
)
text(75,330, "d)", cex = 1.7, xpd = NA, font = 2)
#shrinked####
#dev.off()

#plots with different xlim####
xlims = c(120,150)
#shrinkdis####
par(mfrow=c(2,2))
par(mar = c(5,5,1,1))
par(oma = c(0,0,0,0))
#DAMM
plot(df$timeNA, df$scale, cex = 0.5, cex.axis = 1.3, cex.lab = 1.3, xlim = xlims, xlab = "Day of year", ylab = expression("C efflux (mg " ~ m^{-2} ~ hr^{-1} ~ ")"))
points(df$timeNA, df$DAMM, col = 6, cex = 0.5)
legend(110,323,c("Trenched plots","DAMM"),
       pch = c(16,16),
       col=c("black",6),cex=1.3, bty = "n"
)
text(75,330, "a)", cex = 1.7, xpd = NA, font = 2)

#MCNiP
plot(df$timeNA, df$scale, cex = 0.5, cex.axis = 1.3, cex.lab = 1.3, xlim = xlims, xlab = "Day of year", ylab = expression("C efflux (mg " ~ m^{-2} ~ hr^{-1} ~ ")"))
points(df$timeNA, df$mcnip, col = 4, cex = 0.5)
legend(110,323,c("MCNiP"),
       pch = c(16),
       col=c(4),cex=1.3, bty = "n"
)
text(75,330, "b)", cex = 1.7, xpd = NA, font = 2)

#DMC
plot(df$timeNA, df$scale, cex = 0.5, cex.axis = 1.3, cex.lab = 1.3, xlim = xlims, xlab = "Day of year", ylab = expression("C efflux (mg " ~ m^{-2} ~ hr^{-1} ~ ")"))
points(df$timeNA, df$dmc, col = 2, cex = 0.5)
legend(110,323,c("DAMM-MCNiP"),
       pch = c(16),
       col=c(2),cex=1.3, bty = "n"
)
text(75,330, "c)", cex = 1.7, xpd = NA, font = 2)

#NON
plot(df$timeNA, df$scale, cex = 0.5, cex.axis = 1.3, cex.lab = 1.3, xlim = xlims, xlab = "Day of year", ylab = expression("C efflux (mg " ~ m^{-2} ~ hr^{-1} ~ ")"))
points(df$timeNA, df$non, col = 3, cex = 0.5)
legend(110,323,c("DAMM-MCP"),
       pch = c(16),
       col=c(3),cex=1.3, bty = "n"
)
text(75,330, "d)", cex = 1.7, xpd = NA, font = 2)
#shrinkd####

#linregs####
drfit <- lm(DAMM ~ scale, data = df)
dr <- summary(drfit) #R2 = 0.52
a <- summary(drfit)$sigma #RMSE

mrfit <- lm(mcnip ~ scale, data = df)
mr <- summary(mrfit) #R2 = 0.50
b <- summary(mrfit)$sigma

dmrfit <- lm(dmc ~ scale, data = df)
dmr <- summary(dmrfit) #R2 = 0.17
c <- summary(dmrfit)$sigma

dnrfit <- lm(non ~ scale, data = df)
dnr <- summary(dnrfit) #R2 = 0.16
d <- summary(dnrfit)$sigma

oneone <- cbind(c(0,500), c(0,500))

#save linreg plots code####
#pdf(("/Users/rzabramoff/Documents/records/dropbox backup/dropbox.6.1.2015/BU/Dissertation_Research/data 2014/Ch4/linreg_031615.pdf"))

par(mfrow=c(2,2))
par(mar = c(5,5,3,1))
par(oma = c(0,0,0,0))
#DAMM
plot(df$scale, df$DAMM, main = "DAMM", xlim = c(0,400), ylim = c(0,400), cex = 0.5, cex.axis = 1.2, cex.lab = 1.2, xlab = expression("Measured C efflux (mg " ~ m^{-2} ~ hr^{-1} ~ ")"), ylab = expression("Predicted C efflux (mg " ~ m^{-2} ~ hr^{-1} ~ ")"))
abline(drfit)
lines(oneone, lty = 2)
text(-75,475, "a)", cex = 1.7, xpd = NA, font = 2)

#MCNiP
plot(df$scale, df$mcnip, main = "MCNiP", xlim = c(0,400), ylim = c(0,400), cex = 0.5, cex.axis = 1.2, cex.lab = 1.2, xlab = expression("Measured C efflux (mg " ~ m^{-2} ~ hr^{-1} ~ ")"), ylab = expression("Predicted C efflux (mg " ~ m^{-2} ~ hr^{-1} ~ ")"))
abline(mrfit)
lines(oneone, lty = 2)
text(-75,475, "b)", cex = 1.7, xpd = NA, font = 2)

#DMC
plot(df$scale, df$dmc, main = "DAMM-MCNiP", xlim = c(0,400), ylim = c(0,400), cex = 0.5, cex.axis = 1.2, cex.lab = 1.2, xlab = expression("Measured C efflux (mg " ~ m^{-2} ~ hr^{-1} ~ ")"), ylab = expression("Predicted C efflux (mg " ~ m^{-2} ~ hr^{-1} ~ ")"))
abline(dmrfit)
lines(oneone, lty = 2)
text(-75,475, "c)", cex = 1.7, xpd = NA, font = 2)

#NON
plot(df$scale, df$non, main = "DAMM-MCP", xlim = c(0,400), ylim = c(0,400), cex = 0.5, cex.axis = 1.2, cex.lab = 1.2, xlab = expression("Measured C efflux (mg " ~ m^{-2} ~ hr^{-1} ~ ")"), ylab = expression("Predicted C efflux (mg " ~ m^{-2} ~ hr^{-1} ~ ")"))
abline(dnrfit)
lines(oneone, lty = 2)
text(-75,475, "d)", cex = 1.7, xpd = NA, font = 2)

#dev.off()

plot(df$SoilT, df$scale)
plot(df$SoilT, df$DAMM)
plot(df$SoilT, df$dmc)
plot(df$SoilT, df$non)

plot(df$SoilM, df$scale)
plot(df$SoilM, df$DAMM)
plot(df$SoilM, df$dmc)
plot(df$SoilM, df$non)
