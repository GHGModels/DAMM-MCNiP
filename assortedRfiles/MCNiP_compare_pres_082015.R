setwd("/Users/rzabramoff/Documents/records/dropbox backup/dropbox.6.1.2015/BU/Dissertation_Research/data 2014/Ch4")

#upload####
df <-read.csv("/Users/rzabramoff/Documents/records/dropbox backup/dropbox.6.1.2015/BU/Dissertation_Research/data 2014/Ch4/mydataFlux.csv", sep=",",header=T, na.strings="NA")
mcnip <-read.csv("/Users/rzabramoff/Documents/records/dropbox backup/dropbox.6.1.2015/BU/Dissertation_Research/data 2014/Ch4/mcnip_032115.csv", sep=",",header=T, na.strings="NA")
dmc <-read.csv("/Users/rzabramoff/Documents/records/dropbox backup/dropbox.6.1.2015/BU/Dissertation_Research/data 2014/Ch4/dmc_032115.csv", sep=",",header=T, na.strings="NA") #200 year run time
non <-read.csv("/Users/rzabramoff/Documents/records/dropbox backup/dropbox.6.1.2015/BU/Dissertation_Research/data 2014/Ch4/non_032115.csv", sep=",",header=T, na.strings="NA")  #200 year run time

#conversions####
df$mcnip <- mcnip$CMIN*10000*10 #convert from mg/cm3/hr to mg/m2/hr to a 10cm detph
df$dmc <- dmc$CMIN*10000*10
df$non <- non$CMIN*10000*10

#plots####
df$timeNA <- ifelse(is.na(df$scale), NA, df$time)

#shrink dis####
par(mfrow=c(1,1))
par(mar = c(4,4,1,1))
par(oma = c(0,0,0,0))
#Trenched
plot(df$timeNA, df$scale, cex = 0.5, cex.axis = 1.3, cex.lab = 1.3, xlab = "Day of year", ylab = expression("C efflux (mg " ~ m^{-2} ~ hr^{-1} ~ ")"))
#points(df$timeNA, df$DAMM, col = 6, cex = 0.5)
legend(110,323,c("Trenched plots"),
       pch = c(16),
       col=c("black"),cex=1.3, bty = "n"
)

#DAMM
plot(df$timeNA, df$scale, cex = 0.5, cex.axis = 1.3, cex.lab = 1.3, xlab = "Day of year", ylab = expression("C efflux (mg " ~ m^{-2} ~ hr^{-1} ~ ")"))
points(df$timeNA, df$DAMM, col = 6, cex = 0.5)
legend(110,323,c("DAMM"),
       pch = c(16),
       col=c(6),cex=1.3, bty = "n"
)


#MCNiP
plot(df$timeNA, df$scale, cex = 0.5, cex.axis = 1.3, cex.lab = 1.3, xlab = "Day of year", ylab = expression("C efflux (mg " ~ m^{-2} ~ hr^{-1} ~ ")"))
points(df$timeNA, df$mcnip, col = 4, cex = 0.5)
legend(110,323,c("MCNiP"),
       pch = c(16),
       col=c(4),cex=1.3, bty = "n"
)


#DMC
plot(df$timeNA, df$scale, cex = 0.5, cex.axis = 1.3, cex.lab = 1.3, xlab = "Day of year", ylab = expression("C efflux (mg " ~ m^{-2} ~ hr^{-1} ~ ")"))
points(df$timeNA, df$dmc, col = 2, cex = 0.5)
legend(110,323,c("DAMM-MCNiP"),
       pch = c(16),
       col=c(2),cex=1.3, bty = "n"
)


#shrinked####

#calculate residuals vs soilM
fit <- lm(df$dmc ~ df$scale)
summary(fit)

fit2 <- lm(df$DAMM ~ df$scale)
summary(fit2)

# df$add <- df$dmc + 50
# 
# fit3 <- lm(df$add ~ df$scale) 
# summary(fit3)

reg <- lm(fit$residuals ~ df[complete.cases(df[,9]),]$SoilM)
reg2 <- lm(fit2$residuals ~ df[complete.cases(df[,c(9,11)]),]$SoilM)
trip <- df[complete.cases(df[,c(9,11)]),]

#plot residuals vs soilM for DAMM-MCNiP and DAMM####
par(mfrow = c(1,2))
plot(df[complete.cases(df[,9]),]$SoilM, fit$residuals, xlim = c(0.1,1.05), ylim = c(-155,155), ylab = "DAMM-MCNiP Residuals", xlab = expression(theta~ "(" ~ cm^{3} ~ H[2] ~ "O" ~ cm^{-3} ~ "soil)"))
abline(h = 0)
text(-0.05,165, "a)", cex = 1.5, xpd = NA, font = 2)

# plot(df[complete.cases(df[,9]),]$SoilM, fit3$residuals, xlim = c(0.1,1.05), ylim = c(-155,155), ylab = "DAMM-MCNiP Residuals", xlab = expression(theta~ "(" ~ cm^{3} ~ H[2] ~ "O" ~ cm^{-3} ~ "soil)"))
# abline(h = 0)
# text(-0.05,165, "a)", cex = 1.5, xpd = NA, font = 2)

plot(trip$SoilM, fit2$residuals, ylab = "DAMM Residuals", xlim = c(0.1,1.05), ylim = c(-155,155), xlab = expression(theta~ "(" ~ cm^{3} ~ H[2] ~ "O" ~ cm^{-3} ~ "soil)"))
abline(h = 0)
text(-0.05,165, "b)", cex = 1.5, xpd = NA, font = 2)


#more sensitive to temp or soilm?####
summary(lm(df$DAMM ~df$SoilT)) #0.80
summary(lm(df$DAMM ~df$SoilM)) #0.02
summary(lm(df$dmc ~df$SoilT))  #0.48
summary(lm(df$dmc ~df$SoilM)) #0.11
summary(lm(df$dmc ~df$SoilM*df$SoilT))
summary(lm(df$DAMM ~df$SoilM*df$SoilT))

cor(df$DAMM, df$SoilT, use = "complete")
cor(df$DAMM, df$SoilM, use = "complete")
cor(df$dmc, df$SoilT, use = "complete")
cor(df$dmc, df$SoilM, use = "complete")

newdf <- df[complete.cases(df[,9]),]
newdf$residuals <- fit$residuals
hi <- newdf[newdf$SoilT > 13.3,] #above 3rd quart
lo <- newdf[newdf$SoilT < 13.3,] #below 1st quart

par(mfrow = c(1,1))
plot(hi$SoilM, hi$residuals, ylim = c(-150,150), col = 2, ylab = "DAMM-MCNiP Residuals", xlab = expression(theta~ "(" ~ cm^{3} ~ H[2] ~ "O" ~ cm^{-3} ~ "soil)"))
points(lo$SoilM, lo$residuals, col =4)
abline(h = 0)

#newfig 
df$sub <- (df$dmc - df$DAMM)
plot(df$SoilM, df$sub)
