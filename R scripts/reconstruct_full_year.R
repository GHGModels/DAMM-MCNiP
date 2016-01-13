#reconstruct full year of temperature and soil moisture data using 
#linear regressions with Harvard Forest Fisher Met Station 15-minute data

#original df (to append to?)
#origdf <-read.csv("/Users/rzabramoff/Documents/records/dropbox backup/dropbox.6.1.2015/BU/Dissertation_Research/data 2014/Ch4/mydataFlux.csv", sep=",",header=T, na.strings="NA")
origdf <-read.csv("/Users/rzabramoff/Dropbox (Climate)/damm-mcnip/data & excel files/mydataFluxRight.csv", sep=",",header=T, na.strings="NA")

#import fisher met data
df <- read.csv("/Users/rzabramoff/Dropbox (Climate)/damm-mcnip/data & excel files/hf001-10-15min-m.csv")

#select out 2009, variables: temp, prec (and rh?)
df$timeobj <- strptime(as.character(df$datetime), format = "%Y-%m-%dT%H:%M")
df9 <- df[df$timeobj >= "2009-01-01 00:00:00" & df$timeobj <= "2009-12-31 23:45:00",]
df9.trim <- df9[17:35056,] 

#average into hourly for temp and prec
df9.trim$time <- format(df9.trim$timeobj, "%Y-%m-%d %H")
library(data.table)
dt <- data.table(df9.trim)
hf <- dt[,list(s10t = mean(s10t), prec = mean(prec),rh = mean(rh), dewp = mean(dewp), 
               bar = mean(bar), airt = mean(airt), slrr = mean(slrr), parr = mean(parr),
               netr = mean(netr)),by=time]
hf$doy <- strftime(hf$time, format = "%j")
hf$hour <- rep(0:23,365)
newdf <- hf[hf$doy >= 113 & hf$doy <= 303,]
newdf <- newdf[1:4561,]

#linear reg temperature with measured temperature
reg1 <- lm(origdf$SoilT ~ newdf$s10t)
summary(reg1)
reg1$coefficients[1]
reg1$coefficients[2]
plot(origdf$SoilT, newdf$s10t)

#correct fisher data back to orig site
hf$corrected.temp <- reg1$coefficients[1] + reg1$coefficients[2]*hf$s10t
#corrected is a full year of temp
#check by plotting
grizz = rgb(0.5,0.5,0.5)

#plot
#pdf(file = ("/Users/rzabramoff/Dropbox (Climate)/damm-mcnip/figures/measVcorrectedTemp.pdf"), height = 3.74016,
#    width = 4.52756)
par(mar = c(3,3,1,1))
par(oma = c(0,0,0,0))
plot(hf$doy, hf$corrected.temp, ylim = c(0,25), cex = 0.25, cex.axis = 0.75,
     pch = 16, ann = FALSE) 
points(newdf$doy, origdf$SoilT, col = grizz, cex = 0.25, pch = 16)
mtext(side = 1, text = "Day of Year", line = 2, cex = 0.75)
mtext(side = 2, text = expression("Soil Temperature ("~degree~"C)"), line = 2, cex = 0.75)
legend("topright", c("Measured", "Fisher Station (corrected)"),
       col = c(grizz,1), pch = c(16,16), cex = 0.75)
#dev.off()

# ####nope####
#linear reg soil moisture with prec 
 reg2 <- lm(origdf$SoilM ~ newdf$prec + newdf$dewp + newdf$bar + newdf$airt + newdf$slrr + newdf$parr + newdf$netr) 
reg2 <- lm(origdf$SoilM ~ newdf$prec) 
summary(reg2)

# #correct fisher data back to orig site # hf$corrected.soilm <-
# reg2$coefficients[1] + reg2$coefficients[2]*hf$prec +
# reg2$coefficients[3]*hf$dewp+ #   reg2$coefficients[4]*hf$bar +
# reg2$coefficients[5]*hf$airt + reg2$coefficients[6]*hf$slrr + #  
# reg2$coefficients[7]*hf$parr + reg2$coefficients[8]*hf$netr hf$corrected.soilm
# <- reg2$coefficients[1] + reg2$coefficients[2]*hf$prec #corrected is a full
# year of soilm...ish #check by plotting pdf(file = ("/Users/rzabramoff/Dropbox
# (Climate)/damm-mcnip/figures/measVcorrectedTemp.pdf"), height = 3.74016, width
# = 4.52756) par(mar = c(3,3,1,1)) par(oma = c(0,0,0,0)) plot(hf$doy,
# hf$corrected.soilm, ylim = c(0,1), cex = 0.25, cex.axis = 0.75, pch = 16, ann
# = FALSE) points(newdf$doy, origdf$SoilM, col = grizz, cex = 0.25, pch = 16) 
# mtext(side = 1, text = "Day of Year", line = 2, cex = 0.75) mtext(side = 2,
# text = expression("Soil Moisture (" ~ cm^{-3} ~ H[2]~"O" ~ cm^{-3} ~ "soil)"),
# line = 2, cex = 0.75) legend("topright", c("Measured", "Fisher Station
# (corrected)"), col = c(grizz,1), pch = c(16,16), cex = 0.75) dev.off()
# 
# #plot soilm against prec sel = 1 sel2 = 365 plot(hf[hf$doy > sel & hf$doy <
# sel2,]$doy, hf[hf$doy > sel & hf$doy < sel2,,]$corrected.soilm, ylim = c(0,1))
# points(newdf[newdf$doy > sel & newdf$doy < sel2,]$doy, origdf[origdf$DOY > sel
# & origdf$DOY < sel2,]$SoilM, col = 2)
# ####nope####

#soil moisture time
#import stream discharge data
sf <- read.csv("/Users/rzabramoff/Dropbox (Climate)/damm-mcnip/data & excel files/hf070-04-15min.csv")

#select out 2009
sf$timeobj <- strptime(as.character(sf$datetime), format = "%Y-%m-%dT%H:%M")
sf9 <- sf[sf$timeobj >= "2009-01-01 00:00:00" & sf$timeobj <= "2009-12-31 23:45:00",]
sf9.trim <- sf9[5:35044,] 

#average hourly for stages and discharge
sf9.trim$time <- format(sf9.trim$timeobj, "%Y-%m-%d %H")
st <- data.table(sf9.trim)
shf <- st[,list(nb.stg = mean(nb.stg), nl.stg = mean(nl.stg), bl.stg = mean(bl.stg), bu.stg = mean(bu.stg),
                bgs.stg = mean(bgs.stg), bvs.stg = mean(bvs.stg), nb.dis = mean(nb.dis), nl.dis = mean(nl.dis), 
                nt.dis = mean(nt.dis), bl.dis = mean(bl.dis), bu.dis = mean(bu.dis)),by=time]
shf$doy <- strftime(shf$time, format = "%j")
shf$hour <- rep(0:23,365)
newsf <- shf[shf$doy >= 113 & shf$doy <= 303,]
newsf <- newsf[1:4561,]

#preliminary plotting
par(mfrow = c(2,1))
plot(shf$doy, shf$bgs.stg)
plot(origdf$DOY, origdf$SoilM, xlim=c(0,365))
par(mfrow = c(1,1))

#linear reg stage with measured soil moisture
reg2 <- lm(origdf$SoilM ~ newsf$bgs.stg)
summary(reg2)
reg2$coefficients[1]
reg2$coefficients[2]

#correct stage data back to orig site
shf$corrected.soilm <- reg2$coefficients[1] + reg2$coefficients[2]*shf$bgs.stg
#corrected is a full year of soil moisture
#check by plotting

#plot
#pdf(file = ("/Users/rzabramoff/Dropbox (Climate)/damm-mcnip/figures/measVcorrectedSoilmRight.pdf"), height = 3.74016,
#    width = 4.52756)
par(mar = c(3,3,1,1))
par(oma = c(0,0,0,0))
plot(shf$doy, shf$corrected.soilm, ylim = c(0.2,0.3), cex = 0.25, cex.axis = 0.75,
     pch = 16, ann = FALSE) 
points(newsf$doy, origdf$SoilM, col = grizz, cex = 0.25, pch = 16)
mtext(side = 1, text = "Day of Year", line = 2, cex = 0.75)
mtext(side = 2, text = expression("Soil Moisture (" ~ cm^{-3} ~ H[2]~"O" ~ cm^{-3} ~ "soil)"), line = 2, cex = 0.75)
legend("topleft", c("Measured", "Black Gum Swamp (corrected)"),
       col = c(grizz,1), pch = c(16,16), cex = 0.75)
#dev.off()

#append temp and soilm to origdf and export as slightly modified filename
front <- hf[hf$doy < 113,]
back <- hf[hf$doy > 302,]
sfront <- shf[shf$doy < 113,]
sback <- shf[shf$doy > 302,]

gapfilled <- data.frame(matrix(NA, nrow=8759, ncol=3))
names(gapfilled) = c("doy","temp","soilm")
gapfilled$doy <- c(front$doy, newdf$doy, back$doy[2:length(back$doy)])
gapfilled$temp <- c(front$corrected.temp, origdf$SoilT, back$corrected.temp[2:length(back$corrected.temp)])
gapfilled$soilm <- c(sfront$corrected.soilm, origdf$SoilM, sback$corrected.soilm[2:length(sback$corrected.soilm)])

#alternate soilm dataset (use ratio between on and off growing season to modify intercept)
altdf <- read.csv("/Users/rzabramoff/Dropbox (Climate)/damm-mcnip/data & excel files/hf006-01-soil-respiration.csv")
plot(altdf$jd, altdf$vsm)
ends <- altdf[altdf$jd <113 | altdf$jd > 303,]
middle <- altdf[altdf$jd > 113 & altdf$jd < 303,]
mean(middle$vsm, na.rm=T)
mean(ends$vsm, na.rm=T)

write.csv(gapfilled, file = "/Users/rzabramoff/Dropbox (Climate)/damm-mcnip/data & excel files/mydataGapfilledRight.csv")
#go into matlab and propagate - use yellow sticky note totally professional workflow

#plot reconstructed as gray and real as black
#pdf(file = ("/Users/rzabramoff/Dropbox (Climate)/damm-mcnip/figures/reconstructedTemp.pdf"), height = 3.74016,
#    width = 4.52756)
par(mar = c(3,3,1,1))
par(oma = c(0,0,0,0))
plot(front$doy, front$corrected.temp, xlim = c(0,365), ylim = c(0,25), cex = 0.25, cex.axis = 0.75,
     pch = 16, ann = FALSE) 
points(newdf$doy, origdf$SoilT, cex = 0.25, pch = 16, col = grizz)
points(back$doy[2:length(back$doy)], back$corrected.temp[2:length(back$corrected.temp)],
      cex = 0.25, pch = 16)
mtext(side = 1, text = "Day of Year", line = 2, cex = 0.75)
mtext(side = 2, text = expression("Soil Temperature ("~degree~"C)"), line = 2, cex = 0.75)
legend("topright", c("Measured", "Fisher Station (corrected)"),
       col = c(grizz,1), pch = c(16,16), cex = 0.75)
#dev.off()

#pdf(file = ("/Users/rzabramoff/Dropbox (Climate)/damm-mcnip/figures/reconstructedSoilmRight.pdf"), height = 3.74016,
#    width = 4.52756)
par(mar = c(3,3,1,1))
par(oma = c(0,0,0,0))
plot(sfront$doy, sfront$corrected.soilm, xlim = c(0,365), ylim = c(0.2,0.3), cex = 0.25, cex.axis = 0.75,
     pch = 16, ann = FALSE) 
points(newsf$doy, origdf$SoilM, col = grizz, cex = 0.25, pch = 16)
points(sback$doy[2:length(sback$doy)], sback$corrected.soilm[2:length(sback$corrected.soilm)],
       cex = 0.25, pch = 16)
mtext(side = 1, text = "Day of Year", line = 2, cex = 0.75)
mtext(side = 2, text = expression("Soil Moisture (" ~ cm^{-3} ~ H[2]~"O" ~ cm^{-3} ~ "soil)"), line = 2, cex = 0.75)
legend("topleft", c("Measured", "Black Gum Swamp (corrected)"),
       col = c(grizz,1), pch = c(16,16), cex = 0.75)
#dev.off()
