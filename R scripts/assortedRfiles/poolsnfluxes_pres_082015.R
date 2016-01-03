#upload####
df <-read.csv("/Users/rzabramoff/Documents/records/dropbox backup/dropbox.6.1.2015/BU/Dissertation_Research/data 2014/Ch4/mydataFlux.csv", sep=",",header=T, na.strings="NA")

bulk <-read.csv("/Users/rzabramoff/Documents/records/dropbox backup/dropbox.6.1.2015/BU/Dissertation_Research/data 2014/Ch4/bulk_032215.csv", sep=",",header=T, na.strings="NA")
rhizo <-read.csv("/Users/rzabramoff/Documents/records/dropbox backup/dropbox.6.1.2015/BU/Dissertation_Research/data 2014/Ch4/rhizo_032215.csv", sep=",",header=T, na.strings="NA")

#bulk rhizo scale####
#convert mg/cm3/hr to mg/m2/hr
frac = 0.2048506
normal <- NULL
normal <- (rhizo*frac+bulk*(1-frac))*(0.1*100*100*100)

normal$SOC <- normal$SOC*0.001
normal$DOC <- normal$DOC*0.001
normal$MIC <- normal$MIC*0.001
normal$ENZ <- normal$ENZ*0.001
#plots####
df$timeNA <- ifelse(is.na(df$scale), NA, df$time)

#shrinko####
#par(mfrow=c(3,4), mai = c(0.6,0.6,0.3,0.1), mgp = c(2,0.7,0))
par(mfrow = c(1,4))
par(mar = c(4,5.2,3,2))
par(oma = c(0,0,0,0))

#DEP
plot(df$timeNA, normal$DEP, type="p", main = "Depolymerization", 
     col =4, cex = 0.5,  cex.axis = 1.7, cex.lab = 1.7, cex.main = 2,
     xlab = "Day of Year",
     ylab = expression("C flux (mg " ~ m^{-2} ~ hr^{-1} ~ ")"))

#DOC
plot(df$timeNA, normal$DOC, type="p", main = "DOC", col =6, 
     cex = 0.5,  cex.axis = 1.7, cex.lab = 1.7, cex.main = 2, xlab = "Day of year", 
     ylab = expression("C pool (g " ~ m^{-2} ~ ")"))

#CMIN
plot(df$timeNA, normal$CMIN, type="p", main = "C mineralization", 
     col =4, cex = 0.5,  cex.axis = 1.7, cex.lab = 1.7, cex.main = 2, 
     xlab = "Day of Year",
     ylab = expression("C flux (mg " ~ m^{-2} ~ hr^{-1} ~ ")"))

#NMIN
plot(df$timeNA, normal$NMIN, type="p", main = "N mineralization", 
     col =4, cex = 0.5, cex.axis = 1.7, cex.lab = 1.7, cex.main = 2,
     xlab = "Day of Year",
     ylab = expression("N flux (mg " ~ m^{-2} ~ hr^{-1} ~ ")"))

#shrinko####
#par(mfrow=c(3,4), mai = c(0.6,0.6,0.3,0.1), mgp = c(2,0.7,0))
par(mfrow = (c(1,2)))
par(mar = c(4,4,3,2))
par(oma = c(0,0,0,0))
#SoilT
plot(df$timeNA, df$SoilT, type="p", main = "Soil temperature", 
     col =1, cex = 0.5, cex.axis = 1.3, cex.lab = 1.3, cex.main = 1.5,
     xlab = "Day of Year", 
     ylab = expression("Soil temperature (" ~degree~ "C)"))

#SoilM
plot(df$timeNA, df$SoilM, type="p", main = "Soil moisture", 
     col =1, cex = 0.5, cex.axis = 1.3, cex.lab = 1.3, cex.main = 1.5,
     xlab = "Day of Year", 
     ylab = expression(theta~ "(" ~ cm^{3} ~ H[2] ~ "O" ~ cm^{-3} ~ "soil)"))



#endplots####

