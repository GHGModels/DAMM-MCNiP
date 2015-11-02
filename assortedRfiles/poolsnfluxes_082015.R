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

#save plots code 12 panel####
#pdf(("/Users/rzabramoff/Documents/records/dropbox backup/dropbox.6.1.2015/BU/Dissertation_Research/data 2014/Ch4/poolflux_031615.pdf"))

#shrinko####
par(mfrow=c(3,4), mai = c(0.6,0.6,0.3,0.1), mgp = c(2,0.7,0))
#par(mar = c(5,5,5,1))
# par(oma = c(0,0,0,0))
#SoilT
plot(df$timeNA, df$SoilT, type="p", main = "Soil temperature", 
     col =1, cex = 0.5, cex.axis = 1, cex.lab = 1.1, xaxt = "n",
     xlab = "", 
     ylab = expression("Soil temperature (" ~degree~ "C)"))
text(75,21.8, "a)", cex = 1.1, xpd = NA, font = 2)

#SoilM
plot(df$timeNA, df$SoilM, type="p", main = "Soil moisture", 
     col =1, cex = 0.5, cex.axis = 1, cex.lab = 1.1, xaxt = "n",
     xlab = "", 
     ylab = expression(theta~ "(" ~ cm^{3} ~ H[2] ~ "O" ~ cm^{-3} ~ "soil)"))
text(75,1.11, "b)", cex = 1.1, xpd = NA, font = 2)

#EX
plot(df$timeNA, normal$EX, type="p", main = "Root input C", 
     col =1, cex = 0.5, cex.axis = 1, cex.lab = 1.1, xaxt = "n",
     xlab = "",
     ylab = expression("C flux (mg " ~ m^{-2} ~ hr^{-1} ~ ")"))
text(75,7.1, "c)", cex = 1.1, xpd = NA, font = 2)

#NMIN
plot(df$timeNA, normal$NMIN, type="p", main = "N mineralization", 
     col =4, cex = 0.5, cex.axis = 1, cex.lab = 1.1, xaxt = "n",
     xlab = "",
     ylab = expression("N flux (mg " ~ m^{-2} ~ hr^{-1} ~ ")"))
text(75,3, "d)", cex = 1.1, xpd = NA, font = 2)

#DEP
plot(df$timeNA, normal$DEP, type="p", main = "Depolymerization", 
     col =4, cex = 0.5, cex.axis = 1, cex.lab = 1.1,xaxt = "n",
     xlab = "",
     ylab = expression("C flux (mg " ~ m^{-2} ~ hr^{-1} ~ ")"))
text(75,290, "e)", cex = 1.1, xpd = NA, font = 2)

#UPT
plot(df$timeNA, normal$UPT, type="p", main = "C Uptake", col =4,
     cex = 0.5, cex.axis = 1, cex.lab = 1.1, xaxt = "n",
     xlab = "",
     ylab = expression("C flux (mg " ~ m^{-2} ~ hr^{-1} ~ ")"))
text(75,242, "f)", cex = 1.1, xpd = NA, font = 2)

#CMIN
plot(df$timeNA, normal$CMIN, type="p", main = "C mineralization", 
     col =4, cex = 0.5, cex.axis = 1, cex.lab = 1.1, xaxt = "n",
     xlab = "",
     ylab = expression("C flux (mg " ~ m^{-2} ~ hr^{-1} ~ ")"))
text(75,168, "g)", cex = 1.1, xpd = NA, font = 2)

#over
plot(df$timeNA, normal$over, type="p", main = "Overflow C", 
     col =4, cex = 0.5, cex.axis = 1, cex.lab = 1.1, xaxt = "n",
     xlab = "",
     ylab = expression("C flux (mg " ~ m^{-2} ~ hr^{-1} ~ ")"))
text(75,35, "h)", cex = 1.1, xpd = NA, font = 2)

#SOC
plot(df$timeNA, normal$SOC, type="p", main = "SOC", col =6, 
     cex = 0.5, cex.axis = 1, cex.lab = 1.1, xlab = "Day of year", 
     ylab = expression("C pool (g " ~ m^{-2} ~ ")"))
text(65,538.71, "i)", cex = 1.1, xpd = NA, font = 2)

#DOC
plot(df$timeNA, normal$DOC, type="p", main = "DOC", col =6, 
     cex = 0.5, cex.axis = 1, cex.lab = 1.1, xlab = "Day of year", 
     ylab = expression("C pool (g " ~ m^{-2} ~ ")"))
text(75,0.335, "j)", cex = 1.1, xpd = NA, font = 2)

#MIC
plot(df$timeNA, normal$MIC, type="p", main = "Microbial Biomass C",
     col =6, cex = 0.5, cex.axis = 1, cex.lab = 1.1, 
     xlab = "Day of year", 
     ylab = expression("C pool (g " ~ m^{-2} ~ ")"))
text(75,95.08, "k)", cex = 1.1, xpd = NA, font = 2)

#ENZ
plot(df$timeNA, normal$ENZ, type="p", main = "Enzyme C", col =6,
     cex = 0.5, cex.axis = 1, cex.lab = 1.1, xlab = "Day of year", 
     ylab = expression("C pool (g " ~ m^{-2} ~ ")"))
text(75,2.106, "l)", cex = 1.1, xpd = NA, font = 2)

#endplots####
#dev.off()

#NMIN v overflow####
par(mfrow=c(1,1))
par(mar = c(5,5,5,1))
par(oma = c(0,0,0,0))
plot(normal$NMIN, normal$over, type="p", main = "Overflow C vs. N mineralization", 
     col =4, cex = 0.5, cex.axis = 1.3, cex.lab = 1.3,  
     xlab = expression("N flux (mg " ~ m^{-2} ~ hr^{-1} ~ ")"), ylab = expression("C flux (mg " ~ m^{-2} ~ hr^{-1} ~ ")"))

plot(normal$MIC, normal$SOC, type="p", main = "SOC vs. Microbial biomass C", 
     col =4, cex = 0.5, cex.axis = 1.3, cex.lab = 1.3, 
     xlab = expression("C pool (mg " ~ m^{-2} ~ ")"), ylab = expression("C pool (mg " ~ m^{-2} ~ ")"))

plot(df$timeNA, normal$CMIN/normal$NMIN, type="p", main = "Immobilization Index", 
     col =1, cex = 0.5, cex.axis = 1.3, cex.lab = 1.3, ylim=(c(0,10000)), 
     xlab = "Day of year", 
     ylab = expression("C flux / N flux (mg " ~ m^{-2} ~ hr^{-1} ~ ")"))


##soilM v Nmin and overflow (regression sucks cuz lags)
par(mfrow=c(3,1), mai = c(0.6,0.7,0.3,0.1), mgp = c(2,0.7,0))
plot(df$timeNA, df$SoilM, col =1, cex = 1, cex.axis = 1, cex.lab = 1.1, xaxt = "n", main = "Soil Moisture",
     xlab = "", 
     ylab = expression(theta~ "(" ~ cm^{3} ~ H[2] ~ "O" ~ cm^{-3} ~ "soil)"))
plot(df$timeNA, normal$over, col = 4, cex = 1, cex.axis = 1, cex.lab = 1.1, xaxt = "n", main = "Overflow Carbon",
     xlab = "",
     ylab = expression("C flux (mg " ~ m^{-2} ~ hr^{-1} ~ ")"))
plot(df$timeNA, normal$NMIN, col = 6, cex = 1, cex.axis = 1, cex.lab = 1.1, xaxt = "n", main = "N mineralization",
     xlab = "",
     ylab = expression("N flux (mg " ~ m^{-2} ~ hr^{-1} ~ ")"))

plot(df$SoilM, normal$NMIN)
plot(df$SoilM, normal$over)

summary(lm(normal$over~df$SoilM))
summary(lm(normal$NMIN~df$SoilM))

par(mfrow = c(1,1))
plot(normal$CMIN, normal$NMIN)

summary(lm(normal$NMIN~normal$CMIN))

#CUE calcs
#cue <- 1 - rco2/upt
cue <- 1- diff(normal$CMIN)/diff(normal$UPT)
tf <- df[1:4560,] #truncate df so can plot cue [df-1]
torm <- normal[1:4560,] #truncate norm so can plot cue [df-1]
plot(tf$SoilM, cue)
plot(torm$NMIN, cue)
plot(torm$over, cue)

bigdf <- cbind(df,normal)
bigdf$cue <- c(cue, 0.3117215)

library(ggplot2)
ggplot(bigdf, aes(x=CMIN, y=DEP, fill=SoilM)) +
  geom_point(shape=21, size=2.5) +
  scale_fill_gradient(low="blue", high="red")

ggplot(bigdf, aes(x=SoilT, y=NMIN, fill=SoilM)) +
  geom_point(shape=21, size=2.5) +
  scale_fill_gradient(low="blue", high="red")

p1 <- ggplot(bigdf, aes(x=SoilM, y=CMIN, fill=SoilT)) +
  geom_point(shape=21, size=2.5) +
  scale_fill_gradient(low="blue", high="red") +
  ylim(0,350)

p2 <- ggplot(bigdf, aes(x=SoilM, y=scale, fill=SoilT)) +
  geom_point(shape=21, size=2.5) +
  scale_fill_gradient(low="blue", high="red") +
  ylim(0,350)

library(grid)
library(gridExtra)
grid.arrange(p1, p2, ncol = 2)


