
# Dual Arrhenius and Michaelis-Menten (DAMM) kinetics model function coded for the statistical package R. The following code corresponds to the configuration of the
#  DAMM model used for modeling the field data described in the text. Arguments to the model function are the calibrated parameters, soil temperature (soilT, ?C) and
#  soil moisture (soilM, %). The fixed parameter values shown here are specific to the field data used in this study.

mydataFlux<-read.csv("/Users/rzabramoff/Dropbox (Climate)/damm-mcnip/data & excel files/mydataFluxRight.csv", sep=",",header=T, na.strings="NA")

soilT<-mydataFlux$SoilT
soilM<-mydataFlux$SoilM
xAlphaSx<-5.38*10^10 #1.0815e11#
EaSx<-72.26 #61.77#
kMSx<-9.95*10^-7 #8.7e-4# 
DAMM_Cflux <- function(xAlphaSx,EaSx,kMSx,soilT,soilM){
  R <- 8.314472e-3 #kJ K-1 mol-1
  O2airfrac <- 0.209 #L O2 L-1 air
  BD <- 0.80 #bulk density of soil
  PD <- 2.52 #particle density of soil
  porosity <- 1-BD/PD #total porosity
  Sxtot <- 0.048 #C content (g/cm3)
  psx <- 4.14e-4
  Dliq <- 3.17
  Dgas <- 1.67
  kMO2 <- 0.121
  Soildepth <- 10 #effective soil depth in cm
  
  Sx <- Sxtot*psx*Dliq*(soilM)^3
  O2 <- Dgas*O2airfrac*((porosity - soilM)^(4/3))
  MMSx <- Sx/(kMSx+Sx)
  MMO2 <- O2/(kMO2+O2)
  VmaxSx <- xAlphaSx*exp(-EaSx/(R*(soilT+273.15)))
  Resp <- VmaxSx*MMSx*MMO2
  areaCflux <- 10000*Soildepth*Resp
  return(areaCflux)
}

ac = DAMM_Cflux(xAlphaSx,EaSx,kMSx,soilT,soilM)
rc = ac/(10000*10)
one = cbind(c(0,400), c(0,400))
mydataFlux$scale <- mydataFlux$flux*10000*10

plot(mydataFlux$time, mydataFlux$scale, main = "DAMM", xlab = "Day of year", ylab = "C flux (mg/m2/hr)")
lines(mydataFlux$time, ac, col = 2)

plot(mydataFlux$scale, ac, main = "DAMM", xlab = "measured C flux (mg/m2/hr)", ylab = "DAMM C flux (mg/m2/hr)")
lines(one, col = 2)

mydataFlux$DAMM <- ac
newdf$goodwater <- ac

fit <- lm(mydataFlux$scale ~ ac)
summary(fit) ##R2 = 0.55

#write.csv(mydataFlux,"/Users/rzabramoff/Dropbox (Climate)/damm-mcnip/data & excel files/mydataFluxRight.csv")

###run DAMM on 2013 data
kath13<-read.csv("/Users/rzabramoff/Dropbox (Climate)/damm-mcnip/data & excel files/kath13.csv", sep=",",header=T, na.strings="NA")
soilT<-kath13$soilT_T
soilM<-kath13$VSM_T
ac = DAMM_Cflux(xAlphaSx,EaSx,kMSx,soilT,soilM)

plot(kath13$DOYtime, kath13$flux, main = "DAMM", xlab = "Day of year", ylab = "C flux (mg/m2/hr)")
lines(kath13$DOYtime, ac, col = 2)

plot(kath13$flux, ac, main = "DAMM", xlab = "measured C flux (mg/m2/hr)", ylab = "DAMM C flux (mg/m2/hr)")
lines(one, col = 2)
kath13$DAMM13 <- ac

fit <- lm(kath13$flux ~ ac)
summary(fit) ##R2 = 0.63

#write.csv(kath13,"/Users/rzabramoff/Dropbox (Climate)/damm-mcnip/data & excel files/kath13.csv")

###run DAMM on 2014 data
kath14<-read.csv("/Users/rzabramoff/Dropbox (Climate)/damm-mcnip/data & excel files/kath14.csv", sep=",",header=T, na.strings="NA")
soilT<-kath14$soilT_T
soilM<-kath14$VSM_T
ac = DAMM_Cflux(xAlphaSx,EaSx,kMSx,soilT,soilM)

plot(kath14$DOYtime, kath14$flux, main = "DAMM", xlab = "Day of year", ylab = "C flux (mg/m2/hr)")
lines(kath14$DOYtime, ac, col = 2)

plot(kath14$flux, ac, main = "DAMM", xlab = "measured C flux (mg/m2/hr)", ylab = "DAMM C flux (mg/m2/hr)")
lines(one, col = 2)
kath14$DAMM14 <- ac

fit <- lm(kath14$flux ~ ac)
summary(fit) ##R2 = 0.22

#write.csv(kath14,"/Users/rzabramoff/Dropbox (Climate)/damm-mcnip/data & excel files/kath14.csv")

#mydataFlux <- mydataFlux[,-c(1:3,12,15:16)]
mydataFlux$DAMM09 <- mydataFlux$DAMM
mydataFlux$DAMM13 <- newdf$kath13
mydataFlux$DAMM14 <- newdf$kath14
#write.csv(mydataFlux,"/Users/rzabramoff/Dropbox (Climate)/damm-mcnip/data & excel files/mydataFluxRight.csv")


