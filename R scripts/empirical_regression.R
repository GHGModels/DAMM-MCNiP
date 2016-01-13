#develop linear empirical relationship

df <-read.csv("/Users/rzabramoff/Dropbox (Climate)/damm-mcnip/data & excel files/mydataFluxRight.csv", sep=",",header=T, na.strings="NA")

#soil T normal
#soil M has long tail

reg <- lm(df$scale ~ df$SoilT + df$SoilM)
regsum <- summary(reg)
rmse <- summary(reg)$sigma #RMSE
predict <- regsum$coefficients[1] + regsum$coefficients[2]*df$SoilT + regsum$coefficients[3]*df$SoilM

plot(df$DOY, df$scale)
points(df$DOY, predict, col = 2)
df$predictFit09 <- predict

kath13<-read.csv("/Users/rzabramoff/Dropbox (Climate)/damm-mcnip/data & excel files/kath13.csv", sep=",",header=T, na.strings="NA")
predict13 <- regsum$coefficients[1] + regsum$coefficients[2]*kath13$soilT_T + regsum$coefficients[3]*kath13$VSM_T
plot(kath13$DOYtime, kath13$flux)
points(kath13$DOYtime, predict13, col=2)

kath14<-read.csv("/Users/rzabramoff/Dropbox (Climate)/damm-mcnip/data & excel files/kath14.csv", sep=",",header=T, na.strings="NA")
predict14 <- regsum$coefficients[1] + regsum$coefficients[2]*kath14$soilT_T + regsum$coefficients[3]*kath14$VSM_T
plot(kath14$DOYtime, kath14$flux)
points(kath14$DOYtime, predict14, col=2)

kath13$predict13 <- predict13
kath14$predict14 <- predict14

#df1 <- df[,-c(13:14)]

#write.csv(df, file = "/Users/rzabramoff/Dropbox (Climate)/damm-mcnip/data & excel files/mydataFluxRight.csv")
#write.csv(df, file = "/Users/rzabramoff/Dropbox (Climate)/damm-mcnip/data & excel files/mydataFluxP.csv")
#write.csv(kath13, file = "/Users/rzabramoff/Dropbox (Climate)/damm-mcnip/data & excel files/kath13.csv")
#write.csv(kath14, file = "/Users/rzabramoff/Dropbox (Climate)/damm-mcnip/data & excel files/kath14.csv")
