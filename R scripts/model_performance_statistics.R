#model performance statistics
#read in model output DAMM and scale
df <-read.csv("/Users/rzabramoff/Dropbox (Climate)/damm-mcnip/data & excel files/mydataFluxRight.csv")
#read in DAMM-MCNiP ECA
dmc <-read.csv("/Users/rzabramoff/Dropbox (Climate)/damm-mcnip/matlab scripts/DMCECA_var_cmin.csv")
df$dmc <- dmc$DMC*1000
df$data <- df$scale*0.01

#linear regression full dataset
fit1 <- lm(predict/100 ~ data, data = df)
sumfit1 <- summary(fit1) 
rmse1 <- summary(fit1)$sigma #RMSE

fit2 <- lm(DAMM/100 ~ data, data = df)
sumfit2 <- summary(fit2) 
rmse2 <- summary(fit2)$sigma #RMSE

fit3 <- lm(dmc ~ data, data = df)
sumfit3 <- summary(fit3) 
rmse3 <- summary(fit3)$sigma #RMSE

plot(df$data, df$predict/100, xlim=c(0,3), ylim=c(0,3))
plot(df$data, df$dmc, xlim=c(0,3), ylim=c(0,3))

#linear regression partial dataset
fit1 <- lm(predict/100 ~ data, data = df[1:(4561-1032),])
sumfit1 <- summary(fit1) 
rmse1 <- summary(fit1)$sigma #RMSE

fit2 <- lm(DAMM/100 ~ data, data = df[1:(4561-1032),])
sumfit2 <- summary(fit2) 
rmse2 <- summary(fit2)$sigma #RMSE

fit3 <- lm(dmc ~ data, data = df[1:(4561-1032),])
sumfit3 <- summary(fit3) 
rmse3 <- summary(fit3)$sigma #RMSE

#cross-correlation
ccf(df$data, df$dmc, na.action = na.exclude)
ccf(df$data, df$DAMM/100, na.action = na.exclude)
ccf(df$data, df$predict/100, na.action = na.exclude)

#variance #maybe use
var(df$data, na.rm = T)
var(df$dmc)
var(df$DAMM/100, na.rm = T)
var(df$predict/100)

#covariance
cov(df$data, df$dmc, use = "complete.obs")
cov(df$data, df$DAMM/100, use = "complete.obs")
cov(df$data, df$predict/100, use = "complete.obs")

#correlation
cor(df$data, df$dmc, use = "complete.obs")
cor(df$data, df$DAMM/100, use = "complete.obs")
cor(df$data, df$predict/100, use = "complete.obs")

#sum up season C efflux
sel <-df[complete.cases(df),]
sum(sel$data)
sum(sel$predict/100)
sum(sel$DAMM/100)
sum(sel$dmc)

#sum up spring/summer C efflux
sel <-df[complete.cases(df),]
sel <-sel[sel$DOY < 260,]
sum(sel$data)
sum(sel$predict/100)
sum(sel$DAMM/100)
sum(sel$dmc)