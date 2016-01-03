library(rjags)
mydataFlux<-read.csv("C:/Users/rose/Dropbox/current/root research/data 2014/hfdata.csv", sep=",",header=T, na.strings="NA")

FitNorm <- " model {
  sigma ~ dgamma(0.001,0.001)  ## prior residual precision
  AUPTC ~ dnorm(1000000, 0.00001)
  EaUPTC ~ dnorm(50, 0.001)
  kmUPTC ~ dnorm(0.1, 0.001)

for(i in 1:n){			## loop over reps
mu[i] <- 0.69*0.9144*AUPTC*exp(-EaUPTC/(0.008314*(temp[i]+273)))*1.908/(kmUPTC+ 1.908)*(1.67*0.209*((0.6825-soilM[i])^4/3))/(0.121+(1.67*0.209*((1-soilM[i])^4/3)))
CMIN[i] ~ dnorm(mu[i], sigma)
py[i] ~ dnorm(mu[i], sigma)
}
}"
        
data = list(CMIN = mydataFlux$scale, temp = mydataFlux$SoilT, soilM = soilM, n=4561)

init.cond1 <- list()
init.cond1[[1]] = list(sigma=0.1, AUPTC = 5E10, EaUPTC = 48, kmUPTC = 0.3)
init.cond1[[2]] = list(sigma=0.001, AUPTC = 1E11, EaUPTC = 67, kmUPTC = 0.3)
init.cond1[[3]] = list(sigma=1, AUPTC = 1E6, EaUPTC = 72, kmUPTC = 0.3)

init = init.cond1

## compile JAGS model
j.model2   <- jags.model (file = textConnection(FitNorm),
                          data = data,
                          inits = init,
                          n.chains = 3)
## burn-in
b1   <- coda.samples (model = j.model2,
                      variable.names = c("sigma","AUPTC", "EaUPTC","kmUPTC"),
                      n.iter = 101000, n.burnin = 100, n.thin=2)


## diagnostics of the MCMC
bmcmc = b1
plot(bmcmc)  		## mcmc history and density plot
#autocorr.plot(bmcmc)		## autocorrelation
#cumuplot(bmcmc)		## quantile plot
#gelman.plot(bmcmc)		## GRB statistic
summary(bmcmc)	## summary table

plot(density(as.matrix(bmcmc)[,4]))
plot(density(log10(as.matrix(bmcmc)[,4])))

