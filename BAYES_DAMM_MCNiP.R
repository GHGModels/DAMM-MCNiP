library(rjags)

FitNorm <- " model {
sigma ~ dgamma(0.001,0.001)  ## prior residual precision
AUPTC ~ dnorm(1000000, 0.001)
EaUPTC ~ dnorm(50, 0.001)
kmUPTC ~ dnorm(0.1, 0.001)
R <- 0.008314
kmo2 <- 0.121
CN_m <-
r_death <-
DOC_input <- 0.0005
MICtoSON <- 0.5
CN_enz <- 3
r_ECloss <-

for(i in 1:n){  		## loop over data

mu[i] <- (1-CUE)*MIC_C[i]*AUPTC*exp(-EaUPTC/(R*(temp[i]+273)))*DOC[i]/(kmUPTC+ DOC[i])*(1.67*0.209*((0.6825-soilM[i])^4/3))/(kmo2+(1.67*0.209*((1-soilM[i])^4/3)))
CMIN[i] ~ dnorm(mu[i], sigma)
MIC_C[i] ~ MIC_C[i-1] + (CN_m*Growth[i] - DEATH_C[i])
DOC[i] ~ DOC[i-1] + (DOC_input + DECOM_C[i] + DEATH_C[i] + (1-MICtoSON) + (1/CN_enz)*ELOSS[i] - UPT_N[i])
Growth[i] ~ ##need a conditional statement 
DEATH_C[i] ~ r_death*MIC_C[i] 
ELOSS[i] ~ r_ECloss*EC[i]
UPT_N[i] ~ 
EC[i] ~ Ec[i-1] + EPROD[i] - ELOSS[i]
EPROD[i] ~   ##need conditional statement
}
}"
        
data = list(CMIN = mydataFlux$scale, temp = mydataFlux$SoilT, soilM = soilM, n=1)

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
                      n.iter = 10100, n.burnin = 100, n.thin=2)


## diagnostics of the MCMC
bmcmc = b1
plot(bmcmc)  		## mcmc history and density plot
#autocorr.plot(bmcmc)		## autocorrelation
#cumuplot(bmcmc)		## quantile plot
gelman.plot(bmcmc)		## GRB statistic
summary(bmcmc)		## summary table
sumb = summary(bmcmc)