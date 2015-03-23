#ambient
cmin.b <- c(2.0173, 2.7425) #mg/cm3/yr #(bulk, rhizo) #200 year run time, last year double avg soilm
nmin.b <- c(0.0304, 0.0370) 

#warm
cmin.f <- c(3.212, 3.9149) #c(2.0254, 2.7460) #mg/cm3/yr #200 year run time, only last year is warmed
nmin.f <- c(0.0336, 0.0417) #c(0.0293, 0.0359)

#wet
cmin.w <- c(4.8502, 5.4928) #c(2.0254, 2.7460) #mg/cm3/yr #200 year run time, last year double soilm
nmin.w <- c(0.0231, 0.0358) #c(0.0293, 0.0359)
 
#both
cmin.t <- c(8.5742, 8.6548) #c(2.0254, 2.7460) #mg/cm3/yr #200 year run time, last year double soilm
nmin.t <- c(0.0204, 0.0379) #c(0.0293, 0.0359)

#scale
cmin.sb <-cmin.b*(0.01*100*100*100/1000) #g/m2/yr
nmin.sb <-nmin.b*(0.01*100*100*100/1000)
cmin.sf <-cmin.f*(0.01*100*100*100/1000) #g/m2/yr
nmin.sf <-nmin.f*(0.01*100*100*100/1000)
cmin.sw <-cmin.w*(0.01*100*100*100/1000) #g/m2/yr
nmin.sw <-nmin.w*(0.01*100*100*100/1000)
cmin.st <-cmin.t*(0.01*100*100*100/1000) #g/m2/yr
nmin.st <-nmin.t*(0.01*100*100*100/1000)

cmin.vb <- cmin.sb[2]*frac+cmin.sb[1]*(1-frac)
nmin.vb <- nmin.sb[2]*frac+nmin.sb[1]*(1-frac)
cmin.vf <- cmin.sf[2]*frac+cmin.sf[1]*(1-frac)
nmin.vf <- nmin.sf[2]*frac+nmin.sf[1]*(1-frac)
cmin.vw <- cmin.sw[2]*frac+cmin.sw[1]*(1-frac)
nmin.vw <- nmin.sw[2]*frac+nmin.sw[1]*(1-frac)
cmin.vt <- cmin.st[2]*frac+cmin.st[1]*(1-frac)
nmin.vt <- nmin.st[2]*frac+nmin.st[1]*(1-frac)

cballs <- c(cmin.vb, cmin.vf, cmin.vw, cmin.vt)
nballs <- c(nmin.vb, nmin.vf, nmin.vw, nmin.vt)

#barplot
names <- c("Ambient", "Warm", "Wet", "Both") #add expression

pdf(("C:/Users/rose/Dropbox/current/root research/data 2014/Ch4/warmwet_031615.pdf"))
par(mfrow = c(1,2))
par(mar = c(3,5,4,1))
barplot(cballs,names.arg = names, ylab = expression("C flux (g " ~ m^{-2} ~ yr^{-1} ~ ")"), ylim = c(0,100), main = "C mineralization")
barplot(nballs,names.arg = names, ylab = expression("N flux (g " ~ m^{-2} ~ yr^{-1} ~ ")"), ylim = c(0, 1), main = "N mineralization")
dev.off()
