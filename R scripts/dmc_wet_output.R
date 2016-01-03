#ambient
cmin.b <- c(2.0173, 2.7425) #mg/cm3/yr #(bulk, rhizo) #200 year run time, last year double avg soilm
nmin.b <- c(0.0304, 0.0370) 

#wet
cmin.f <- c(4.8502, 5.4928) #c(2.0254, 2.7460) #mg/cm3/yr #200 year run time, last year double soilm
nmin.f <- c(0.0231, 0.0358) #c(0.0293, 0.0359)
 
#scale
cmin.sb <-cmin.b*(0.01*100*100*100/1000) #g/m2/yr
nmin.sb <-nmin.b*(0.01*100*100*100/1000)
cmin.sf <-cmin.f*(0.01*100*100*100/1000) #g/m2/yr
nmin.sf <-nmin.f*(0.01*100*100*100/1000)

cmin.vb <- cmin.sb[2]*frac+cmin.sb[1]*(1-frac)
nmin.vb <- nmin.sb[2]*frac+nmin.sb[1]*(1-frac)
cmin.vf <- cmin.sf[2]*frac+cmin.sf[1]*(1-frac)
nmin.vf <- nmin.sf[2]*frac+nmin.sf[1]*(1-frac)

cballs <- c(cmin.vb, cmin.vf)
nballs <- c(nmin.vb, nmin.vf)

#barplot
names <- c("Ambient", "2x soil moisture") #add expression
cee <- rep("C mineralization", 2)
enn <- rep("N mineralization", 2)
c.rr <- as.data.frame(c(cmin.vb, cmin.vf))
n.rr <- as.data.frame(c(nmin.vb, nmin.vf))
names(n.rr) <- names(c.rr) 
secs <- rbind(c.rr,n.rr)
mins <- as.data.frame(rbind(cbind(cee,names), cbind(enn,names)))
mins <- cbind(mins,secs)

pdf(("C:/Users/rose/Dropbox/current/root research/data 2014/Ch4/wet_031615.pdf"))
par(mfrow = c(1,2))
par(mar = c(3,5,4,1))
barplot(cballs,names.arg = names, ylab = expression("C flux (g " ~ m^{-2} ~ yr^{-1} ~ ")"), ylim = c(0,50), main = "C mineralization")
barplot(nballs,names.arg = names, ylab = expression("N flux (g " ~ m^{-2} ~ yr^{-1} ~ ")"), ylim = c(0, 1), main = "N mineralization")
dev.off()
