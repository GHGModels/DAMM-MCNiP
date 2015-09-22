dat <- read.csv("/Users/rzabramoff/Documents/records/dropbox backup/dropbox.6.1.2015/BU/Dissertation_Research/data 2014/Ch4/decompvbioc.csv", header=TRUE, sep=",", na.strings="NA", dec=".", strip.white=TRUE)

#convert
dat$exscale <- dat$ex*(0.1*100*100*100)/(1000) #conv from mg/cm3/hr to mg/m2/hr
dat$dec <- dat$DECOMPc*(0.1*100*100*100)/(1000) 
dat$bio <- dat$BioC*(0.1*100*100*100)/(1000) 
dat$soc <- dat$SOC*(0.1*100*100*100)/(1000) 

dat$DB <- dat$dec/dat$bio
dat$Root.DOCN <- as.factor (dat$Root.DOCN) 

#plots####
#pdf(file="C:/Users/rose/Dropbox/current/root research/data 2014/Ch4/decompvbioc.pdf",width=7,height=5)

#shrinkdis####
par(mfrow=c(1,2))
par(mar = c(5,5,2,1))
par(oma = c(0,0,0,0))
#SOC
plot(dat[dat$Root.DOCN == 1,]$exscale, dat[dat$Root.DOCN == 1,]$soc, type="l", main = "SOC", col =2, cex = 0.5, cex.axis = 1, cex.lab = 1, xlab = expression("Root input (g " ~ m^{-2} ~ hr^{-1} ~ ")"), ylab = expression("SOC pool (g " ~ m^{-2} ~ ")"))
points(dat[dat$Root.DOCN == 3,]$exscale, dat[dat$Root.DOCN == 3,]$soc, type="l", col = 3, cex = 0.5)
points(dat[dat$Root.DOCN == 7,]$exscale, dat[dat$Root.DOCN == 7,]$soc, type="l", col = 4, cex = 0.5)
points(dat[dat$Root.DOCN == 10,]$exscale, dat[dat$Root.DOCN == 10,]$soc, type="l", col = 5, cex = 0.5)
points(dat[dat$Root.DOCN == 27.6,]$exscale, dat[dat$Root.DOCN == 27.6,]$soc, type="l", col = 6, cex = 0.5)
points(dat[dat$Root.DOCN == 100,]$exscale, dat[dat$Root.DOCN == 100,]$soc, type="l", col = 8, cex = 0.5)
legend("topright",c("1","3","7","10","27.6","100"),
       pch = c(16,16,16,16,16,16),
       col=c(2,3,4,5,6,8),cex=1, title = "C:N"
)
text(-10,601, "a)", cex = 1.1, xpd = NA, font = 2)


#DB
plot(dat[dat$Root.DOCN == 1,]$exscale, dat[dat$Root.DOCN == 1,]$DB, type="l", ylim = c(2,3.5), main = "Microbial efficiency", col =2, cex = 0.5, cex.axis = 1, cex.lab = 1, xlab = expression("Root input (g " ~ m^{-2} ~ hr^{-1} ~ ")"), ylab = expression("Depolymerization/MB (g C mg M" ~ B^{-1} ~ hr^{-1} ~ ")"))
points(dat[dat$Root.DOCN == 3,]$exscale, dat[dat$Root.DOCN == 3,]$DB, type="l", col = 3, cex = 0.5)
points(dat[dat$Root.DOCN == 7,]$exscale, dat[dat$Root.DOCN == 7,]$DB, type="l", col = 4, cex = 0.5)
points(dat[dat$Root.DOCN == 10,]$exscale, dat[dat$Root.DOCN == 10,]$DB, type="l", col = 5, cex = 0.5)
points(dat[dat$Root.DOCN == 27.6,]$exscale, dat[dat$Root.DOCN == 27.6,]$DB, type="l", col = 6, cex = 0.5)
points(dat[dat$Root.DOCN == 100,]$exscale, dat[dat$Root.DOCN == 100,]$DB, type="l", col = 8, cex = 0.5)
legend("topright",c("1","3","7","10","27.6","100"),
       pch = c(16,16,16,16,16,16),
       col=c(2,3,4,5,6,8),cex=1, title = "C:N"
)
text(-10,3.65, "b)", cex = 1.1, xpd = NA, font = 2)

#shrinkd####

#dev.off()

pdf(file="C:/Users/rose/Dropbox/current/root research/data 2014/Ch4/decompvbiocnew.pdf",width=7,height=5)

#shrinkdis####
par(mfrow=c(2,2))
par(mar = c(5,5,2,1))
par(oma = c(0,0,0,0))
#SOC
plot(dat[dat$Root.DOCN == 1,]$exscale, dat[dat$Root.DOCN == 1,]$soc, type="l", main = "SOC", col =2, xlim = c(0,130), cex = 0.5, cex.main = 1, cex.axis = 1, cex.lab =1, lwd = 2, xlab = expression("Root input (g " ~ m^{-2} ~ hr^{-1} ~ ")"), ylab = expression("SOC pool (g " ~ m^{-2} ~ ")"))
points(dat[dat$Root.DOCN == 3,]$exscale, dat[dat$Root.DOCN == 3,]$soc, type="l", col = 3, cex = 0.5, lwd = 2)
points(dat[dat$Root.DOCN == 7,]$exscale, dat[dat$Root.DOCN == 7,]$soc, type="l", col = 4, cex = 0.5, lwd = 2)
points(dat[dat$Root.DOCN == 10,]$exscale, dat[dat$Root.DOCN == 10,]$soc, type="l", col = 5, cex = 0.5, lwd = 2)
points(dat[dat$Root.DOCN == 27.6,]$exscale, dat[dat$Root.DOCN == 27.6,]$soc, type="l", col = 6, cex = 0.5, lwd = 2)
points(dat[dat$Root.DOCN == 100,]$exscale, dat[dat$Root.DOCN == 100,]$soc, type="l", col = 8, cex = 0.5, lwd = 2)
legend("topright",c("1","3","7","10","27.6","100"),
       pch =
         c(16,16,16
               ,16,16,16),
       col=c(2,3,4,5,6,8),cex=0.8, title = "C:N"
)
text(-10,601, "a)", cex = 1.3, xpd = NA, font = 2)

#DB
plot(dat[dat$Root.DOCN == 1,]$exscale, dat[dat$Root.DOCN == 1,]$DB, type="l", ylim = c(2,3.5), lwd = 2, main = "Microbial efficiency", col =2, cex = 0.5, cex.main = 1, cex.axis = 1, cex.lab = 1, xlab = expression("Root input (g " ~ m^{-2} ~ hr^{-1} ~ ")"), ylab = expression("Depolymerization/MB (g C mg M" ~ B^{-1} ~ hr^{-1} ~ ")"))
points(dat[dat$Root.DOCN == 3,]$exscale, dat[dat$Root.DOCN == 3,]$DB, type="l", col = 3, cex = 0.5, lwd = 2)
points(dat[dat$Root.DOCN == 7,]$exscale, dat[dat$Root.DOCN == 7,]$DB, type="l", col = 4, cex = 0.5, lwd = 2)
points(dat[dat$Root.DOCN == 10,]$exscale, dat[dat$Root.DOCN == 10,]$DB, type="l", col = 5, cex = 0.5, lwd = 2)
points(dat[dat$Root.DOCN == 27.6,]$exscale, dat[dat$Root.DOCN == 27.6,]$DB, type="l", col = 6, cex = 0.5, lwd = 2)
points(dat[dat$Root.DOCN == 100,]$exscale, dat[dat$Root.DOCN == 100,]$DB, type="l", col = 8, cex = 0.5, lwd = 2)
text(-10,3.65, "b)", cex = 1.3, xpd = NA, font = 2)

#B
plot(dat[dat$Root.DOCN == 3,]$soc, dat[dat$Root.DOCN == 3,]$bio, type="l", ylim = c(70,130), lwd = 2, xlim = c(180,575), main = "Microbial biomass", col =3, cex = 0.5,cex.main = 1,  cex.axis = 1, cex.lab = 1, ylab = expression("Microbial biomass (g " ~ m^{-2} ~ ")"), xlab = expression("SOC pool (g " ~ m^{-2} ~ ")"))
points(dat[dat$Root.DOCN == 1,]$soc, dat[dat$Root.DOCN == 1,]$bio, type="l", col =2, cex = 0.5, lwd = 2)
points(dat[dat$Root.DOCN == 7,]$soc, dat[dat$Root.DOCN == 7,]$bio, type="l", col = 4, cex = 0.5, lwd = 2)
points(dat[dat$Root.DOCN == 10,]$soc, dat[dat$Root.DOCN == 10,]$bio, type="l", col = 5, cex = 0.5, lwd = 2)
points(dat[dat$Root.DOCN == 27.6,]$soc, dat[dat$Root.DOCN == 27.6,]$bio, type="l", col = 6, cex = 0.5, lwd = 2)
points(dat[dat$Root.DOCN == 100,]$soc, dat[dat$Root.DOCN == 100,]$bio, type="l", col = 8, cex = 0.5, lwd = 2)
text(150,136, "c)", cex = 1.3, xpd = NA, font = 2)
#shrinted####

dev.off()
# #D
# plot(dat[dat$Root.DOCN == 3,]$soc, dat[dat$Root.DOCN == 3,]$dec, type="l", ylim = c(250,275), main = "Depolymerization", col =3, cex = 0.5, cex.axis = 1, cex.lab = 1, ylab = expression("Depolymerization (g " ~ m^{-2} ~ hr^{-1} ~ ")"), xlab = expression("SOC pool (g " ~ m^{-2} ~ ")"))
# points(dat[dat$Root.DOCN == 1,]$soc, dat[dat$Root.DOCN == 1,]$dec, type="l", col =2, cex = 0.5)
# points(dat[dat$Root.DOCN == 7,]$soc, dat[dat$Root.DOCN == 7,]$dec, type="l", col = 4, cex = 0.5)
# points(dat[dat$Root.DOCN == 10,]$soc, dat[dat$Root.DOCN == 10,]$dec, type="l", col = 5, cex = 0.5)
# points(dat[dat$Root.DOCN == 27.6,]$soc, dat[dat$Root.DOCN == 27.6,]$dec, type="l", col = 6, cex = 0.5)
# points(dat[dat$Root.DOCN == 100,]$soc, dat[dat$Root.DOCN == 100,]$dec, type="l", col = 8, cex = 0.5)
# text(601,247.5, "b)", cex = 1.1, xpd = NA, font = 2)
