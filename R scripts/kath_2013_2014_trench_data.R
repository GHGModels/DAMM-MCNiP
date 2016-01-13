kath <- read.csv("/Users/rzabramoff/Desktop/2013_2014HFtrenchedsummary_for Rose.csv", header=TRUE, sep=",", na.strings="NA", dec=".", strip.white=TRUE)

#kath$timeobj <- strptime(as.character(kath$DOYtime), format = "%j")
kath$timeobj <- strptime(as.character(kath$DateTime), format = "%m/%d/%y %H:%M")
kath$time <- format(kath$timeobj, "%Y-%m-%d %H")
sel <- kath[kath$Treatment == "T",]
#plot(sel$timeobj,sel$flux)
short <- na.omit(sel)
short$flux <- short$flux #mg m-2 hr-1 to mg cm-3 hr-1
thirteen <- short[short$Year == "2013",]
fourteen <- short[short$Year == "2014",]
#plot(thirteen$timeobj, thirteen$flux)
#plot(fourteen$timeobj, fourteen$flux)

#average thirteen and fourteen hourly
library(data.table)
dt <- data.table(thirteen)
d13 <- dt[,list(DOYtime = mean(DOYtime), flux = mean(flux), soilT_T = mean(soilT_T), 
               VSM_T = mean(VSM_T)), by=time]
d13$doy <- strftime(d13$time, format = "%j")
dt <- data.table(fourteen)
d14 <- dt[,list(DOYtime = mean(DOYtime), flux = mean(flux), soilT_T = mean(soilT_T), 
                VSM_T = mean(VSM_T)), by=time]
d14$doy <- strftime(d14$time, format = "%j")

write.csv(d13, file = "/Users/rzabramoff/Dropbox (Climate)/damm-mcnip/data & excel files/kath13.csv")
write.csv(d14, file = "/Users/rzabramoff/Dropbox (Climate)/damm-mcnip/data & excel files/kath14.csv")

# plot(kath$DOYtime , kath$flux) #gC m-2 day-1
# 
# plot(kath$soil.T_C, kath$flux)
# points(kath$soilT_T, kath$flux, col = 2)
# 
# plot(kath$VSM_C, kath$flux)
# points(kath$VSM_T, kath$flux, col = 2)

# library(rgl)
# plot3d(kath$soil.T_C, kath$VSM_C, kath$flux, lwd=1)
# 
# library(plot3D)
# scatter3D(z = kath$flux, x = kath$soil.T_C, y = kath$VSM_C, pch = 1, cex = 0.2, 
#           theta = 0, phi = 70, ticktype = "detailed",
#           xlab = "temp", ylab = "soilm", zlab = "flux", clab = "flux", colkey = list(length = 0.8, width = 0.4),            
#           main = "flux")


