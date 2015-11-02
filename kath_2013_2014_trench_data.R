kath <- read.csv("/Users/rzabramoff/Desktop/2013_2014HFtrenchedsummary_for Rose.csv", header=TRUE, sep=",", na.strings="NA", dec=".", strip.white=TRUE)
plot(kath$DOYtime , kath$flux) #gC m-2 day-1

plot(kath$soil.T_C, kath$flux)
points(kath$soilT_T, kath$flux, col = 2)

plot(kath$VSM_C, kath$flux)
points(kath$VSM_T, kath$flux, col = 2)

library(rgl)
plot3d(kath$soil.T_C, kath$VSM_C, kath$flux, lwd=1)

library(plot3D)
scatter3D(z = kath$flux, x = kath$soil.T_C, y = kath$VSM_C, pch = 1, cex = 0.2, 
          theta = 0, phi = 70, ticktype = "detailed",
          xlab = "temp", ylab = "soilm", zlab = "flux", clab = "flux", colkey = list(length = 0.8, width = 0.4),            
          main = "flux")


