library(plyr)
library(zoo)

df <- read.csv("/Users/rzabramoff/Documents/records/dropbox backup/dropbox.6.1.2015/BU/Dissertation_Research/data 2014/Ch4/DAMM-MCNiP/2009_LPH_data_for_Sudeep_DAMM_input.csv", header=TRUE, sep=",", na.strings="NA", dec=".", strip.white=TRUE) 
df$soilM = as.numeric(as.character(df$soilM))
sud <- df[df$Treatment == "T",]
#mydataFlux <-read.table("C:/Users/rose/Dropbox/current/root research/data 2013/Marc_Harvard_Forest_soil_temperature.txt", header=TRUE, sep="\t", na.strings="NA", dec=".", strip.white =TRUE)
sud$TIME <- as.character(sud$TIME)
sud$TIME <- format(sud$TIME, format = "%H:%M", tz = "")
sud$HOUR= ifelse(substr(x=sud$TIME, start=2,stop=2)==":", substr(x=sud$TIME, start=1,stop=1), substr(x=sud$TIME, start=1,stop=2))

##something wrong with merge, losing some values!!
fix.dt <- data.frame(DOYTIME = round(seq(112,303, by = 1/48), digits = 4))
sud$DOYTIME <- round(sud$DOYTIME, digits = 4)
fix.table <- function(sud,fix.dt){
  return( merge(sud,fix.dt, by.x = "DOYTIME", by.y = "DOYTIME", all = TRUE))
}

df1 <- fix.table(sud,fix.dt)

###Gapfill df1
#interpolate soil temp and moisture
df <- df1[-c(1:48),]
df$SoilT.int <- na.approx(df$SoilT, rule=2)
df$SoilM.int <- na.approx(df$soilM, rule=2)
df$DOY.int <- floor(df$DOYTIME)
#interpolate hour
day <- read.csv("/Users/rzabramoff/Documents/records/dropbox backup/dropbox.6.1.2015/BU/Dissertation_Research/data 2014/Ch4/DAMM-MCNiP/oneday.csv", header=TRUE, sep=",", na.strings="NA", dec=".", strip.white=TRUE) 
day1 <- rep(day$col, times=191)
day1[9169] <- "0:00"
day1 <- day1[-c(1:48)]
df$HOURlong.int <- day1
df$HOUR.int <- ifelse(substr(x=df$HOURlong.int, start=2,stop=2)==":", substr(x=df$HOURlong.int, start=1,stop=1), substr(x=df$HOURlong.int, start=1,stop=2))


plot(df$DOYTIME, df$SoilT.int)
plot(df$DOYTIME, df$SoilM.int)
plot(df$DOYTIME, df$FLUX)

ht<- ddply(df, c("DOY.int", "HOUR.int"), summarise, N=length(SoilT.int), mean = mean(SoilT.int, na.rm=TRUE), sd=sd(SoilT.int, na.rm=TRUE), se=sd/sqrt(N))      
colnames(ht) = c("DOY", "hour", "N", "SoilT", "SoilT.sd", "SoilT.se")
hm<- ddply(df, c("DOY.int", "HOUR.int"), summarise, N=length(SoilM.int), mean = mean(SoilM.int, na.rm=TRUE), sd=sd(SoilM.int, na.rm=TRUE), se=sd/sqrt(N))      
colnames(hm) = c("DOY", "hour", "N", "SoilM", "SoilM.sd", "SoilM.se")
flux<- ddply(df, c("DOY.int", "HOUR.int"), summarise, N=length(FLUX), mean = mean(FLUX, na.rm=TRUE), sd=sd(FLUX, na.rm=TRUE), se=sd/sqrt(N))      
colnames(flux) = c("DOY", "hour", "N", "flux", "flux.sd", "flux.se")

plot(ht$SoilT)
plot(hm$SoilM)
plot(flux$flux)

measurement <- cbind(ht,hm)
hf <- measurement[,-c(3,5:9,11:12)]
hf$SoilM <- hf$SoilM/100       #soilM is in cm3 H20/cm3 soil ... soilT is in degC ... flux is in mgC/cm3/h assuming a soil depth of 10cm
hf$flux <- flux$flux*(1/(0.1*100*100*100))
hf$time <- hf$DOY + as.numeric(hf$hour)/48
hf$hour <- as.numeric(hf$hour)  
hformat <- hf[,c(3:4)]

dold <-read.csv("/Users/rzabramoff/Documents/records/dropbox backup/dropbox.6.1.2015/BU/Dissertation_Research/data 2014/Ch4/mydataFlux.csv")
dold$SoilMCrazy = dold$SoilM
dold$SoilM = hformat$SoilM
write.csv(dold,"/Users/rzabramoff/Dropbox (Climate)/damm-mcnip/data & excel files/mydataFluxRight.csv") 

