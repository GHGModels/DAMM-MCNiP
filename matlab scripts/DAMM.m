%DAMM Stand-alone

%read in data
F = xlsread('hformat.xlsx');
soilT = rep(F(:,1),Yr);
soilM = rep(F(:,2),Yr);

%define parameters
 xAlphaSx = 1.0815e11; %5.38*10^10
 EaSx = 61.77; %72.26
 kMSx = 8.7e-4; %9.95*10^-7 

ac = DAMM_Cflux(xAlphaSx,EaSx,kMSx,soilT,soilM);
rc = ac./(10000*10);
one = cbind(c(0,400), c(0,400)); %fix syntax

%%%%
mydataFlux$scale <- mydataFlux$flux*10000*10
  
plot(mydataFlux$time, mydataFlux$scale, main = "DAMM", xlab = "Day of year", ylab = "C flux (mg/m2/hr)")
lines(mydataFlux$time, ac, col = 2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%define function %fix syntax
function [areaCflux]=DAMM_Cflux(xAlphaSx,EaSx,kMSx,soilT,soilM)
	R = 8.314472e-3; %kJ K-1 mol-1
	O2airfrac = 0.209; %L O2 L-1 air
	BD = 0.80; %bulk density of soil
	PD = 2.52; %particle density of soil
	porosity = 1-BD/PD; %total porosity
	Sxtot = 0.048; %C content (g/cm3)
	psx = 4.14e-4;
	Dliq = 3.17;
	Dgas = 1.67;
	kMO2 = 0.121;
	Soildepth = 10; %effective soil depth in cm

	Sx = Sxtot*psx*Dliq*(soilM)^3;
  O2 = Dgas*O2airfrac*((porosity - soilM)^(4/3));
	MMSx = Sx/(kMSx+Sx);
	MMO2 = O2/(kMO2+O2);
	VmaxSx = xAlphaSx*exp(-EaSx/(R*(soilT+273.15)));
	Resp = VmaxSx*MMSx*MMO2;
	areaCflux = 10000*Soildepth*Resp;
end
