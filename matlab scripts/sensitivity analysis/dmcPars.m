%define parameter inputs to dmc
function par=dmcPars()

par.Yr = 50;
par.F = xlsread('hformat.xlsx');
par.T = rep(par.F(:,1),par.Yr);
par.soilM = rep(par.F(:,2),par.Yr);
par.Nt =  length(par.T); % number of timesteps model will run for
par.BD = 0.8; %bulk density in g/cm3
par.PD = 2.52; %particle density in g/cm3
par.sat = 0.5; %1-BD/PD;

%set maximum soilM value, if greater than porosity = 0.68 (i.e. there is more water than space for it), the model will crash  
for i = 1:length(par.soilM)
 if par.soilM(i) >= par.sat;
                par.soilM(i) = par.sat;
            else 
               par.soilM(i) = par.soilM(i);
 end
end

%used parameter values from Allison et al. 2010 unless otherwise noted, N
%these parameter values are the default values for the base model
par.dt  = 0.1; %timestep interval units = hours
par.E= 0.00026; %exudate value <- can modify
par.CN_ex= 27.6;

par.R = 0.008314; %gas constant used in Arrhenius equation (kJ mol-1 degree-1)
par.CN_s = 27.6;%C:N of SOM as given by Schimel & Weintraub 2003
par.CN_m = 10;%C:N of microbes
par.CN_enz = 3;%C:N of enzymes
par.p = 0.5; %fraction of  C  initally allocated to enz production
par.q = 0.5; %fraction of N initally  allocated to enz production
par.a = 0.5; %fraction of enzyme pool acting on SOC pool(1-a = fraction of enz pool acting on SON pool)

par.r_death = 0.00015;%Microbial biomass turnover rate, hours-1
par.r_ECloss = 0.001;%enzyme pool turnover rate hours-1
par.MIC_to_SOC  = 0.5; %proportion of dead microbial biomass that re-enters SOC pool, (1-MICtoSOC = proporation of microbial biomass that re-enters DOC pool)
par.MIC_to_SON = 0.5; %proportion of dead microbial biomass that re-enters SON pool, (1-MICtoSOC = proporation of microbial biomass that re-enters DON pool)

par.DOC_input = 0.0005;%0.0005;%external C input into DOC pool(e.g. root exudates,root turnover) mg C cm-3 soil hour-1
par.DON_input = par.DOC_input/par.CN_s;%external N input into DON pool(e.g.root turnover)mg N cm-3 soil hour-1 
par.Litter_C = 0.0005; %external C input into SOC pool(e.g. leaf litter, FWD) mg C cm-3 soil hour-1
par.Litter_N = par.Litter_C/27.6;%external N input into SON pool (e.g.leaf litter, FWD) mg N cm-3 soil hour-1

%%not sure if this is the correct stoichiometry for litter inputs

%Vmax determined by Arrheinus equation, Km by a linear relationship with temp
%uptake kinetic temperature relationship parameters %%from Davidson et al 2012
par.A_UPT_C= 1.0815E11; %Arrhenius constant for uptake C vmax unit(mg DOC cm-3 soil hours-1)
par.Ea_UPT_C= 61.77; %activation energy for arrhenius equation(kJ mol-1) 
par.km_UPT_C = 0.3; %can make this temperature sensitive using code below if desired (mg cm-3)

%Depolymerization kinetic temperature relationship parameters %%from Davidson et al 2012
par.Vmax_0      = 1.0815E11; %Arrhenius constant mg SOM cm-3 soil hours-1
par.Ea_up       = 61.77; %activation energy for arrhenius equation, kJ mol-1 
par.Km_C = 0.0025;%can make this temperature sensitive using code below if desired
par.CUE = 0.31;

%DAMM constants
par.Km_O2 = 0.121; %cm3 O2/cm3 air
par.Dgas = 1.67; 
par.O2airfrac = 0.209; % L O2/ L air 
par.BD = 0.8; %bulk density in g/cm3
par.PD = 2.52; %particle density in g/cm3
par.porosity = 1 - par.BD/par.PD; 
par.frac = 0.000414;%26.2467; %0.000414;
par.Dliq = 3.17;

end
