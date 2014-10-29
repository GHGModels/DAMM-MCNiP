Yr = 3000;
T = 20;
soilM = 0.3; %0.229; %initial soil moisture in cm3 H20/ cm3 soil

%Code runs base model under normal and warmed temperatures
tic
%used parameter values from Allison et al. 2010 unless otherwise noted, N
%these parameter values are the default values for the base model
dt  = 0.1; %timestep interval units = hours
Nt =  8760*Yr;%number of timesteps model will run for
E=0; 
R = 0.008314; %gas constant used in Arrhenius equation (kJ mol-1 degree-1)
CN_s = 27.6;%C:N of SOM as given by Schimel & Weintraub 2003
CN_m = 10;%C:N of microbes
CN_enz = 3;%C:N of enzymes
p = 0.5; %fraction of  C  initally allocated to enz production
q = 0.5; %fraction of N initally  allocated to enz production
a = 0.5; %fraction of enzyme pool acting on SOC pool(1-a = fraction of enz pool acting on SON pool)

r_death = 0.00015;%Microbial biomass turnover rate, hours-1
r_ECloss = 0.001;%enzyme pool turnover rate hours-1
MIC_to_SOC  = 0.5; %proportion of dead microbial biomass that re-enters SOC pool, (1-MICtoSOC = proporation of microbial biomass that re-enters DOC pool)
MIC_to_SON = 0.5; %proportion of dead microbial biomass that re-enters SON pool, (1-MICtoSOC = proporation of microbial biomass that re-enters DON pool)
DOC_input = 0.0005;%0.0005;%external C input into DOC pool(e.g. root exudates,root turnover) mg C cm-3 soil hour-1
DON_input = DOC_input/CN_s;%external N input into DON pool(e.g.root turnover)mg N cm-3 soil hour-1 
Litter_C = 0.0005; %external C input into SOC pool(e.g. leaf litter, FWD) mg C cm-3 soil hour-1
Litter_N = Litter_C/CN_s;%external N input into SON pool (e.g.leaf litter, FWD) mg N cm-3 soil hour-1, 
%%not sure if this is the correct stoichiometry for litter inputs

%Vmax determined by Arrheinus equation, Km by a linear relationship with temp
%uptake kinetic temperature relationship parameters
A_UPT_C  = 100000000; %Arrhenius constant for uptake C vmax unit(mg DOC cm-3 soil hours-1)
Ea_UPT_C= 48;%activation energy for arrhenius equation( kJ mol-1)
b_UPT_C = 0.1;%0.01; %intercept, mg cm-3 degree-1
m_UPT_C = 0.01;%0.1;%slope, mg cm-3

%Depolymerization kinetic temperature relationship parameters
Vmax_0      = 100000000;%;%Arrhenius constantmg SOM cm-3 soil hours-1
Ea_up       =  48; %activation energy for arrhenius equation, kJ mol-1 
Km_0        = 0.0015;%500; %intercept, mg cm-3
Km_slope    = 0.00005;%5; %slope,mg cm-3 degree-1

%Coefficients for linear relationship between temperature and CUE
%b_CUE = 0.63;%intercept, mg C mg-1 soil
%m_CUE=-0.016;%slope, degree-1

%DAMM constants
Km_O2 = 0.121; %cm3 O2/cm3 air
Dgas = 1.67; 
O2airfrac = 0.209; % L O2/ L air 
BD = 0.8; %bulk density in g/cm3
PD = 2.52; %particle density in g/cm3
porosity = 1 - BD/PD; 
frac = 0.000414;%26.2467; %0.000414;
Dliq = 3.17;

%Variables, 
%%this section of code creates empty matrices so that values at each timestep can be saved during model runs, I save all the variables,
%including intermidiates such as uptake C (UPT_C) or less important variables, such as Death_C(microbial turnover) so that I can understand all model
%dynamics, but if the code is running slow, you can  create zero matrices just for the variables you're intersted in. 

%important pools
MIC_C = zeros(Nt,1); %Microbial biomass C
MIC_N =  zeros(Nt,1);%microbial biomass N
EC = zeros(Nt,1); %Enzyme pool
SOC = zeros(Nt,1); %SOC
SON = zeros(Nt,1); %SON
DOC = zeros(Nt,1); %DOC
DON = zeros(Nt,1); %DON
CMIN = zeros(Nt,1); %C mineralized
NMIN = zeros(Nt,1); %N mineralized
O2 = zeros(Nt,1); %concentration of O2
sol_SOC = zeros(Nt,1); %SOC that can be solubilized
sol_SON = zeros(Nt,1); %SON that can be solubilized
avail_SOC = zeros(Nt,1); %available SOC
avail_SON = zeros(Nt,1); %available SON

%other variables
DECOM_C = zeros(Nt,1); % C depolymerized by enzymes
DECOM_N =zeros(Nt,1);% N depolymerized by enzymes
UPT_C = zeros(Nt,1); %C taken up by microbes
UPT_N =zeros(Nt,1); %N taken up by microbes
Enz_C = zeros(Nt,1); %Amount of C taken up by microbes that is available to use for enzyme production
Enz_N = zeros(Nt,1); %Amount of N taken up by microbes that is available to use for enzyme production
EPROD = zeros(Nt,1); %enzymes produced
ELOSS = zeros(Nt,1); %enzymes that turned over
Growth_C = zeros(Nt,1); %Amount of C taken up by microbes that is available to use for microbial  growth
Growth_N = zeros(Nt,1);%Amount of N taken up by microbes that is available to use for microbial  growth
Growth = zeros(Nt,1); %Amount of new biomass grown 
DEATH_C = zeros(Nt,1); %Amount of microbial biomass C that turned over
DEATH_N = zeros(Nt,1);%Amount of microbial biomass N that turned over
overflow_C = zeros(Nt,1);%microbial overflow C metabolism
 
%Equilibrium Initial conditions, 
%determined after spinning up model for 5000 years, these are the inital values for each of the pools at time=0. 
MIC_C(1,:) = 1.1957; %0.5;
MIC_N(1,:) = 0.1196; %0.05;
SOC(1,:) = 144.5986; %100;
SON(1,:) = 5.4413; %3.6232; 
DOC(1,:) = 0.00091631; %0.5; 
DON(1,:)= 0.00049421; %0.0333; 
EC(1,:)= 0.0325; %0.01;  

for i = 1:Nt
%this section will calculate O2 concentration at time i
O2(i) = Dgas * O2airfrac * ((porosity-soilM)^(4/3));

%this section will calculate available substrate for depolymerization and
%uptake
sol_SOC(i) = frac*SOC(i);
sol_SON(i) = frac*SON(i);
avail_SOC(i) = Dliq*(soilM^3)*sol_SOC(i); %maybe need to initialize this pool -try it first
avail_SON(i) = Dliq*(soilM^3)*sol_SON(i);
    
%this section of code will calculate vmax  Km and CUE at 20C and 25C. 
% Equations for kinetic temperature relationships
%uptake kinetics(base model assumes C and N kinetics are equal)
vmax_UPT_C = A_UPT_C .* exp(-Ea_UPT_C./(R.*(T+273))); %temp sensitive according to arrhenius
km_UPT_C = b_UPT_C + m_UPT_C * T; %linear function of temp
CUE = 0.31; %b_CUE + m_CUE * T; %carbon use efficiency, linear function of temp

vmax_upt_N = vmax_UPT_C;
km_upt_N = km_UPT_C; 

%depolymerization kinetics (base model assumes C and N kinetics are equal)
Vmax_C = Vmax_0 * exp(-Ea_up./(R.*(T + 273)));
Km_C = 0.0025;%(Km_slope * T) + Km_0; %0.0025;

Vmax_N = Vmax_C; 
Km_N =  Km_C;

%This section of code calculates the changes in pool sizes over model time
%using a series of differential equations
            UPT_C(i) = MIC_C(i) .* vmax_UPT_C *DOC(i) ./ (km_UPT_C + DOC(i)); %microbial C uptake,michaelis-menton dynamics
            CMIN(i) =  UPT_C(i) .* (1-CUE); %C mineralization
            
            UPT_N(i) = MIC_N(i) .* vmax_upt_N *DON(i) ./ (km_upt_N + DON(i)); %microbial N uptake michaelis-menton dynamics
            DEATH_C(i) = r_death .* MIC_C(i); %microbial C turnover, first order process
            DEATH_N(i) = r_death .* MIC_N(i); %microbial N turnover,first order process
            
      %Resource Allocation
            %enzyme production
            Enz_C(i) = p.*(CUE.*UPT_C(i));%amount of C available for enzyme production after C allocated to C mineralization
            Enz_N(i) = q.*UPT_N(i);%amount of N available for enzyme production
            if (Enz_C(i)/CN_enz) >=Enz_N(i) %enz production N limited case
                EPROD(i) = Enz_N(i);%number of enzymes produced = N cost for enzyme production(1 enzyme made for every 1 N used)
            else %C limited case
                EPROD(i) = Enz_C(i)/CN_enz;%enz production C limited
            end
            %Amount of C and N available for growth dependent on how much is left over after enz production 
            Growth_C(i) = (1-p).*(UPT_C(i).*CUE)+ Enz_C(i) - CN_enz*EPROD(i);%Available C for biomass growth
            Growth_N(i) = (1-q).*UPT_N(i)+ Enz_N(i) - EPROD(i);%available N for biomass growth
            %micrbial growth
            if Growth_C(i)/CN_m >= Growth_N(i) %microbes N limited case
                Growth(i) = Growth_N(i);%amount of microbial N biomass equals cost of N to produce microbes 
            else
                Growth(i) = Growth_C(i)/CN_m;%microbes C limited case
                
            end
            %after enz prod and growth, N or C left over gets mineralized
            overflow_C(i) = Growth_C(i) - CN_m* Growth(i);%extra C after microbes produced goes to overflow metabolism
            NMIN(i) = Growth_N(i) - Growth(i);% N mineralization
            %Microbial biomass pools
            MIC_C(i+1)  =  MIC_C(i) + dt*(CN_m*Growth(i) - DEATH_C(i));%Microbial biomass C pool, Growth multiplied by C:N of microbes because Growth = N cost to build microbe, but C cost is greater
            MIC_N(i+1) = MIC_N(i) + dt*(Growth(i) - DEATH_N(i));%Microbial biomass N pool
            
            %Enzyme Pool & Turnover
            ELOSS(i) = r_ECloss * EC(i);%enzyme turnover
            EC(i+1) = EC(i) +  dt * (EPROD(i) - ELOSS(i));%enzyme pool
            
            %Depolymerization inputs,derived from Allison et al 2010.
            DECOM_C(i) = Vmax_C .* a*EC(i) .* (avail_SOC(i) ./(Km_C + avail_SOC(i)))*O2(i)/(Km_O2 + O2(i)); %depolymerization of SOC by enzymes
            DECOM_N(i) = Vmax_N.*(1-a)*EC(i).*(avail_SON(i)./(Km_N + avail_SON(i)))*O2(i)/(Km_O2 + O2(i)); %depolymerization of SON by enzymes
            
            %SOM pools
            SOC(i+1) = SOC(i) + dt * (Litter_C + DEATH_C(i) * MIC_to_SOC - DECOM_C(i));
            SON(i+1) = SON(i) + dt * (Litter_N + DEATH_N(i) * MIC_to_SON - DECOM_N(i));
 
            %Dissolved C&N pools
            DOC(i+1) = DOC(i) + dt * (DOC_input + E + DECOM_C(i) + DEATH_C(i)*(1-MIC_to_SOC) + (CN_enz/(1+CN_enz)).*ELOSS(i) - UPT_C(i));
            DON(i+1) = DON(i) + dt * (DON_input+ E/CN_s + DECOM_N(i) + DEATH_N(i) * (1-MIC_to_SON) + (1/CN_enz).*ELOSS(i)- UPT_N(i));
            %turnover from enzymes deposited here. Because enzymes are not split up into separate pools within the model, I calculate the
            %amount of C and N created by the turnover of an enzyme by using the C:N of the enzymes

end

%save vars.mat

%will create figures for each pool over time 

% figure
% plot(MIC_C,'LineWidth',3)
% legend('seasonal model')
% title('Microbial Biomass')
% xlabel('timesteps')
% ylabel('mg C/cm^3 soil') 
% 
% figure
% plot(EC,'LineWidth',3)
% legend('seasonal model')
% title('Enzyme pool')
% xlabel('timesteps')
% ylabel('mg C/cm^3 soil')

figure
plot(SOC,'LineWidth',3)
legend('seasonal model')
title('SOC')
xlabel('timesteps')
ylabel('mg C/cm^3 soil')

% figure
% plot(SON,'LineWidth',3)
% legend('seasonal model')
% title('SON')
% xlabel('timesteps')
% ylabel('mg N/cm^3 soil')
% 
% figure
% plot(DOC,'LineWidth',3)
% legend('seasonal model')
% title('DOC')
% xlabel('timesteps')
% ylabel('mg C/cm^3 soil')

% figure
% plot(DON,'LineWidth',3)
% legend('seasonal model')
% title('DON')
% xlabel('timesteps')
% ylabel('mg N/cm^3 soil')
% 
% figure
% plot(CMIN,'LineWidth',3)
% legend('seasonal model')
% title('Soil respiration')
% xlabel('timesteps')
% ylabel('mg C/cm^3 soil/timestep')

% figure
% plot(NMIN,'LineWidth',3)
% legend('seasonal model')
% title('N mineralization')
% xlabel('timesteps')
% ylabel('mg N /cm^3 soil/timestep')
% % 
% figure
% plot(DECOM_C,'LineWidth',3)
% legend('seasonal model')
% title('C decomposition')
% xlabel('timesteps')
% ylabel('mg C /cm^3 soil/timestep')
% % 
% figure
% plot(UPT_C,'LineWidth',3)
% legend('seasonal model')
% title('C uptake')
% xlabel('timesteps')
% ylabel('mg C /cm^3 soil/timestep')
% 
% figure
% plot(Growth_C,'LineWidth',3)
% legend('seasonal model')
% title('C growth')
% xlabel('timesteps')
% ylabel('mg C /cm^3 soil/timestep')
% 
% figure
% plot(overflow_C,'LineWidth',3)
% legend('seasonal model')
% title('overflow C')
% xlabel('timesteps')
% ylabel('mg C /cm^3 soil/timestep')




toc

resp = sum(CMIN(Nt-8760:Nt,1));  %in mg/cm3/yr
nmin = sum(NMIN(Nt-8760:Nt,1));  %in mg/cm3/yr
decomp = sum(DECOM_C(Nt-8760:Nt,1));
avail = avail_SOC(Nt,1);
upt = sum(UPT_C(Nt-8760:Nt,1));
micc = MIC_C(Nt,1);
micn = MIC_N(Nt,1);
soc = SOC(Nt,1);
son = SON(Nt,1);
doc = DOC(Nt,1);
don = DON(Nt,1);
ec = EC(Nt,1);
resp
nmin



