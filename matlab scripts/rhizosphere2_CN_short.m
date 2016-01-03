function [] = rhizosphere2_CN(rootDOC,deep,deep2,filename,filename2,path);

       %Code runs base model under normal and added exudation

%used parameter values from Allison et al. 2010 unless otherwise noted, N
%these parameter values are the default values for the base model
dt  = 0.1; %timestep interval units = hours
Nt = 2000000;%788400;%number of timesteps model will run for
DOCexp=-0.05; %exponent multiplier to decline inputs exponentially with depth

T = 20; % Temperature degree C, 20C is default, 25C is warmed
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
DOC_input = 0.0005;%external C input into DOC pool(e.g. root exudates,root turnover) mg C cm-3 soil hour-1
DON_input = 0.0005/CN_s;%external N input into DON pool(e.g.root turnover)mg N cm-3 soil hour-1 
Litter_C =  0.0005; %external C input into SOC pool(e.g. leaf litter, FWD) mg C cm-3 soil hour-1
Litter_N = 0.0005/CN_s;%external N input into SON pool (e.g.leaf litter, FWD) mg N cm-3 soil hour-1, 
%%not sure if this is the correct stoichiometry for litter inputs

%Vmax determined by Arrheinus equation, Km by a linear relationship with temp

%uptake kinetic temperature relationship parameters
A_UPT_C  = 100000000; %Arrhenius constant for uptake C vmax unit(mg DOC cm-3 soil hours-1)
Ea_UPT_C= 48;%activation energy for arrhenius equation( kJ mol-1)
b_UPT_C = 0.1;%0.01; %intercept, mg cm-3 degree-1
m_UPT_C = 0.01;%0.1;%slope, mg cm-3

%Depolymerization kinetic temperature relationship parameters
Vmax_0      = 100000000;%;%Arrhenius constantmg SOM cm-3 soil hours-1
Ea_up       =  48;   %activation energy for arrhenius equation, kJ mol-1
Km_0        = 500; %intercept, mg cm-3
Km_slope    = 5; %slope,mg cm-3 degree-1

%Coefficients for linear relationship between temperature and CUE
b_CUE = 0.63;%intercept, mg C mg-1 soil
m_CUE=-0.016;%slope, degree-1


%Variables, 
%%this section of code creates empty matrices so that values at each timestep can be saved during model runs, I save all the variables,
%including intermidiates such as uptake C (UPT_C) or less important variables, such as Death_C(microbial turnover) so that I can understand all model
%dynamics, but if the code is running slow, you can  create zero matrices just for the variables you're intersted in. 

%important pools
BioC = zeros(Nt, length(deep)); %Microbial biomass C
MIC_N =  zeros(Nt, length(deep));%microbial biomass N
EC = zeros(Nt, length(deep)); %Enzyme pool
SOC = zeros(Nt, length(deep)); %SOC
SON = zeros(Nt, length(deep)); %SON
DOC = zeros(Nt, length(deep)); %DOC
DON = zeros(Nt, length(deep)); %DON
CO2all = zeros(Nt, length(deep)); %C mineralized
MIN = zeros(Nt, length(deep)); %N mineralized

%other variables
DECOMPc = zeros(Nt, length(deep)); % C depolymerized by enzymes
DECOM_N =zeros(Nt, length(deep));% N depolymerized by enzymes
UPT_C = zeros(Nt, length(deep)); %C taken up by microbes
UPT_N =zeros(Nt, length(deep)); %N taken up by microbes
EnzC = zeros(Nt, length(deep)); %Amount of C taken up by microbes that is available to use for enzyme production
Enz_N = zeros(Nt, length(deep)); %Amount of N taken up by microbes that is available to use for enzyme production
EPROD = zeros(Nt, length(deep)); %enzymes produced
ELOSS = zeros(Nt, length(deep)); %enzymes that turned over
Growth_C = zeros(Nt, length(deep)); %Amount of C taken up by microbes that is available to use for microbial  growth
Growth_N = zeros(Nt, length(deep));%Amount of N taken up by microbes that is available to use for microbial  growth
Growth = zeros(Nt, length(deep)); %Amount of new biomass grown 
DEATH_C = zeros(Nt, length(deep)); %Amount of microbial biomass C that turned over
DEATH_N = zeros(Nt, length(deep));%Amount of microbial biomass N that turned over
overflow_C = zeros(Nt, length(deep));%microbial overflow C metabolism
 
%Equilibrium Initial conditions, 
%determined after spinning up model for 2e6timesteps, these are the inital values for each of the pools at time=0. 
BioC(1,:) =  1.9356;% 0.5;% 1.9461;
MIC_N(1,:) =  0.1936;% 0.05;% 0.1946;
SOC(1,:) =  82.7878;% 100;% 80.9226;
SON(1,:) =   3.6320;% 3.6232;% 3.6321;
DOC(1,:) =  0.00073085;% 0.5;% 0.0007377; %0.0049;
DON(1,:)=    0.00042707;% 0.0333;% 0.00043935;  %0.0029;
EC(1,:)=    0.0381;% 0.01;% 0.0393;

%this section of code will calculate vmax  Km and CUE at 20C and 25C. 

% Equations for kinetic temperature relationships
%uptake kinetics(base model assumes C and N kinetics are equal)
vmax_UPT_C = A_UPT_C .* exp(-Ea_UPT_C./(R.*(T+273))); %temp sensitive according to arrhenius
km_UPT_C = b_UPT_C + m_UPT_C * T; %linear function of temp
CUE = b_CUE + m_CUE * T; %carbon use efficiency, linear function of temp

vmax_upt_N = vmax_UPT_C;
km_upt_N = km_UPT_C; 

%depolymerization kinetics (base model assumes C and N kinetics are equal)
Vmax_C = Vmax_0 * exp(-Ea_up./(R.*(T + 273)));
Km_C = (Km_slope * T) + Km_0;

Vmax_N = Vmax_C; 
Km_N =  Km_C;

%This section of code calculates the changes in pool sizes over model time
%using a series of differential equations
for j=1:length(deep)
for i = 1:Nt
    
            DOCmult = exp(DOCexp*deep(j));
            DONmult = DOCmult/CN_s; %increase CN from CN_s to 100
        
            UPT_C(i,j) = BioC(i,j) .* vmax_UPT_C * DOC(i,j) ./ (km_UPT_C + DOC(i,j)); %microbial C uptake,michaelis-menton dynamics
            CO2all(i,j) =  UPT_C(i,j) .* (1-CUE); %C mineralization
            
            UPT_N(i,j) = MIC_N(i,j) .* vmax_upt_N * DON(i,j) / (km_upt_N + DON(i,j)); %microbial N uptake michaelis-menton dynamics
            DEATH_C(i,j) = r_death .* BioC(i,j); %microbial C turnover, first order process
            DEATH_N(i,j) = r_death .* MIC_N(i,j); %microbial N turnover,first order process
            
      %Resource Allocation
            %enzyme production
            EnzC(i,j) = p.*(CUE.*UPT_C(i,j));%amount of C available for enzyme production after C allocated to C mineralization
            Enz_N(i,j) = q.*UPT_N(i,j);%amount of N available for enzyme production
            if (EnzC(i,j)/CN_enz) >=Enz_N(i,j) %enz production N limited case
                EPROD(i,j) = Enz_N(i,j);%number of enzymes produced = N cost for enzyme production(1 enzyme made for every 1 N used)
            else %C limited case
                EPROD(i,j) = EnzC(i,j)/CN_enz;%enz production C limited
            end
            %Amount of C and N available for growth dependent on how much is left over after enz production 
            Growth_C(i,j) = (1-p).*(UPT_C(i,j).*CUE)+ EnzC(i,j) - CN_enz*EPROD(i,j);%Available C for biomass growth
            Growth_N(i,j) = (1-q).*UPT_N(i,j)+ Enz_N(i,j) - EPROD(i,j);%available N for biomass growth
            %micrbial growth
            if Growth_C(i,j)/CN_m >= Growth_N(i,j) %microbes N limited case
                Growth(i,j) = Growth_N(i,j);%amount of microbial N biomass equals cost of N to produce microbes 
            else
                Growth(i,j) = Growth_C(i,j)/CN_m;%microbes C limited case
                
            end
            %after enz prod and growth, N or C left over gets mineralized
            overflow_C(i,j) = Growth_C(i,j) - CN_m* Growth(i,j);%extra C after microbes produced goes to overflow metabolism
            MIN(i,j) = Growth_N(i,j) - Growth(i,j);% N mineralization
            %Microbial biomass pools
            BioC(i+1,j)  =  BioC(i,j) + dt*(CN_m*Growth(i,j) - DEATH_C(i,j));%Microbial biomass C pool, Growth multiplied by C:N of microbes because Growth = N cost to build microbe, but C cost is greater
            MIC_N(i+1,j) = MIC_N(i,j) + dt*(Growth(i,j) - DEATH_N(i,j));%Microbial biomass N pool
            
            %Enzyme Pool & Turnover
            ELOSS(i,j) = r_ECloss * EC(i,j);%enzyme turnover
            EC(i+1,j) = EC(i,j) +  dt * (EPROD(i,j) - ELOSS(i,j));%enzyme pool
            
            %Depolymerization inputs,derived from Allison et al 2010.
            DECOMPc(i,j) = Vmax_C .* a*EC(i,j) .* (SOC(i,j) ./(Km_C + SOC(i,j)));%depolymerization of SOC by enzymes
            DECOM_N(i,j) = Vmax_N.*(1-a)*EC(i,j).*(SON(i,j)./(Km_N+SON(i,j)));%depolymerization of SON by enzymes
            
            %SOM pools
            SOC(i+1,j) = SOC(i,j) + dt * (Litter_C*DOCmult + DEATH_C(i,j) * MIC_to_SOC - DECOMPc(i,j));
            SON(i+1,j) = SON(i,j) + dt * (Litter_N*DOCmult + DEATH_N(i,j) *MIC_to_SON - DECOM_N(i,j));
            
            %Dissolved C&N pools
            DOC(i+1,j) = DOC(i,j) + dt * (DOC_input*DOCmult+ rootDOC*DOCmult + DECOMPc(i,j) + DEATH_C(i,j)*(1-MIC_to_SOC) + (CN_enz/(1+CN_enz)).*ELOSS(i,j) - UPT_C(i,j));
            DON(i+1,j) = DON(i,j) + dt * (DON_input*DOCmult+ rootDOC*DONmult + DECOM_N(i,j) + DEATH_N(i,j) * (1-MIC_to_SON) + (1/CN_enz).*ELOSS(i,j)- UPT_N(i,j));
            %turnover from enzymes deposited here. Because enzymes are not split up into separate pools within the model, I calculate the
            %amount of C and N created by the turnover of an enzyme by using the C:N of the enzymes

end
end
        
%year means
% save('C:\Users\rose\Desktop\test.mat');
DECOMPc_mean=zeros(228,deep2);
for i=1:228
    DECOMPc_mean(i,1:deep2)=mean(DECOMPc(i*8760-8759:i*8760,1:deep2),1);
end

%year means
min_mean=zeros(228,deep2);
for i=1:228
    min_mean(i,1:deep2)=mean(MIN(i*8760-8759:i*8760,1:deep2),1);
end

%year means
BioC_mean=zeros(228,deep2);
for i=1:228
    BioC_mean(i,1:deep2)=mean(BioC(i*8760-8759:i*8760,1:deep2),1);
end

%year means
SOC_mean=zeros(228,deep2);
for i=1:228
    SOC_mean(i,1:deep2)=mean(SOC(i*8760-8759:i*8760,1:deep2),1);
end


DECOMPc_mean_last_row = DECOMPc_mean(end,:);
min_mean_last_row = min_mean(end,:);
BioC_mean_last_row = BioC_mean(end,:);
SOC_mean_last_row = SOC_mean(end,:);




%save([path,filename]);
save([path,filename2],'BioC_mean_last_row','min_mean_last_row','DECOMPc_mean_last_row','SOC_mean_last_row');