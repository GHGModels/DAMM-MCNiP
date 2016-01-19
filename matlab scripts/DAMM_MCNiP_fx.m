function DAMM_MCNiP_fx(Yr, scal, imp, n, A, Ea, frac, cue, amp, refin)

switch n 
    case '2900' %badwater
        F = xlsread('hformat.xlsx');
        T = rep(F(:,1),Yr);
        soilM = rep(F(:,2).*scal+imp,Yr);
    case '2009'
        F = xlsread('../data & excel files/mydataGapfilledRight.xlsx');
        T = rep(F(:,3),Yr);
        soilM = rep(F(:,4).*scal+imp,Yr);
        Fluxes = xlsread('../data & excel files/mydataFlux.xlsx');
    case '2013'
        F = xlsread('../data & excel files/kath13.xlsx');
        T = rep(F(:,7),Yr);
        soilM = rep(F(:,8).*scal+imp,Yr); %x2 for some reason
        Fluxes = F(:,6);
    case '2014'
        F = xlsread('../data & excel files/kath14.xlsx');
        T = rep(F(:,7),Yr);
        soilM = rep(F(:,8).*scal+imp,Yr); %x2 for some reason
        Fluxes = F(:,6);
end

load('initsMM.mat')

BD = 0.8; %bulk density in g/cm3
PD = 2.52; %particle density in g/cm3
sat = 0.5; %1-BD/PD;

%set maximum soilM value, if greater than porosity = 0.68 (i.e. there is more water than space for it), the model will crash  
for i = 1:length(soilM)
 if soilM(i) >= sat;
                soilM(i) = sat;
            else 
               soilM(i) = soilM(i);
 end
end

%start clock
tic
%used parameter values from Allison et al. 2010 unless otherwise noted, N
%these parameter values are the default values for the base model
dt  = 0.1; %timestep interval units = hours
Nt =  length(T); % number of timesteps model will run for
E = 0; %0.00026; %exudate value <- can modify
CN_ex= 27.6; %C:N of exudates

R = 0.008314; %gas constant used in Arrhenius equation (kJ mol-1 degree-1)
CN_s = 27.6;%C:N of SOM as given by Schimel & Weintraub 2003
CN_l = 50;%C:N of litter
CN_m = 10;%C:N of microbes
CN_enz = 3;%C:N of enzymes
p = 0.5; %fraction of  C  initally allocated to enz production
q = 0.5; %fraction of N initally  allocated to enz production
a = 0.5; %fraction of enzyme pool acting on SOC pool(1-a = fraction of enz pool acting on SON pool)

r_death = 0.00015;%Microbial biomass turnover rate, hours-1
r_ECloss = 0.001;%enzyme pool turnover rate hours-1
MIC_to_SOC  = 0.5; %proportion of dead microbial biomass that re-enters SOC pool, (1-MICtoSOC = proporation of microbial biomass that re-enters DOC pool)
MIC_to_SON = 0.5; %proportion of dead microbial biomass that re-enters SON pool, (1-MICtoSOC = proporation of microbial biomass that re-enters DON pool)

%input seasonal parameters
A1=amp/2;         %seasonal amplitude
A2=0;                %daily amplitude
w1=2.*pi/Nt;
w2=2.*pi;
ref=refin/2;
t=1:Nt;
DOC_input = ref+A1.*sin(w1.*t-pi/2)+A2.*sin(w2.*t);%0.0005;%external C input into DOC pool(e.g. root exudates,root turnover) mg C cm-3 soil hour-1
DON_input = DOC_input/CN_s;%external N input into DON pool(e.g.root turnover) mg N cm-3 soil hour-1 
Litter_C = DOC_input; %external C input into SOC pool(e.g. leaf litter, FWD) mg C cm-3 soil hour-1
Litter_N = Litter_C/CN_l;%external N input into SON pool (e.g.leaf litter, FWD) mg N cm-3 soil hour-1

%Vmax determined by Arrheinus equation, Km by a linear relationship with temp
%uptake kinetic temperature relationship parameters %%from Davidson et al 2012
A_UPT_C= A; %2.43E10; %1.0815E11; %Arrhenius constant for uptake C vmax unit(mg DOC cm-3 soil hours-1)
Ea_UPT_C= 61.77; %59.19; %61.77; %activation energy for arrhenius equation( kJ mol-1) 
km_UPT_C = 0.3; %can make this temperature sensitive using code below if desired

%Depolymerization kinetic temperature relationship parameters %%from Davidson et al 2012
A_C = A_UPT_C; %2.43E10; %1.0815E11; %Arrhenius constant mg SOM cm-3 soil hours-1
Ea_C = Ea; %Ea_UPT_C; %59.19; %61.77; %activation energy for arrhenius equation, kJ mol-1 
Km_C = 0.0025;%can make this temperature sensitive using code below if desired

%DAMM constants
Km_O2 = 0.121; %cm3 O2/cm3 air
Dgas = 1.67; 
O2airfrac = 0.209; % L O2/ L air 
porosity = 1 - BD/PD; 
frac = frac; %0.000414;
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
%determined after spinning up model for 1000 years, these are the inital values for each of the pools at time=0. Commented values are spin up defaults.
%%good water & full year
MIC_C(1,:) = initsMM(1);
MIC_N(1,:) = initsMM(2);
SOC(1,:) = initsMM(3);
SON(1,:) = initsMM(4);
DOC(1,:) = initsMM(5);
DON(1,:)= initsMM(6);
EC(1,:)= initsMM(7); 

for i = 1:Nt
%this section will calculate O2 concentration at time i
O2(i) = Dgas * O2airfrac * ((porosity-soilM(i))^(4/3));

%this section will calculate available substrate and enzymes at reaction site for depolymerization
sol_SOC(i) = Dliq*(soilM(i)^3)*frac*SOC(i);
sol_SON(i) = Dliq*(soilM(i)^3)*frac*SON(i);
    
%this section of code will calculate vmax  Km and CUE.
% Equations for kinetic temperature relationships
%uptake kinetics(base model assumes C and N kinetics are equal)
vmax_UPT_C = A_UPT_C .* exp(-Ea_UPT_C./(R.*(T(i)+273))); %temp sensitive according to arrhenius
vmax_upt_N = vmax_UPT_C;
%km_UPT_C = b_UPT_C + m_UPT_C * T(i); %linear function of temp option
km_upt_N = km_UPT_C;
CUE = cue; %0.31; %CUE is set
%b_CUE + m_CUE * T(i); %linear function of temp option

%depolymerization kinetics (base model assumes C and N kinetics are equal)
Vmax_C = A_C * exp(-Ea_C./(R.*(T(i) + 273)));
%Km_C = 0.0025;%(Km_slope * T(i)) + Km_0; %linear function of temp option 

Vmax_N = Vmax_C; 
Km_N =  Km_C;

%This section of code calculates the changes in pool sizes over model time
%using a series of differential equations

            UPT_C(i) = MIC_C(i) .* vmax_UPT_C *(DOC(i) ./ (km_UPT_C + DOC(i))) .* O2(i)./(Km_O2 + O2(i)); %microbial C uptake,michaelis-menton dynamics
            CMIN(i) =  UPT_C(i) .* (1-CUE); %C mineralization
            
            UPT_N(i) = MIC_N(i) .* vmax_upt_N *(DON(i) ./ (km_upt_N + DON(i))) .* O2(i)./(Km_O2 + O2(i)); %microbial N uptake michaelis-menton dynamics
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
            DECOM_C(i) = Vmax_C .* a*EC(i) .*sol_SOC(i) ./(Km_C + sol_SOC(i)); %depolymerization of SOC by enzymes
            DECOM_N(i) = Vmax_N.*(1-a)*EC(i).*sol_SON(i)./(Km_N + sol_SON(i)); %depolymerization of SON by enzymes
            
            %SOM pools
            SOC(i+1) = SOC(i) + dt * (Litter_C(i) + DEATH_C(i) * MIC_to_SOC - DECOM_C(i));
            SON(i+1) = SON(i) + dt * (Litter_N(i) + DEATH_N(i) * MIC_to_SON - DECOM_N(i));
 
            %Dissolved C&N pools
            DOC(i+1) = DOC(i) + dt * (DOC_input(i) + E + DECOM_C(i) + DEATH_C(i)*(1-MIC_to_SOC) + (CN_enz/(1+CN_enz)).*ELOSS(i) - UPT_C(i));
            DON(i+1) = DON(i) + dt * (DON_input(i) + E/CN_ex + DECOM_N(i) + DEATH_N(i) * (1-MIC_to_SON) + (1/CN_enz).*ELOSS(i)- UPT_N(i));
            %turnover from enzymes deposited here. Because enzymes are not split up into separate pools within the model, I calculate the
            %amount of C and N created by the turnover of an enzyme by using the C:N of the enzymes

end

save(['vars_DMC_',n,'.mat'])

end
