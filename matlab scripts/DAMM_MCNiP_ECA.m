Yr = 100;
F = xlsread('hformat.xlsx');
T = rep(F(:,1),Yr);
soilM = rep(F(:,2),Yr);
Fluxes = xlsread('../data & excel files/mydataFlux.xlsx');

% %expo = 0.1:0.01:0.6;
% %insert2 = mean(soilM) + mean(soilM).^(1./expo);
% %plot(yr1(314:500)); %where a nice soilm peak is...(7.75 days)
% insert1 = linspace(0.15,0.5,50); %length 50 lin up
% insert2 = logspace(log10(0.5),log10(0.15),142); %length 142 exp down %figure out!
% insert = [insert1';insert2'];
% %t=1:360;
% %insert = F(314:500,2); %mean(soilM) + 0.3.*sin(2.*pi/365.*t-pi/2);
% padding1 = rep(min(insert),round(length(F(:,1)).*0.5));
% padding2 = rep(min(insert),round((length(F(:,1)).*0.5)-(length(insert)+1)));
% soilM = vertcat(padding1,insert,padding2);
% soilM = rep(soilM,Yr);
% T = rep(mean(F(:,1)),length(soilM));
% %plot(soilM(Nt-4561:Nt))

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
E = 0.00026; %exudate value <- can modify
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

DOC_input = 0.0005;%0.0005;%external C input into DOC pool(e.g. root exudates,root turnover) mg C cm-3 soil hour-1
DON_input = DOC_input/CN_s;%external N input into DON pool(e.g.root turnover) mg N cm-3 soil hour-1 
Litter_C = 0.0005; %external C input into SOC pool(e.g. leaf litter, FWD) mg C cm-3 soil hour-1
Litter_N = Litter_C/CN_l;%external N input into SON pool (e.g.leaf litter, FWD) mg N cm-3 soil hour-1

%Vmax determined by Arrheinus equation, Km by a linear relationship with temp
%uptake kinetic temperature relationship parameters %%from Davidson et al 2012
A_UPT_C= 1.0815E11; %Arrhenius constant for uptake C vmax unit(mg DOC cm-3 soil hours-1)
Ea_UPT_C= 61.77; %activation energy for arrhenius equation( kJ mol-1) 
km_UPT_C = 0.3; %can make this temperature sensitive using code below if desired

%Depolymerization kinetic temperature relationship parameters %%from Davidson et al 2012
A_C      = 1.0815E11; %Arrhenius constant mg SOM cm-3 soil hours-1
Ea_C       = 61.77; %activation energy for arrhenius equation, kJ mol-1 
Km_C = 0.0025;%can make this temperature sensitive using code below if desired

%DAMM constants
Km_O2 = 0.121; %cm3 O2/cm3 air
Dgas = 1.67; 
O2airfrac = 0.209; % L O2/ L air 
porosity = 1 - BD/PD; 
frac = 0.000414;
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
%determined after spinning up model for 5000 years, these are the inital values for each of the pools at time=0. Commented values are spin up defaults.
MIC_C(1,:) = 1.1957; %0.5;
MIC_N(1,:) = 0.1196; %0.05;
SOC(1,:) = 144.5986; %100;
SON(1,:) = 5.4413; %3.6232; 
DOC(1,:) = 0.00091631; %0.5; 
DON(1,:)= 0.00049421; %0.0333; 
EC(1,:)= 0.0325; %0.01;  

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
CUE = 0.31; %CUE is set
%b_CUE + m_CUE * T(i); %linear function of temp option

%depolymerization kinetics (base model assumes C and N kinetics are equal)
Vmax_C = A_C * exp(-Ea_C./(R.*(T(i) + 273)));
%Km_C = 0.0025;%(Km_slope * T(i)) + Km_0; %linear function of temp option 

Vmax_N = Vmax_C; 
Km_N =  Km_C;

%This section of code calculates the changes in pool sizes over model time
%using a series of differential equations

            UPT_C(i) = MIC_C(i) .* vmax_UPT_C .* DOC(i) .* O2(i) / ...
                (km_UPT_C .* (1 + (MIC_C(i)./km_UPT_C) + (DOC(i)./km_UPT_C) + (O2(i)./Km_O2))); %microbial C uptake,ECA dynamics
            CMIN(i) =  UPT_C(i) .* (1-CUE); %C mineralization
            
            UPT_N(i) = MIC_N(i) .* vmax_upt_N .* DON(i) .* O2(i) /...
                (km_upt_N .* (1 + (MIC_N(i)./km_upt_N) + (DON(i)./km_upt_N) + (O2(i)./Km_O2))); %microbial N uptake ECA dynamics
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
            DECOM_C(i) = Vmax_C .* a * EC(i) .* sol_SOC(i) ./ (Km_C + sol_SOC(i) + EC(i)); %ECA depolymerization of SOC by enzymes
            DECOM_N(i) = Vmax_N .*(1-a)* EC(i) .* sol_SON(i) ./ (Km_N + sol_SON(i) + EC(i)); %ECA depolymerization of SON by enzymes
            
            %SOM pools
            SOC(i+1) = SOC(i) + dt * (Litter_C + DEATH_C(i) * MIC_to_SOC - DECOM_C(i));
            SON(i+1) = SON(i) + dt * (Litter_N + DEATH_N(i) * MIC_to_SON - DECOM_N(i));
 
            %Dissolved C&N pools
            DOC(i+1) = DOC(i) + dt * (DOC_input + E + DECOM_C(i) + DEATH_C(i)*(1-MIC_to_SOC) + (CN_enz/(1+CN_enz)).*ELOSS(i) - UPT_C(i));
            DON(i+1) = DON(i) + dt * (DON_input+ E/CN_ex + DECOM_N(i) + DEATH_N(i) * (1-MIC_to_SON) + (1/CN_enz).*ELOSS(i)- UPT_N(i));
            %turnover from enzymes deposited here. Because enzymes are not split up into separate pools within the model, I calculate the
            %amount of C and N created by the turnover of an enzyme by using the C:N of the enzymes

end

save vars.mat

%will create figures for each pool over time 
% % 
figure
plot(MIC_C,'LineWidth',3)
legend('seasonal model')
title('Microbial Biomass')
xlabel('timesteps')
ylabel('mg C/cm^3 soil') 
xlim([Nt-4561 Nt])
%hx = graph2d.constantline(length(T)-2281, 'LineStyle',':', 'Color',[.7 .7 .7]);
%changedependvar(hx,'x');
% % 
figure
plot(EC,'LineWidth',3)
legend('seasonal model')
title('Enzyme pool')
xlabel('timesteps')
ylabel('mg C/cm^3 soil')
xlim([Nt-4561 Nt])

figure
plot(SOC,'LineWidth',3)
legend('seasonal model')
title('SOC')
xlabel('timesteps')
ylabel('mg C/cm^3 soil')
xlim([Nt-4561 Nt])
% 
% figure
% plot(SON,'LineWidth',3)
% legend('seasonal model')
% title('SON')
% xlabel('timesteps')
% ylabel('mg N/cm^3 soil')
% xlim([Nt-4561 Nt])
% 
figure
plot(DOC,'LineWidth',3)
legend('seasonal model')
title('DOC')
xlabel('timesteps')
ylabel('mg C/cm^3 soil')
xlim([Nt-4561 Nt])
% 
% figure
% plot(DON,'LineWidth',3)
% legend('seasonal model')
% title('DON')
% xlabel('timesteps')
% ylabel('mg N/cm^3 soil')
% xlim([Nt-4561 Nt])
% % 
figure
plot(CMIN,'LineWidth',3)
hold on;
plot(Nt-4560:Nt,Fluxes(:,9)./(10000*10)); %sync up the x-axis
legend('seasonal model')
title('Soil respiration')
xlabel('timesteps')
ylabel('mg C/cm^3 soil/timestep')
xlim([Nt-4561 Nt])
% 
figure
plot(NMIN,'LineWidth',3)
legend('seasonal model')
title('N mineralization')
xlabel('timesteps')
ylabel('mg N /cm^3 soil/timestep')
xlim([Nt-4561 Nt])

figure
plot(DECOM_C,'LineWidth',3)
legend('seasonal model')
title('C decomposition')
xlabel('timesteps')
ylabel('mg C /cm^3 soil/timestep')
xlim([Nt-4561 Nt])
% 
% figure
% plot(DECOM_N,'LineWidth',3)
% legend('seasonal model')
% title('N decomposition')
% xlabel('timesteps')
% ylabel('mg N /cm^3 soil/timestep')
% xlim([Nt-4561 Nt])
% 
figure
plot(UPT_C,'LineWidth',3)
legend('seasonal model')
title('C uptake')
xlabel('timesteps')
ylabel('mg C /cm^3 soil/timestep')
xlim([Nt-4561 Nt])
% 
% figure
% plot(UPT_N,'LineWidth',3)
% legend('seasonal model')
% title('N uptake')
% xlabel('timesteps')
% ylabel('mg N /cm^3 soil/timestep')
% xlim([Nt-4561 Nt])
% 
% figure
% plot(DEATH_C,'LineWidth',3)
% legend('seasonal model')
% title('Microbial turnover C')
% xlabel('timesteps')
% ylabel('mg C /cm^3 soil/timestep')
% xlim([Nt-4561 Nt])
% 
% figure
% plot(DEATH_N,'LineWidth',3)
% legend('seasonal model')
% title('Microbial turnover N')
% xlabel('timesteps')
% ylabel('mg N /cm^3 soil/timestep')
% xlim([Nt-4561 Nt])

figure
plot(Growth_C,'LineWidth',3)
legend('seasonal model')
title('C growth')
xlabel('timesteps')
ylabel('mg C /cm^3 soil/timestep')
xlim([Nt-4561 Nt])
% 
% figure
% plot(overflow_C,'LineWidth',3)
% legend('seasonal model')
% title('overflow C')
% xlabel('timesteps')
% ylabel('mg C /cm^3 soil/timestep') 
% xlim([Nt-4561 Nt])
% 
figure
plot(soilM,'LineWidth',3)
legend('seasonal model')
title('SoilM')
xlabel('timesteps')
ylabel('mg C /cm^3 soil/timestep')
xlim([Nt-4561 Nt])

% figure
% plot(T,'LineWidth',3)
% legend('seasonal model')
% title('temperature')
% xlabel('timesteps')
% ylabel('mg C /cm^3 soil/timestep')
% xlim([Nt-4561 Nt])
% 
% figure
% plot(solSOC,'LineWidth',3)
% legend('seasonal model')
% title('sol_SOC')
% xlabel('timesteps')
% ylabel('mg C /cm^3 soil/timestep') 
% xlim([Nt-4561 Nt])

figure
plot(O2,'LineWidth',3)
legend('seasonal model')
title('O2')
xlabel('timesteps')
ylabel('mg C /cm^3 soil/timestep')
xlim([Nt-4561 Nt])
% 

%end clock
toc

%annual estimates at end of run
% resp = sum(CMIN(Nt-4561:Nt,1));%in mg/cm3/yr
% nmin = sum(NMIN(Nt-4561:Nt,1));
% decomp = sum(DECOM_C(Nt-4561:Nt,1));
% 
% resp
% nmin
%   headers = ['CMIN'];
%   heads = cellstr(headers);
%   csvwrite_with_headers('cmin_nov615_roots.csv',CMIN(Nt-4560:Nt),heads)

% write files
% headers = ['SOC'];
% heads = cellstr(headers);
% csvwrite_with_headers('DMC_soc.csv',SOC(Nt-4560:Nt),heads)
% 
% headers = ['CMIN'];
% heads = cellstr(headers);
% csvwrite_with_headers('DMC_cmin.csv',CMIN(Nt-4560:Nt),heads)




