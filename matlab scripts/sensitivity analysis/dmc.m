%damm-mcnip model
function out = dmc(par)

%par = dmcPars();

%Variables, 
%%this section of code creates empty matrices so that values at each timestep can be saved during model runs, I save all the variables,
%including intermidiates such as uptake C (UPT_C) or less important variables, such as Death_C(microbial turnover) so that I can understand all model
%dynamics, but if the code is running slow, you can  create zero matrices just for the variables you're intersted in. 

%important pools
MIC_C = zeros(par.Nt,1); %Microbial biomass C
MIC_N =  zeros(par.Nt,1);%microbial biomass N
EC = zeros(par.Nt,1); %Enzyme pool
SOC = zeros(par.Nt,1); %SOC
SON = zeros(par.Nt,1); %SON
DOC = zeros(par.Nt,1); %DOC
DON = zeros(par.Nt,1); %DON
CMIN = zeros(par.Nt,1); %C mineralized
NMIN = zeros(par.Nt,1); %N mineralized
O2 = zeros(par.Nt,1); %concentration of O2
sol_SOC = zeros(par.Nt,1); %SOC that can be solubilized
sol_SON = zeros(par.Nt,1); %SON that can be solubilized

%other variables
DECOM_C = zeros(par.Nt,1); % C depolymerized by enzymes
DECOM_N =zeros(par.Nt,1);% N depolymerized by enzymes
UPT_C = zeros(par.Nt,1); %C taken up by microbes
UPT_N =zeros(par.Nt,1); %N taken up by microbes
Enz_C = zeros(par.Nt,1); %Amount of C taken up by microbes that is available to use for enzyme production
Enz_N = zeros(par.Nt,1); %Amount of N taken up by microbes that is available to use for enzyme production
EPROD = zeros(par.Nt,1); %enzymes produced
ELOSS = zeros(par.Nt,1); %enzymes that turned over
Growth_C = zeros(par.Nt,1); %Amount of C taken up by microbes that is available to use for microbial  growth
Growth_N = zeros(par.Nt,1);%Amount of N taken up by microbes that is available to use for microbial  growth
Growth = zeros(par.Nt,1); %Amount of new biomass grown 
DEATH_C = zeros(par.Nt,1); %Amount of microbial biomass C that turned over
DEATH_N = zeros(par.Nt,1);%Amount of microbial biomass N that turned over
overflow_C = zeros(par.Nt,1);%microbial overflow C metabolism

%Equilibrium Initial conditions, 
%determined after spinning up model for 5000 years, these are the inital values for each of the pools at time=0. Commented values are spin up defaults.
MIC_C(1,:) = 1.1957; %0.5;
MIC_N(1,:) = 0.1196; %0.05;
SOC(1,:) = 144.5986; %100;
SON(1,:) = 5.4413; %3.6232; 
DOC(1,:) = 0.00091631; %0.5; 
DON(1,:)= 0.00049421; %0.0333; 
EC(1,:)= 0.0325; %0.01;  

for i = 1:par.Nt
%this section will calculate O2 concentration at time i
O2(i) = par.Dgas * par.O2airfrac * ((par.porosity-par.soilM(i))^(4/3));

%this section will calculate available substrate and enzymes at reaction site for depolymerization
sol_SOC(i) = par.Dliq*(par.soilM(i)^3)*par.frac*SOC(i);
sol_SON(i) = par.Dliq*(par.soilM(i)^3)*par.frac*SON(i);
    
%this section of code will calculate vmax  Km and CUE.
% Equations for kinetic temperature relationships
%uptake kinetics(base model assumes C and N kinetics are equal)
vmax_UPT_C = par.A_UPT_C .* exp(-par.Ea_UPT_C./(par.R.*(par.T(i)+273))); %temp sensitive according to arrhenius
vmax_upt_N = vmax_UPT_C;
km_upt_N = par.km_UPT_C;

%depolymerization kinetics (base model assumes C and N kinetics are equal)
Vmax_C = par.Vmax_0 * exp(-par.Ea_up./(par.R.*(par.T(i) + 273)));
Vmax_N = Vmax_C; 
Km_N =  par.Km_C;

%This section of code calculates the changes in pool sizes over model time
%using a series of differential equations

            UPT_C(i) = MIC_C(i) .* vmax_UPT_C *(DOC(i) ./ (par.km_UPT_C + DOC(i)))*O2(i)/(par.Km_O2 + O2(i)); %microbial C uptake,michaelis-menton dynamics
            CMIN(i) =  UPT_C(i) .* (1-par.CUE); %C mineralization
            
            UPT_N(i) = MIC_N(i) .* vmax_upt_N *(DON(i) ./ (km_upt_N + DON(i)))*O2(i)/(par.Km_O2 + O2(i)); %microbial N uptake michaelis-menton dynamics
            DEATH_C(i) = par.r_death .* MIC_C(i); %microbial C turnover, first order process
            DEATH_N(i) = par.r_death .* MIC_N(i); %microbial N turnover,first order process
            
      %Resource Allocation
            %enzyme production
            Enz_C(i) = par.p.*(par.CUE.*UPT_C(i));%amount of C available for enzyme production after C allocated to C mineralization
            Enz_N(i) = par.q.*UPT_N(i);%amount of N available for enzyme production
            if (Enz_C(i)/par.CN_enz) >=Enz_N(i) %enz production N limited case
                EPROD(i) = Enz_N(i);%number of enzymes produced = N cost for enzyme production(1 enzyme made for every 1 N used)
            else %C limited case
                EPROD(i) = Enz_C(i)/par.CN_enz;%enz production C limited
            end
            %Amount of C and N available for growth dependent on how much is left over after enz production 
            Growth_C(i) = (1-par.p).*(UPT_C(i).*par.CUE)+ Enz_C(i) - par.CN_enz*EPROD(i);%Available C for biomass growth
            Growth_N(i) = (1-par.q).*UPT_N(i)+ Enz_N(i) - EPROD(i);%available N for biomass growth
            %micrbial growth
            if Growth_C(i)/par.CN_m >= Growth_N(i) %microbes N limited case
                Growth(i) = Growth_N(i);%amount of microbial N biomass equals cost of N to produce microbes 
            else
                Growth(i) = Growth_C(i)/par.CN_m;%microbes C limited case
                
            end
            %after enz prod and growth, N or C left over gets mineralized
            overflow_C(i) = Growth_C(i) - par.CN_m* Growth(i);%extra C after microbes produced goes to overflow metabolism
            NMIN(i) = Growth_N(i) - Growth(i);% N mineralization
            %Microbial biomass pools
            MIC_C(i+1)  =  MIC_C(i) + par.dt*(par.CN_m*Growth(i) - DEATH_C(i));%Microbial biomass C pool, Growth multiplied by C:N of microbes because Growth = N cost to build microbe, but C cost is greater
            MIC_N(i+1) = MIC_N(i) + par.dt*(Growth(i) - DEATH_N(i));%Microbial biomass N pool
            
            %Enzyme Pool & Turnover
            ELOSS(i) = par.r_ECloss * EC(i);%enzyme turnover
            EC(i+1) = EC(i) +  par.dt * (EPROD(i) - ELOSS(i));%enzyme pool
            
            %Depolymerization inputs,derived from Allison et al 2010.
            DECOM_C(i) = Vmax_C .* par.a*EC(i) .*sol_SOC(i) ./(par.Km_C + sol_SOC(i)); %depolymerization of SOC by enzymes
            DECOM_N(i) = Vmax_N.*(1-par.a)*EC(i).*sol_SON(i)./(Km_N + sol_SON(i)); %depolymerization of SON by enzymes
            
            %SOM pools
            SOC(i+1) = SOC(i) + par.dt * (par.Litter_C + DEATH_C(i) * par.MIC_to_SOC - DECOM_C(i));
            SON(i+1) = SON(i) + par.dt * (par.Litter_N + DEATH_N(i) * par.MIC_to_SON - DECOM_N(i));
 
            %Dissolved C&N pools
            DOC(i+1) = DOC(i) + par.dt * (par.DOC_input + par.E + DECOM_C(i) + DEATH_C(i)*(1-par.MIC_to_SOC) + (par.CN_enz/(1+par.CN_enz)).*ELOSS(i) - UPT_C(i));
            DON(i+1) = DON(i) + par.dt * (par.DON_input+ par.E/par.CN_ex + DECOM_N(i) + DEATH_N(i) * (1-par.MIC_to_SON) + (1/par.CN_enz).*ELOSS(i)- UPT_N(i));
            %turnover from enzymes deposited here. Because enzymes are not split up into separate pools within the model, I calculate the
            %amount of C and N created by the turnover of an enzyme by using the C:N of the enzymes

end

%annual estimates at end of run
soc = sum(SOC(par.Nt-4561:par.Nt,1));
son = sum(SON(par.Nt-4561:par.Nt,1));
doc = sum(DOC(par.Nt-4561:par.Nt,1));
don = sum(DON(par.Nt-4561:par.Nt,1));
mbc = sum(MIC_C(par.Nt-4561:par.Nt,1));
mbn = sum(MIC_N(par.Nt-4561:par.Nt,1));
ec = sum(EC(par.Nt-4561:par.Nt,1));
nmin = sum(NMIN(par.Nt-4561:par.Nt,1));
resp = sum(CMIN(par.Nt-4561:par.Nt,1));%in mg/cm3/yr

out = [soc son doc don mbc mbn ec nmin resp];

end

