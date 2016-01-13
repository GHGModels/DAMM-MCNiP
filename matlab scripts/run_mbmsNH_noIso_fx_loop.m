function run_mbmsNH_noIso_fx_loop(cni,reps,continueOn,dormOn,mref,Tref,surfinit)

[status,results]=system('pwd');
sstrs=strsplit(results,'/microbial_model_trunk');

matfldir=[sstrs{1},'/microbial_model_trunk/C_plus_NH/data'];

system(['mkdir -p ', matfldir]);

     if continueOn == 1
     iofile=[matfldir,'/mbmsNH_noIso.mat'];   
     else
     iofile=[matfldir,'/mbmsNH_noIso.mat'];   
         %iofile=[matfldir,'/mbmsNH_noIso_',num2str(reps),'reps_mref',num2str(mref),'_Tref',num2str(Tref),'_surfinit',num2str(surfinit),'.mat'];
     end

% Set number of pools
n_polymers=1;   % cellulose-hemicellulose and ligin
n_monomers=1;   % dissolved organic carbon (compounds of varying complexity)
n_enzymes=1;    % use a generic enzyme to minimize model parameters (BG/Phenox activities reflect corresponding degradation rates) 
n_microbep=1;   % microbes can directly uptake DOC 
n_micc=n_microbep;
n_surfaces=1;
n_co2=1;
n_enzymes_ads=n_enzymes*n_surfaces; % consider every surface has all enzymes/monomers (e.g., 3enzymes x 2surfaces network: enz1-sur1,enz1-sur2,enz2-sur1,...,enz3-sur2)
n_monomers_ads=n_monomers*n_surfaces;
n_all=n_microbep+n_micc+n_enzymes+n_polymers+n_monomers+n_surfaces+n_co2+n_enzymes_ads+n_monomers_ads;

% Set key parameters those results in microbe biomass of 3.89% total organic carbon
cpar0=[300/365*.8 300/365*.2 4d-4*57.653333 1d-1*0.536542 2d-2*546.713826 ...
    1.5d-5*875.965050 0.002215393623082 0.102525103775396 1d-5*611.427356 2d-2*120.666087 ...
    .001 .006 .01 .006 5d-6*386.294367]';
global cpar;            % a vector of considered calibrating parameters
cpar=[
    {repmat(cpar0(1)/n_polymers,1,n_polymers)};       % 1 x npolymers, polymer input, 1/day
    {repmat(cpar0(2)/n_monomers,1,n_monomers)};       % 1 x nmonomers, monomer input, 1/day
    {cpar0(3)};       % 1 x nmicrobes, microbial maintenance rate (1/day)
    {cpar0(4)};       % 1 x nmicrobes, reserve turnover rate, 1/day
    {repmat(cpar0(5),1,n_monomers)};       % 1 x nmonomers, maximum doc uptake rate                   (1/day)
    
    {cpar0(6)};       % 1 x nmicrobes, microbial death rate       (1/day)
    {cpar0(7)};       % 1 x nenzymes, maximum enzyme production rate     (1/day) [inverted from analytic solution]
    {cpar0(8)};       % 1 x nmicrobes, maximum microbial growth rate (1/day) [inverted from analytic solution]
    {cpar0(9)};       % 1 x nenzymes, enzyme decay rate          (1/day)
    {repmat(cpar0(10),1,n_polymers)};      % 1 x npolymers, maximum som degradation rate              (1/day)
    
    {repmat(cpar0(11),1,n_surfaces)};      % 1 x nsurfaces, maximum enzyme adsorption rate (1/day)
    {repmat(cpar0(12),1,n_surfaces)};      % 1 x nsurfaces, decay rate of adsorbed enzyme (1/day)
    {repmat(cpar0(13),1,n_monomers)};      % 1 x nmonomers, maximum monomer adsorption rate (1/day)
    {repmat(cpar0(14),1,n_monomers)};      % 1 x nmonomers, decay rate of adsorbed monomer (1/day)
    {cpar0(15)};      % 1 x nmicrobes, enzyme production rate used for analytic solution only (1/day)
    ];

%
global vid;             % id index
global par_mic;         % each par represents a microbe
global par_enz;         % each par represents a enzyme
global par_surface;     % each par represents a mineral surface
global par_ss;          % parameter structure for substrate
global par_water;       % paraters for water module
global input;           % substrate input structure

% Set id index - C cycling
vid.microbep=1:n_microbep; 
vid.micc=vid.microbep(end)+(1:n_microbep);
vid.surfaces=vid.micc(end)+(1:n_surfaces);
vid.monomers=vid.surfaces(end)+(1:n_monomers);
vid.monomers_ads=vid.monomers(end)+(1:n_monomers_ads);
vid.polymers=vid.monomers_ads(end)+(1:n_polymers);
vid.enzymes=vid.polymers(end)+(1:n_enzymes);
vid.enzymes_ads=vid.enzymes(end)+(1:n_enzymes_ads);
vid.co2=vid.enzymes_ads(end)+(1:n_co2);
vid.cue=vid.co2(end)+(1:n_microbep); %track cue
% Set id index - NH cycling
vid.miccNH=vid.cue(end)+(1:n_microbep);
vid.monomersNH=vid.miccNH(end)+(1:n_monomers);
vid.monomers_adsNH=vid.monomersNH(end)+(1:n_monomers_ads);
vid.polymersNH=vid.monomers_adsNH(end)+(1:n_polymers);
vid.nmin=vid.polymersNH(end)+(1:n_co2);
vid.nue=vid.nmin(end)+(1:n_microbep); 
% Set id index - dormant pools
vid.microbed=vid.nue(end)+(1:n_microbep); 
vid.miccd=vid.microbed(end)+(1:n_microbep); 
vid.miccNd=vid.miccd(end)+(1:n_microbep); 

% Set external input (should change to time series)
input.polymers=cpar{1}.*ones(1,n_polymers);
input.monomers=cpar{2}.*ones(1,n_monomers);
% Set external input NH
input.polymersNH=input.polymers./cni;
input.monomersNH=input.monomers./cni;

% if continueOn == 1
%     %load mbmsNH_noIso_spinup_200_reps_CN50_soilmDeath_2of2.mat
%     load mbmsNH_noIso_spinup_200_reps_CN50_2of2.mat
%     %just make x a list of YOUT_ctrl(end,:)...
%     x = YOUT_ctrl(end,:);
% else
%     if exist('surfinit')
%        surfinit=surfinit;
%     else
%         surfinit=1000;
%     end
% Set initial states
x(vid.microbep)=repmat(20/n_microbep,1,n_microbep);
x(vid.micc)=repmat(5/n_microbep,1,n_microbep);
x(vid.surfaces)=repmat(surfinit/n_surfaces,1,n_surfaces);
x(vid.monomers)=repmat(20/n_monomers,1,n_monomers);
x(vid.monomers_ads)=reshape(repmat(.2*x(vid.monomers)/n_surfaces,n_surfaces,1),n_surfaces*n_monomers,1);
x(vid.polymers)=repmat(200/n_polymers,1,n_polymers);
x(vid.enzymes)=repmat(.1/n_enzymes,1,n_enzymes);
x(vid.enzymes_ads)=reshape(repmat(.2*x(vid.enzymes)/n_surfaces,n_surfaces,1),n_surfaces*n_enzymes,1);
x(vid.co2)=0;
x(vid.cue)=repmat(0,1, n_microbep); %set initial cue
% Set initial states NH
x(vid.miccNH)=repmat(0.5/n_microbep,1,n_microbep);
x(vid.monomersNH)=repmat(1/n_monomers,1,n_monomers);
x(vid.monomers_adsNH)=reshape(repmat(0.01*x(vid.monomers)/n_surfaces,n_surfaces,1),n_surfaces*n_monomers,1);
x(vid.polymersNH)=repmat(10/n_polymers,1,n_polymers);
x(vid.nmin)=0;
x(vid.nue)=repmat(0,1, n_microbep); %set initial nue
% Set initial states - dormant pools
x(vid.microbed)=repmat(20/n_microbep,1,n_microbep);
x(vid.miccd)=repmat(5/n_microbep,1,n_microbep); 
x(vid.miccNd)=repmat(0.5/n_microbep,1,n_microbep);
%end

% Set key parameters (same for each function group of microbe/surface, might change to consider different groups)
par_mic=set_microbe_par_default();          % microbe-related parameters
for i=2:n_microbep 
    par_mic(i)=set_microbe_par_default(); 
end;    

par_enz=set_enzyme_par_default();           % enzyme-related parameters
for i=2:n_enzymes
    par_enz(i)=set_enzyme_par_default();
end
par_surface=set_msurface_par_default();     % mineral surface-related parameters
for i=2:n_surfaces
    par_surface(i)=set_msurface_par_default();
end;

par_ss=set_substrate_par_default();         % substrate-related parameters
%par_mic.kappa_micb = par_mic.kappa_micb*par_ss.soilm; %added 01/014/2016
%(scenario: soil moisture increases reserve yield to metabolism)
%par_mic.mr_micb = par_mic.mr_micb*2/par_ss.soilm; %added 01/014/2016
%(scenario: drying increases maintenance demand imposing a C cost to producing osmolytes)

par_water=set_water_par_default();         % water-related parameters
par_water.manzoni = dormOn;

% % Make a copy of those temperature dependent parameters
par_surface_ref=par_surface;
par_mic_ref=par_mic;
par_enz_ref=par_enz;
%par_mic.kappa_micb=1d10;

% Define physiological temperature response curve
% if exist('Tref')
%    Tref=Tref;
% else
%     Tref=290;          %reference temperature
% end
    Tbase=290;
    [temp0,T_fact00]=get_microbe_physiology_Tcurve(Tref);   % active enzyme fraction in total enzyme vs temperaure

% Define activation energyies
TEa=set_Ea_default();

% %Set up seasonal cycle params
% if mref > -1
%    mref=mref;
% else
%     mref=0.49;
% end
                A1=0.3;
                w1=2.*pi/365;

% Set up and run the model
%     if exist('tempin')
%     yrs = reps*(length(tempin)-1)/365;
%     tempi = repmat(tempin,reps);
%     else
    yrs = reps;
%     end
%     
%     if exist('moisin')
%     yrs = reps*(length(moisin)-1)/365;
%     moisi = repmat(moisin,reps);
%     else
%     yrs = reps;
%     end
%     
%     if exist('dt')
%         dt=dt;
%     else
        dt=1; % days
%     end

tend= 365*yrs;
kend=fix(tend/dt);
TOUT_ctrl=zeros(kend+1,1);
YOUT_ctrl=zeros(kend+1,length(x)); %TEMP solution
YOUT_ctrl(1,:)=x;
TEMP=zeros(kend,1);
for kk = 1 : kend

%     % -- Incorporate temperature effects on parameters
    t=(kk-0.5)*dt;

%     if exist('tempi')
%     temp= tempi(kk) + 273.15;
%     else
    temp=Tref+10.*sin(w1.*t-pi/2);
%     end
%     if exist('moisin')
%     par_ss.soilm=moisi(kk);
%     for i=1:n_enzymes
%     par_enz(i).soilm=moisi(kk);
%     end
%     else
    par_ss.soilm=mref;%+A1.*sin(w1.*t-1.5.*pi);
    for i=1:n_enzymes
    par_enz(i).soilm=mref;%+A1.*sin(w1.*t-pi/2);
    end
%     end
                %par_ss.soilm=mref+A1.*sin(w1.*t-pi/2); %unimodal monomers
                %par_ss.soilm=mref+A1.*sin(w1.*t-1.5.*pi); %mediterranean (reverse)
                %par_ss.hd=mref+A1.*sin(w1.*t-pi/2); %unimodal monomers
                %par_ss.hd=mref+A1.*sin(w1.*t-1.5.*pi); %mediterranean (reverse)
              %for i=1:n_enzymes
                %par_enz(i).soilm=mref+A1.*sin(w1.*t-pi/2); %unimodal enzymes
                %par_enz(i).soilm=mref+A1.*sin(w1.*t-1.5.*pi); %mediterranean (reverse)
                %par_enz(i).he=mref+A1.*sin(w1.*t-pi/2); %unimodal enzymes
                %par_enz(i).he=mref+A1.*sin(w1.*t-1.5.*pi); %mediterranean (reverse)
              %end
              %input.polymers=0.5; %2xmonmer pulse%temp=temp+4;
    TEMP(kk)=temp;
    T_fact=interp1(temp0,T_fact00,temp);    % get active enzyme fraction at given temperature
    fref0=temp/Tref;    % non-equilibrium but not related with enzyme activities
    fref=T_fact*(temp/Tref);  % non-equilibrium processes which are also related with enzyme activities
    Tinv=1/temp-1/Tref;
    
    for i=1:n_microbep
        % equilibrium processes - Arrhenius equation
        par_mic(i).Kaff_monomer=par_mic_ref(i).Kaff_monomer.*exp(-TEa.Ea_Kaff_monomer_micb(i,:)*Tinv);
        par_mic(i).mr_micb=par_mic_ref(i).mr_micb.*exp(-TEa.Ea_mr_micb(i)*Tinv);
        
        % non-equilibrium and enzyme processes - Eyring's transition state theory and Murphy's equations
        par_mic(i).Vmax_micb=par_mic_ref(i).Vmax_micb*fref.*exp(-TEa.Ea_vmax_micb(i,:)*Tinv); % maxium DOC uptake rate by microbe
        par_mic(i).kappa_micb=par_mic_ref(i).kappa_micb*fref.*exp(-TEa.Ea_kappa_micb(i)*Tinv); % reserve turnover rate 
    end
    for i=1:n_enzymes
        % equilibrium processes
        par_enz(i).Kaff_ee_polymer=par_enz_ref(i).Kaff_ee_polymer.*exp(-TEa.Ea_Kaff_ee_polymer_micb(i,:)*Tinv);
        par_enz(i).Kaff_ee_msurf=par_enz_ref(i).Kaff_ee_msurf.*exp(-TEa.Ea_Kaff_ee_msurf(i,:)*Tinv);
        
        % non-equilibrium and enzyme processes
        par_enz(i).Vmax_ee=par_enz_ref(i).Vmax_ee*fref.*exp(-TEa.Ea_vmax_ee(i,:)*Tinv); % maximum SOC degradation rate by enzyme

        % non-equilibrium and non-enzyme processes (mineral adsorption processes)
        par_enz(i).Vmax_ads_enzyme=par_enz_ref(i).Vmax_ads_enzyme*fref0.*exp(-TEa.Ea_Vmax_ads_enzyme(i,:)*Tinv);
    end
    for i=1:n_surfaces
        % equilibrium processes
        par_surface(i).Kaff_monomer=par_surface_ref(i).Kaff_monomer.*exp(-TEa.Ea_Kaff_monomer_msurf(i,:)*Tinv);

        % non-equilibrium and non-enzyme processes
        par_surface(i).Vmax_ads_monomer=par_surface_ref(i).Vmax_ads_monomer*fref0.*exp(-TEa.Ea_Vmax_ads_monomer(i,:)*Tinv);
    end

%     % -- Run core program
    
    YOUT_ctrl(kk+1,:)=adptmbbks1(@many_bug_many_substrate_noIso_NH,YOUT_ctrl(kk,:),TOUT_ctrl(kk)+dt/2,dt);
    TOUT_ctrl(kk+1)=TOUT_ctrl(kk)+dt;
%     if mod(kk,100)==1 disp(kk); end;

end
        save(iofile,'YOUT_ctrl','TOUT_ctrl','vid');

end

