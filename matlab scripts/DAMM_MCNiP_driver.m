%DAMM_MCNiP driver

Yr = 1;
scal = [1 1.6708 1.7177]; %mean 0.356 for all..  
imp = 0;
n = [2009 2013 2014];
A =  2.43E10; %1.0815E11; %5.38*10^10;
Ea = 59.19; % 61.77; %72.26; 
frac = 0.000414; %0.001; %
cue = 0.31;
amp = 0.0005; %0.001; %aka 230 gC/m2 for 190 days
            %0.0005 aka 115 gC/m2 
              %183 gC/m2 (oak litterfall)
refin = 0.001;

for i = 1:length(n)
DAMM_MCNiP_ECA_fx(Yr, scal(i), imp, num2str(n(i)), A, Ea, frac, cue, amp, refin)
end

for i = 1:length(n)
DAMM_MCNiP_fx(Yr, scal(i), imp, num2str(n(i)), A, Ea, frac, cue, amp, refin)
end

dmc = xlsread('/Users/rzabramoff/Dropbox (Climate)/damm-mcnip/data & excel files/mydataFluxRight.xlsx');
kath13 = xlsread('/Users/rzabramoff/Dropbox (Climate)/damm-mcnip/data & excel files/kath13.xlsx');
kath14 = xlsread('/Users/rzabramoff/Dropbox (Climate)/damm-mcnip/data & excel files/kath14.xlsx');

damm09 = dmc(:,13);
predict09 = dmc(:,16);
damm13 = kath13(:,10);
predict13 = kath13(:,11);
damm14 = kath14(:,10);
predict14 = kath14(:,11);

figure
subplot(3,2,1);
Nt = 8759;
xa = [1:8759]./24;
xb = [Nt-6072:Nt-1512]./24;
xc = [2592:2592+(4502-1)]./24;
xd = [2904:2904+(4616-1)]./24;
load('vars_DMCECA_2009')
        plot(xa,CMIN,'LineWidth',3)
        hold on;
        plot(xb,Fluxes(:,9)./(10000*10));
        hold on;
        plot(xb,damm09./(10000*10));
        hold on;
        plot(xb,predict09./(10000*10));
        legend('dmc','data','damm','reg')
        title('ECA 2009')
        xlabel('Day of Year')
        ylabel('mg C cm^-^3 soil hr^-^1')
        xlim([(Nt-6072)/24 (Nt-1512)/24])
subplot(3,2,3);
load('vars_DMCECA_2013')
        plot(xc,CMIN,'LineWidth',3)
        hold on;
        plot(xc,Fluxes./(10000*10));
        hold on;
        plot(xc,damm13./(10000*10));
        hold on;
        plot(xc,predict13./(10000*10));
        title('2013')
        xlabel('Day of Year')
        ylabel('mg C cm^-^3 soil hr^-^1')
subplot(3,2,5);
load('vars_DMCECA_2014')
        plot(xd,CMIN,'LineWidth',3)
        hold on;
        plot(xd,Fluxes./(10000*10));
        hold on;
        plot(xd,damm14./(10000*10));
        hold on;
        plot(xd,predict14./(10000*10));
        title('2014')
        xlabel('Day of Year')
        ylabel('mg C cm^-^3 soil hr^-^1')
subplot(3,2,2);
load('vars_DMC_2009')
        plot(xa,CMIN,'LineWidth',3)
        hold on;
        plot(xb,Fluxes(:,9)./(10000*10));
        hold on;
        plot(xb,damm09./(10000*10));
        hold on;
        plot(xb,predict09./(10000*10));
        title('MM 2009')
        xlabel('Day of Year')
        ylabel('mg C cm^-^3 soil hr^-^1')
        xlim([(Nt-6072)/24 (Nt-1512)/24])
subplot(3,2,4);
load('vars_DMC_2013')
        plot(xc,CMIN,'LineWidth',3)
        hold on;
        plot(xc,Fluxes./(10000*10));
        hold on;
        plot(xc,damm13./(10000*10));
        hold on;
        plot(xc,predict13./(10000*10));
        title('2013')
        xlabel('Day of Year')
        ylabel('mg C cm^-^3 soil hr^-^1')
subplot(3,2,6);
load('vars_DMC_2014')
        plot(xd,CMIN,'LineWidth',3)
        hold on;
        plot(xd,Fluxes./(10000*10));
        hold on;
        plot(xd,damm14./(10000*10));
        hold on;
        plot(xd,predict14./(10000*10));
        title('2014')
        xlabel('Day of Year')
        ylabel('mg C cm^-^3 soil hr^-^1')

% headers = ['DMC'];
% heads = cellstr(headers);
% %csvwrite_with_headers('DMCECA_var_cmin.csv',CMIN(Nt-6072:Nt-1512),heads)
% csvwrite_with_headers('DMCMM_var_cmin.csv',CMIN,heads)
