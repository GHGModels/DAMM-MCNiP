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
subplot(2,2,1);
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
        %hold on;
        %plot(xb,predict09./(10000*10));
        legend('dmc','data','damm')
        title('ECA 2009')
        xlabel('Day of Year')
        ylabel('mg C cm^-^3 soil hr^-^1')
        xlim([(Nt-6072)/24 (Nt-1512)/24])
subplot(2,2,2);
load('vars_DMCECA_2013')
        plot(xc,CMIN,'LineWidth',3)
        hold on;
        plot(xc,Fluxes./(10000*10));
        hold on;
        plot(xc,damm13./(10000*10));
        %hold on;
        %plot(xc,predict13./(10000*10));
        title('2013')
        xlabel('Day of Year')
        ylabel('mg C cm^-^3 soil hr^-^1')
subplot(2,2,3);
load('vars_DMCECA_2014')
        plot(xd,CMIN,'LineWidth',3)
        hold on;
        plot(xd,Fluxes./(10000*10));
        hold on;
        plot(xd,damm14./(10000*10));
        %hold on;
        %plot(xd,predict14./(10000*10));
        title('2014')
        xlabel('Day of Year')
        ylabel('mg C cm^-^3 soil hr^-^1')
        
        
        
 headers = ['DMC'];
 heads = cellstr(headers);
 load('vars_DMCECA_2009')
 csvwrite_with_headers('DMCECA_2009.csv',CMIN(Nt-6072:Nt-1512),heads)
 load('vars_DMCECA_2013')
 csvwrite_with_headers('DMCECA_2013.csv',CMIN,heads)
 load('vars_DMCECA_2014')
 csvwrite_with_headers('DMCECA_2014.csv',CMIN,heads)