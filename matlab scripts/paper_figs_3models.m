%these are all ECA

dmc = xlsread('/Users/rzabramoff/Dropbox (Climate)/damm-mcnip/data & excel files/mydataFluxRight.xlsx');
kath13 = xlsread('/Users/rzabramoff/Dropbox (Climate)/damm-mcnip/data & excel files/kath13.xlsx');
kath14 = xlsread('/Users/rzabramoff/Dropbox (Climate)/damm-mcnip/data & excel files/kath14.xlsx');

damm09 = dmc(:,13);
predict09 = dmc(:,16);
damm13 = kath13(:,10);
predict13 = kath13(:,11);
damm14 = kath14(:,10);
predict14 = kath14(:,11);

brewcolor1 = [31/255 120/255 180/255];
brewcolor3 = [178/255 223/255 138/255];
brewcolor2 = [0 0 0];
brewcolor4 = [0.5 0.5 0.5];

%plotting
figure
subplot(2,2,1);
Nt = 8759;
xa = [1:8759]./24;
xb = [Nt-6072:Nt-1512]./24;
xc = [2592:2592+(4502-1)]./24;
xd = [2904:2904+(4616-1)]./24;
load('vars_DMCECA_2009')
        plot(xb,Fluxes(:,9).*1000/(10000*10),'LineStyle','none',...
    'Marker','.','Color',brewcolor4,'MarkerSize',10);
        hold on;        
        plot(xa,CMIN*1000,'LineStyle','none',...
    'Marker','.','Color',brewcolor1,'MarkerSize',10)
        hold on;
        plot(xb,damm09.*1000/(10000*10),'LineStyle','none',...
    'Marker','.','Color',brewcolor3,'MarkerSize',10);
        hold on;
    %    plot(xb,predict09.*1000/(10000*10),'LineStyle','none',...
    %'Marker','.','Color',brewcolor2,'MarkerSize',10);
        legend('dmc','data','damm')
        title('2009')
        xh = xlabel('Day of Year');
        yh = ylabel('C efflux (gC cm^-^3 hr^-^1)');
        set([xh,yh],'fontweight','bold');
        set(gca,'FontSize',15)
        legend({'Data','DAMM-MCNiP','DAMM'},'Location','northeast','FontSize',10);
        set(findall(gca, 'Type', 'Line'),'LineWidth',2);
        xlim([(Nt-6072)/24 (Nt-1512)/24])
        jerryrig = Nt;
        ylim([0 5])
subplot(2,2,2);
load('vars_DMCECA_2013')
        plot(xc,Fluxes.*1000/(10000*10),'LineStyle','none',...
    'Marker','.','Color',brewcolor4,'MarkerSize',10);
        hold on;
        plot(xc,CMIN.*1000,'LineStyle','none',...
    'Marker','.','Color',brewcolor1,'MarkerSize',10)
        hold on;
        plot(xc,damm13.*1000/(10000*10),'LineStyle','none',...
    'Marker','.','Color',brewcolor3,'MarkerSize',10);
        hold on;
    %    plot(xc,predict13.*1000/(10000*10),'LineStyle','none',...
    %'Marker','.','Color',brewcolor2,'MarkerSize',10);
        title('2013')
        xh = xlabel('Day of Year');
        yh = ylabel('C efflux (gC cm^-^3 hr^-^1)');
        set([xh,yh],'fontweight','bold');
        set(gca,'FontSize',15)
        set(findall(gca, 'Type', 'Line'),'LineWidth',2);
        xlim([(jerryrig-6072)/24 (jerryrig-1512)/24])
        ylim([0 5])
subplot(2,2,3);
load('vars_DMCECA_2014')
        plot(xd,Fluxes.*1000/(10000*10),'LineStyle','none',...
    'Marker','.','Color',brewcolor4,'MarkerSize',10);
        hold on;
        plot(xd,CMIN.*1000,'LineStyle','none',...
    'Marker','.','Color',brewcolor1,'MarkerSize',10)
        hold on;
        plot(xd,damm14.*1000/(10000*10),'LineStyle','none',...
    'Marker','.','Color',brewcolor3,'MarkerSize',10);
        hold on;
   %     plot(xd,predict14.*1000/(10000*10),'LineStyle','none',...
   % 'Marker','.','Color',brewcolor2,'MarkerSize',10);
        title('2014')
        xh = xlabel('Day of Year');
        yh = ylabel('C efflux (gC cm^-^3 hr^-^1)');
        set([xh,yh],'fontweight','bold');
        set(gca,'FontSize',15)
        set(findall(gca, 'Type', 'Line'),'LineWidth',2)
        xlim([(jerryrig-6072)/24 (jerryrig-1512)/24])
        ylim([0 5])

        
        
        
        
% %regressions
% seq = linspace(0,3.5);
% figure
% data = Fluxes(:,9).*1000/(10000*10);
% plot(0:4,0:4,'Color',brewcolor4,'LineStyle','--'); %1:1
% hold on;
% %plot(seq,seq*0.66935+0.32771,'Color',brewcolor1,'LineStyle','--'); %damm-mcnip
% %hold on;
%  plot(seq,seq*0.44871+0.62835,'Color',brewcolor1,'LineStyle','--'); %damm-mcnip
%  hold on;
% plot(seq,seq*0.85022+0.22646,'Color',brewcolor2,'LineStyle','--'); %damm
% hold on;
% plot(seq,seq*0.557823+0.354917,'Color',brewcolor3,'LineStyle','--'); %empirical
% hold on;
% %plot(data,CMIN.*1000,'LineStyle','none',...
% %    'Marker','.','Color',brewcolor1,'MarkerSize',10)
%  plot(data,CMIN(Nt-6072:Nt-1512).*1000,'LineStyle','none',...
%      'Marker','.','Color',brewcolor1,'MarkerSize',10)
% hold on;
% plot(data,Fluxes(:,11).*1000/(10000*10),'LineStyle','none',...
%     'Marker','.','Color',brewcolor2,'MarkerSize',10); 
% hold on;
% plot(data,Fluxes(:,12).*1000/(10000*10),'LineStyle','none',...
%     'Marker','.','Color',brewcolor3,'MarkerSize',10); 
% xh = xlabel('Observed C efflux (gC cm^-^3 hr^-^1)');
% yh = ylabel('Predicted C efflux (gC cm^-^3 hr^-^1)');
% set([xh,yh],'fontweight','bold');
% set(gca,'FontSize',15)
% legend({'1:1','DAMM-MCNiP','DAMM','Regression'},'Location','northeast','FontSize',15);
% set(findall(gca, 'Type', 'Line'),'LineWidth',2);
% xlim([0 3.5]);
% ylim([0 3.5]);
%  
% %sum fluxes to annual, but need to remove those rows where data is NA
% sum(Fluxes(:,9).*1000/(10000*10)); %data
% sum(CMIN(Nt-4560:Nt).*1000);
% sum(Fluxes(:,11).*1000/(10000*10));
% sum(Fluxes(:,12).*1000/(10000*10));
% 
% %temperature and moisture
% figure
% [ax,p1,p2] = plotyy(doytime,Fluxes(:,5),doytime,Fluxes(:,6));
% xh = xlabel(ax(1),'Day of Year'); % label x-axis
% yh1 = ylabel(ax(1),'Temperature (^oC)'); % label left y-axis
% yh2 =ylabel(ax(2),'Volumetric Moisture (cm^-^3 H_2O cm^-^3 soil)'); % label right y-axis
% p1.LineWidth = 2;
% p2.LineWidth = 2;
% p1.Color = 'r';
% p2.Color = 'b';
% set([xh,yh1,yh2],'fontweight','bold');
% set([ax,xh,yh1,yh2],'FontSize',15);
% set([xh,yh1,yh2],'Color','k');
% set(ax(1),'YColor','k');
% set(ax(2),'YColor','k');
% legend({'Soil Temperature 10 cm','Soil Moisture 2-8 cm'},'Location','northeast','FontSize',15);
% %set(gcf,'PaperSize',[11 9]);
% ylim(ax(1),[0 20]);
% ylim(ax(2),[0 1.5]);
% xlim(ax(1),[125 135]);
% xlim(ax(2),[125 135]);
% 
