%make some DAMM figures! har har..
%1 year runs for everybody
%load('vars_DMCECA.mat')

%headers = ['DMC'];
%heads = cellstr(headers);
%csvwrite_with_headers('DMCECA_cmin.csv',CMIN(Nt-6072:Nt-1512),heads)

%Col 5 = SoilT
%Col 6 = SoilM
%Col 8 = fractional DOY
%Col 9 = Flux data
%Col 11 = DAMM
%Col 12 = Linear Reg
brewcolor1 = [31/255 120/255 180/255];
brewcolor2 = [178/255 223/255 138/255];
brewcolor3 = [0 0 0];
brewcolor4 = [0.5 0.5 0.5];

Fluxes = xlsread('../data & excel files/mydataFluxP.xlsx');
doytime = linspace(113,303,4561);

%timeseries fluxes
figure
plot(doytime,Fluxes(:,9).*1000/(10000*10),'LineStyle','none',...
    'Marker','.','Color',brewcolor4,'MarkerSize',10); 
hold on;
%plot(doytime,CMIN.*1000,'LineStyle','none',...
%    'Marker','.','Color',brewcolor1,'MarkerSize',10)
plot(doytime,CMIN(Nt-6072:Nt-1512).*1000,'LineStyle','none',...
    'Marker','.','Color',brewcolor1,'MarkerSize',10)
hold on;
plot(doytime,Fluxes(:,11).*1000/(10000*10),'LineStyle','none',...
    'Marker','.','Color',brewcolor2,'MarkerSize',10); 
hold on;
plot(doytime,Fluxes(:,12).*1000/(10000*10),'LineStyle','none',...
    'Marker','.','Color',brewcolor3,'MarkerSize',10); 
xh = xlabel('Day of Year');
yh = ylabel('C efflux (gC cm^-^3 hr^-^1)');
set([xh,yh],'fontweight','bold');
set(gca,'FontSize',15)
legend({'Data','DAMM-MCNiP','DAMM','Regression'},'Location','northeast','FontSize',15);
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
xlim([100 315]);

%regressions
seq = linspace(0,3.5);
figure
data = Fluxes(:,9).*1000/(10000*10);
plot(0:4,0:4,'Color',brewcolor4,'LineStyle','--'); %1:1
hold on;
%plot(seq,seq*0.66935+0.32771,'Color',brewcolor1,'LineStyle','--'); %damm-mcnip
%hold on;
 plot(seq,seq*0.44871+0.62835,'Color',brewcolor1,'LineStyle','--'); %damm-mcnip
 hold on;
plot(seq,seq*0.85022+0.22646,'Color',brewcolor2,'LineStyle','--'); %damm
hold on;
plot(seq,seq*0.557823+0.354917,'Color',brewcolor3,'LineStyle','--'); %empirical
hold on;
%plot(data,CMIN.*1000,'LineStyle','none',...
%    'Marker','.','Color',brewcolor1,'MarkerSize',10)
 plot(data,CMIN(Nt-6072:Nt-1512).*1000,'LineStyle','none',...
     'Marker','.','Color',brewcolor1,'MarkerSize',10)
hold on;
plot(data,Fluxes(:,11).*1000/(10000*10),'LineStyle','none',...
    'Marker','.','Color',brewcolor2,'MarkerSize',10); 
hold on;
plot(data,Fluxes(:,12).*1000/(10000*10),'LineStyle','none',...
    'Marker','.','Color',brewcolor3,'MarkerSize',10); 
xh = xlabel('Observed C efflux (gC cm^-^3 hr^-^1)');
yh = ylabel('Predicted C efflux (gC cm^-^3 hr^-^1)');
set([xh,yh],'fontweight','bold');
set(gca,'FontSize',15)
legend({'1:1','DAMM-MCNiP','DAMM','Regression'},'Location','northeast','FontSize',15);
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
xlim([0 3.5]);
ylim([0 3.5]);

%regressions doy = 100 - 260
figure
data = Fluxes(1:4561-1032,9).*1000/(10000*10);
plot(0:4,0:4,'Color',brewcolor4,'LineStyle','--'); %1:1
hold on;
plot(seq,seq*0.63587+0.38922,'Color',brewcolor1,'LineStyle','--'); %damm-mcnip
hold on;
plot(seq,seq*0.84658+0.25077,'Color',brewcolor2,'LineStyle','--'); %damm
hold on;
plot(seq,seq*0.54153+0.38335,'Color',brewcolor3,'LineStyle','--'); %empirical
hold on;
plot(data,CMIN(Nt-6072:Nt-2544).*1000,'LineStyle','none',...
    'Marker','.','Color',brewcolor1,'MarkerSize',10)
hold on;
plot(data,Fluxes(1:4561-1032,11).*1000/(10000*10),'LineStyle','none',...
    'Marker','.','Color',brewcolor2,'MarkerSize',10); 
hold on;
plot(data,Fluxes(1:4561-1032,12).*1000/(10000*10),'LineStyle','none',...
    'Marker','.','Color',brewcolor3,'MarkerSize',10); 
xh = xlabel('Observed C efflux (gC cm^-^3 hr^-^1)');
yh = ylabel('Predicted C efflux (gC cm^-^3 hr^-^1)');
set([xh,yh],'fontweight','bold');
set(gca,'FontSize',15)
%legend({'1:1','DAMM-MCNiP','DAMM','Regression'},'Location','northeast','FontSize',15);
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
ylim([0 3.5]);

%just dmc
figure
plot(doytime,CMIN(Nt-6072:Nt-1512).*1000,'LineStyle','none',...
    'Marker','.','Color',brewcolor1,'MarkerSize',10); 
hold on;
plot(doytime,Fluxes(:,9).*1000/(10000*10),'LineStyle','none',...
    'Marker','.','Color',[0.5 0.5 0.5],'MarkerSize',10); 
xh = xlabel('Day of Year');
yh = ylabel('C efflux (gC cm^-^3 hr^-^1)');
set([xh,yh],'fontweight','bold');
set(gca,'FontSize',15)
legend({'Model','Data'},'Location','northeast','FontSize',15);
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
 
%sum fluxes to annual, but need to remove those rows where data is NA
sum(Fluxes(:,9).*1000/(10000*10)); %data
sum(CMIN(Nt-4560:Nt).*1000);
sum(Fluxes(:,11).*1000/(10000*10));
sum(Fluxes(:,12).*1000/(10000*10));

%timeseries fluxes close-up %make separate graph for temp and soilm
figure
plot(doytime,Fluxes(:,9).*1000/(10000*10),'LineStyle','-','Color',brewcolor4); hold on;
plot(doytime,CMIN(Nt-6072:Nt-1512).*1000,'LineStyle','-','Color',brewcolor1); hold on;
plot(doytime,Fluxes(:,11).*1000/(10000*10),'LineStyle','-','Color',brewcolor2); hold on;
plot(doytime,Fluxes(:,12).*1000/(10000*10),'LineStyle','-','Color',brewcolor3);
xh = xlabel('Day of Year'); yh = ylabel('C efflux (gC cm^-^3 hr^-^1)');
set([xh,yh],'fontweight','bold'); set(gca,'FontSize',15);
legend({'Data','DAMM-MCNiP','DAMM','Regression'},'Location','northeast','FontSize',15);
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
%xlim([125 135]);
%xlim([165 175]);
%xlim([230 250]);
%xlim([160 180]);
%xlim([100 260]);
xlim([100 315]);
%rectangle('Position',[0 260 315-260 3.5],'FaceColor',[0.5 0.5 0.5])
area([260 315],[3.5 3.5],'FaceColor',[0.5 0.5 0.5],'FaceAlpha',.3,'EdgeAlpha',.3)

%temperature and moisture
figure
[ax,p1,p2] = plotyy(doytime,Fluxes(:,5),doytime,Fluxes(:,6));
xh = xlabel(ax(1),'Day of Year'); % label x-axis
yh1 = ylabel(ax(1),'Temperature (^oC)'); % label left y-axis
yh2 =ylabel(ax(2),'Volumetric Moisture (cm^-^3 H_2O cm^-^3 soil)'); % label right y-axis
p1.LineWidth = 2;
p2.LineWidth = 2;
p1.Color = 'r';
p2.Color = 'b';
set([xh,yh1,yh2],'fontweight','bold');
set([ax,xh,yh1,yh2],'FontSize',15);
set([xh,yh1,yh2],'Color','k');
set(ax(1),'YColor','k');
set(ax(2),'YColor','k');
legend({'Soil Temperature 10 cm','Soil Moisture 2-8 cm'},'Location','northeast','FontSize',15);
%set(gcf,'PaperSize',[11 9]);
ylim(ax(1),[0 20]);
ylim(ax(2),[0 1.5]);
xlim(ax(1),[125 135]);
xlim(ax(2),[125 135]);

