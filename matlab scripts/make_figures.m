%make damm-mcnip figures script

load('/Users/rzabramoff/Documents/MATLAB/microbial_model_trunk/safe_R1.1/vbsa_damm_mcnip_ECA_N3000.mat')
%plot the main effects of a sensitivity analyis on N=3000 samples
X_Labnew = {'kmDep','eaDep','O2frac','Cfrac','dLiq','dGas','kmO2','p','q','death',...
    'micToSom','aUpt','kmUpt','CUE','eaUpt','a','enzLoss','aDep'};
sensArr = table(num2cell(Si)', X_Labnew');
newp = (sortrows(sensArr,-1));
Si = cell2mat(table2cell(newp(:,1)))';
X_Labels = table2cell(newp(:,2));

figure
boxplot1(Si,X_Labels)
set(gca,'FontSize',20,'XTickLabelRotation',45)

load('/Users/rzabramoff/Documents/MATLAB/microbial_model_trunk/safe_R1.1/glue_damm_mcnip_ECA_N3000.mat')
figure
scatter_plots(X(:,2),Y,[],'RMSE',X_Labels(1,2))
xlabel('eaDep (kJ mol^-^1)')
ylabel('RMSE')
%scatter_plots(X,Y,[],'RMSE',X_Labels,idx)

load('/Users/rzabramoff/Documents/MATLAB/microbial_model_trunk/safe_R1.1/vbsa_damm_mcnip_ECA_dropEa_N3000.mat')
%plot the total effects of a sensitivity analyis on N=3000 samples
X_Labnew = {'kmcn','acn','o2airfrac','frac','dliq','dgas','kmo2','p','q','rdeath',...
    'mictosoc','aupt','kmupt','cue','eaupt','a','renzloss'}; %if cut eaDep
sensArr = table(num2cell(STi)', X_Labnew');
newp = (sortrows(sensArr,-1));
STi = cell2mat(table2cell(newp(:,1)))';
X_Labels = table2cell(newp(:,2));

figure
boxplot1(STi,X_Labels)
set(gca,'FontSize',20,'XTickLabelRotation',45)

