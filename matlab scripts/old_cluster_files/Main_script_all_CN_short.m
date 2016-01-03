%RHIZO 
tic
path = 'C:\Users\rose\Documents\MATLAB\'; %Update this if needed
rootDOC = [0.0009306]; %update if needed
deepB = 1:5;
deep2 = 5; 
nb_loops = 20; %20 loops * 5 depths at a time = 100 depths
for xx = 1:length(rootDOC)
    deep = deepB-deep2;
    DECOMPc_mean_last_row_merged = [];
    min_mean_last_row_merged = [];
    BioC_mean_last_row_merged = [];
    SOC_mean_last_row_merged = [];
    
    filename3 = ['Allison_rhizosphere_Input',int2str(xx),'_all_depths_means.mat'];

    
    for x = 1:nb_loops
        deep = deep+deep2;
        %Compute and save everything
        filename = ['Allison_rhizosphere_Input',int2str(xx),'_depths_',int2str(deep(1)),'-',int2str(deep(end)),'_all_data.mat'];
        filename2 = ['Allison_rhizosphere_Input',int2str(xx),'_depths_',int2str(deep(1)),'-',int2str(deep(end)),'_means.mat'];
        rhizosphere2_CN_short(rootDOC(xx),deep,deep2,filename,filename2,path);
          
        load([path,filename2]);
        DECOMPc_mean_last_row_merged = [DECOMPc_mean_last_row_merged DECOMPc_mean_last_row];
        min_mean_last_row_merged = [min_mean_last_row_merged min_mean_last_row];
        BioC_mean_last_row_merged = [BioC_mean_last_row_merged BioC_mean_last_row];
        SOC_mean_last_row_merged = [SOC_mean_last_row_merged SOC_mean_last_row];
        
        save([path,filename3],'BioC_mean_last_row_merged','min_mean_last_row_merged','DECOMPc_mean_last_row_merged','SOC_mean_last_row_merged');
    end
end

load 'Allison_rhizosphere_Input1_all_depths_means.mat' %load 'Allison_rhizosphere_Input1_all_depths_means.mat'

DECOMPc = DECOMPc_mean_last_row_merged.';
MIN = min_mean_last_row_merged.';
BioC = BioC_mean_last_row_merged.';
SOC = SOC_mean_last_row_merged.';
C = cat(2,DECOMPc,MIN,BioC,SOC);

headers = ['DECrhizo '; 'MINrhizo '; 'BioCrhizo';'SOCrhizo '];
heads = cellstr(headers);

csvwrite_with_headers('rhizo.csv',C,heads)
xlswrite('rhizo.xls',C)

%BULK
path = 'C:\Users\rose\Documents\MATLAB\'; %Update this if needed
DOCexpB = -0.05; %Update this if needed
deepB = 1:5;
deep2 = 5;
nb_loops = 20; %20 loops * 5 depths at a time = 100 depths

for xx = 1:length(DOCexpB)
    deep = deepB-deep2;
    DECOMPc_mean_last_row_merged = [];
    min_mean_last_row_merged = [];
    BioC_mean_last_row_merged = [];
    SOC_mean_last_row_merged = [];
    
    filename3 = ['Allison_bulk_Input',int2str(xx),'_all_depths_means.mat'];
    
    for x = 1:nb_loops
        deep = deep+deep2;
        %Compute and save everything
        filename = ['Allison_bulk_Input',int2str(xx),'_depths_',int2str(deep(1)),'-',int2str(deep(end)),'_all_data.mat'];
        filename2 = ['Allison_bulk_Input',int2str(xx),'_depths_',int2str(deep(1)),'-',int2str(deep(end)),'_means.mat'];
        bulk2_CN_short(DOCexpB(xx),deep,deep2,filename,filename2,path);
        
        %Merge and save final means
        load([path,filename2]);
        DECOMPc_mean_last_row_merged = [DECOMPc_mean_last_row_merged DECOMPc_mean_last_row];
        min_mean_last_row_merged = [min_mean_last_row_merged min_mean_last_row];
        BioC_mean_last_row_merged = [BioC_mean_last_row_merged BioC_mean_last_row];
        SOC_mean_last_row_merged = [SOC_mean_last_row_merged SOC_mean_last_row];
        
        save([path,filename3],'BioC_mean_last_row_merged','min_mean_last_row_merged','DECOMPc_mean_last_row_merged','SOC_mean_last_row_merged');
    end
end

load 'Allison_bulk_Input1_all_depths_means.mat' %load 'Allison_rhizosphere_Input1_all_depths_means.mat'

DECOMPc = DECOMPc_mean_last_row_merged.';
MIN = min_mean_last_row_merged.';
BioC = BioC_mean_last_row_merged.';
SOC = SOC_mean_last_row_merged.';
C = cat(2,DECOMPc,MIN,BioC,SOC);

headers = ['DECbulk  '; 'MINbulk  '; 'BioCbulk ';'SOCbulk  '];
heads = cellstr(headers);

csvwrite_with_headers('bulk.csv',C,heads)
xlswrite('bulk.xls',C)

%import data
bulk = xlsread ('bulk.xls');
rhizo = xlsread ('rhizo.xls');
vol = xlsread ('med_rhizo_volume.xlsx'); 

%DECOMPc, MIN, BioC, SOC
effectDECOMPc=zeros(100,1);
rhizoDECOMPc=zeros(100,1);
bulkDECOMPc=zeros(100,1);
effectMIN=zeros(100,1);
rhizoMIN=zeros(100,1);
bulkMIN=zeros(100,1);
effectBioC=zeros(100,1);
rhizoBioC=zeros(100,1);
bulkBioC=zeros(100,1);
effectSOC=zeros(100,1);
rhizoSOC=zeros(100,1);
bulkSOC=zeros(100,1);

for i=1:100;
rhizoDECOMPc(i) = rhizo(i,1);
bulkDECOMPc(i) = bulk(i,1);
effectDECOMPc(i) = log(rhizoDECOMPc(i)/bulkDECOMPc(i));
rhizoMIN(i) = rhizo(i,2);
bulkMIN(i) = bulk(i,2);
effectMIN(i) = log(rhizoMIN(i)/bulkMIN(i));
rhizoBioC(i) = rhizo(i,4);
bulkBioC(i) = bulk(i,4);
effectBioC(i) = log(rhizoBioC(i)/bulkBioC(i));
rhizoSOC(i) = rhizo(i,5);
bulkSOC(i) = bulk(i,5);
effectSOC(i) = log(rhizoSOC(i)/bulkSOC(i));
end

aDECOMPc = mean(effectDECOMPc(1:15));
amin = mean(effectMIN(1:15));
abioc = mean (effectBioC(1:15));
asoc = mean (effectSOC(1:15));

%E = cat(2,effectDECOMPc,effectMIN,effectCO2,effectBioC,effectEnzC);
a = cat(2,aDECOMPc, amin, abioc, asoc);
headers = ['DECeffect '; 'MINeffect '; 'BioCeffect'; 'SOCeffect '];
heads = cellstr(headers);
csvwrite_with_headers('effect1.csv',a,heads)

load 'Allison_rhizosphere_Input2_all_depths_means.mat'

DECOMPc = DECOMPc_mean_last_row_merged.';
MIN = min_mean_last_row_merged.';
BioC = BioC_mean_last_row_merged.';
SOC = SOC_mean_last_row_merged.';
C = cat(2,DECOMPc,MIN,BioC,SOC);

headers = ['DECrhizo '; 'MINrhizo '; 'BioCrhizo'; 'SOCrhizo '];
heads = cellstr(headers);

csvwrite_with_headers('rhizo2.csv',C,heads)
xlswrite('rhizo2.xls',C)

%import data
bulk = xlsread ('bulk.xls');
rhizo = xlsread ('rhizo2.xls');
vol = xlsread ('med_rhizo_volume.xlsx'); 

%DECOMPc, MIN, CO2, BioC, EnzC
effectDECOMPc=zeros(100,1);
rhizoDECOMPc=zeros(100,1);
bulkDECOMPc=zeros(100,1);
effectMIN=zeros(100,1);
rhizoMIN=zeros(100,1);
bulkMIN=zeros(100,1);
effectBioC=zeros(100,1);
rhizoBioC=zeros(100,1);
bulkBioC=zeros(100,1);
effectSOC=zeros(100,1);
rhizoSOC=zeros(100,1);
bulkSOC=zeros(100,1);

for i=1:100;
rhizoDECOMPc(i) = rhizo(i,1);
bulkDECOMPc(i) = bulk(i,1);
effectDECOMPc(i) = log(rhizoDECOMPc(i)/bulkDECOMPc(i));
rhizoMIN(i) = rhizo(i,2);
bulkMIN(i) = bulk(i,2);
effectMIN(i) = log(rhizoMIN(i)/bulkMIN(i));
rhizoBioC(i) = rhizo(i,4);
bulkBioC(i) = bulk(i,4);
effectBioC(i) = log(rhizoBioC(i)/bulkBioC(i));
rhizoSOC(i) = rhizo(i,5);
bulkSOC(i) = bulk(i,5);
effectSOC(i) = log(rhizoSOC(i)/bulkSOC(i));
end

aDECOMPc = mean(effectDECOMPc(1:15));
amin = mean(effectMIN(1:15));
abioc = mean (effectBioC(1:15));
asoc = mean (effectSOC(1:15));

%E = cat(2,effectDECOMPc,effectMIN,effectCO2,effectBioC,effectEnzC);
a = cat(2,aDECOMPc, amin, abioc, asoc);
headers = ['DECeffect '; 'MINeffect '; 'BioCeffect'; 'SOCeffect '];
heads = cellstr(headers);
csvwrite_with_headers('effect2.csv',a,heads)

load 'Allison_rhizosphere_Input3_all_depths_means.mat'

DECOMPc = DECOMPc_mean_last_row_merged.';
MIN = min_mean_last_row_merged.';
BioC = BioC_mean_last_row_merged.';
SOC = SOC_mean_last_row_merged.';
C = cat(2,DECOMPc,MIN,BioC,SOC);

headers = ['DECrhizo '; 'MINrhizo '; 'BioCrhizo'; 'SOCrhizo '];
heads = cellstr(headers);

csvwrite_with_headers('rhizo3.csv',C,heads)
xlswrite('rhizo3.xls',C)

%import data
bulk = xlsread ('bulk.xls');
rhizo = xlsread ('rhizo3.xls');
vol = xlsread ('med_rhizo_volume.xlsx'); 

%DECOMPc, MIN, CO2, BioC, EnzC
effectDECOMPc=zeros(100,1);
rhizoDECOMPc=zeros(100,1);
bulkDECOMPc=zeros(100,1);
effectMIN=zeros(100,1);
rhizoMIN=zeros(100,1);
bulkMIN=zeros(100,1);
effectBioC=zeros(100,1);
rhizoBioC=zeros(100,1);
bulkBioC=zeros(100,1);
effectSOC=zeros(100,1);
rhizoSOC=zeros(100,1);
bulkSOC=zeros(100,1);

for i=1:100;
rhizoDECOMPc(i) = rhizo(i,1);
bulkDECOMPc(i) = bulk(i,1);
effectDECOMPc(i) = log(rhizoDECOMPc(i)/bulkDECOMPc(i));
rhizoMIN(i) = rhizo(i,2);
bulkMIN(i) = bulk(i,2);
effectMIN(i) = log(rhizoMIN(i)/bulkMIN(i));
rhizoBioC(i) = rhizo(i,4);
bulkBioC(i) = bulk(i,4);
effectBioC(i) = log(rhizoBioC(i)/bulkBioC(i));
rhizoSOC(i) = rhizo(i,5);
bulkSOC(i) = bulk(i,5);
effectSOC(i) = log(rhizoSOC(i)/bulkSOC(i));
end

aDECOMPc = mean(effectDECOMPc(1:15));
amin = mean(effectMIN(1:15));
abioc = mean (effectBioC(1:15));
asoc = mean (effectSOC(1:15));

%E = cat(2,effectDECOMPc,effectMIN,effectCO2,effectBioC,effectEnzC);
a = cat(2,aDECOMPc, amin, abioc, asoc);
headers = ['DECeffect '; 'MINeffect '; 'BioCeffect'; 'SOCeffect '];
heads = cellstr(headers);
csvwrite_with_headers('effect3.csv',a,heads)

toc
