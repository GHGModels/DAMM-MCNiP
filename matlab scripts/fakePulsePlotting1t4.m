%set up data values
leg = {'DOC 1', 'DOC 2','DOC 3','DOC 4'};
ti = {'DOC pool'};
filename{1} = 'vars1.mat';
filename{2} = 'vars2.mat';
filename{3} = 'vars3.mat';
filename{4} = 'vars4.mat';

%plot DOC timeseries after 1:4 fake water pulses
figure
for h = 1:4
    if h == 4
        load(filename{h})
        cols = EC(Nt-4561:Nt,:);
        plot((1:4562)./24,cols(:,1),'LineWidth',3)
    else
        load(filename{h})
        cols = EC(Nt-4561:Nt,:);
        plot((1:4562)./24,cols(:,1),'LineWidth',3)
        hold on;
    end
end
    legend(leg)
    title(ti)
    xlabel('Day of Year')
    ylabel('DOC (mg C cm^-^3)')
    lim=axis; 
    p=patch([2281/24 2281/24 (2281+192)/24 (2281+192)/24],[lim(3) lim(4) lim(4) lim(3)],'r');
      set(p,'FaceAlpha',0.5,'FaceColor',[0.8 0.8 0.8],'EdgeColor','none');
    xlim([2200/24 2500/24]);

%calc sums as barplots for 1:4 fake water pulses
for h = 1:4
    load(filename{h})
    cols = EC(Nt-4561:Nt,:);
    sums(:,h) = sum(cols(:,1));
end

for h = 1:4
    load(filename{h})
    cols = CMIN(Nt-4561:Nt,:);
    sums(:,h) = sum(cols(:,1));
end

sums

    
    