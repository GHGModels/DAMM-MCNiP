%loop to set up scaling coefficients
colsfull{1} = [SOC MIC_C EC];
colsfull{2} = [overflow_C DECOM_C UPT_C CMIN Growth_C soilM];
colsfull{3} = [DOC DON];
colsfull{4} = [DECOM_N UPT_N NMIN Growth_N soilM];
colsfull{5} = [O2 DOC(1:Nt,:) DON(1:Nt,:) soilM];
leg{1} = {'SOC', 'MIC C','ENZ C'};
leg{2} = {'overflow C','decomposition C','uptake C','C min','Growth C','soilM'};
leg{3} = {'DOC', 'DON'};
leg{4} = {'decomposition N','uptake N','N min','Growth N','soilM'};
leg{5} = {'O2','DOC','DON','soilM'};
ti{1} = {'C pools'};
ti{2} = {'C fluxes'};
ti{3} = {'CN pools'};
ti{4} = {'N fluxes'};
ti{5} = {'aqueous vars'};
for h = 1:5
    cols = colsfull{h}(Nt-4561:Nt,:);
    globmean = mean(SOC);
    newcol = zeros(size(cols));
        for i = 1:size(cols,2)
        tempmean = mean(cols(:,i));
        newcol(:,i) = cols(:,i).*globmean./tempmean; %fix this shit right here
        end
        
        %if h == 2
        %newcol(:,1) = newcol(:,1)*0.1;
        %end

    figure
    for i = 1:size(newcol,2)
        if i ~= size(newcol,2)
    plot(newcol(:,i),'LineWidth',3)
    hold on;
        else
    plot(newcol(:,i),'LineWidth',3)
        end
    end
    legend(leg{h})
    title(ti{h})
    xlabel('timesteps')
    ylabel('normalized to SOC range, mg C/cm^3')
    lim=axis; 
    p=patch([2281 2281 2281+51 2281+51],[lim(3) lim(4) lim(4) lim(3)],'r');
      set(p,'FaceAlpha',0.5,'FaceColor',[0.8 0.8 0.8],'EdgeColor','none');
    xlim([2270 2400]);
end
