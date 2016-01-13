figure
plot(MIC_C,'LineWidth',3)
legend('seasonal model')
title('Microbial Biomass')
xlabel('timesteps')
ylabel('mg C/cm^3 soil') 

figure
plot(EC,'LineWidth',3)
legend('seasonal model')
title('Enzyme pool')
xlabel('timesteps')
ylabel('mg C/cm^3 soil')

figure
plot(SOC,'LineWidth',3)
legend('seasonal model')
title('SOC')
xlabel('timesteps')
ylabel('mg C/cm^3 soil')

figure
plot(SON,'LineWidth',3)
legend('seasonal model')
title('SON')
xlabel('timesteps')
ylabel('mg N/cm^3 soil')

figure
plot(DOC,'LineWidth',3)
legend('seasonal model')
title('DOC')
xlabel('timesteps')
ylabel('mg C/cm^3 soil')

figure
plot(DON,'LineWidth',3)
legend('seasonal model')
title('DON')
xlabel('timesteps')
ylabel('mg N/cm^3 soil')

figure
plot(CMIN,'LineWidth',3)
hold on;
plot(Nt-6072:Nt-1512,Fluxes(:,9)./(10000*10));
legend('seasonal model')
title('Soil respiration')
xlabel('timesteps')
ylabel('mg C/cm^3 soil/timestep')
xlim([Nt-6072 Nt-1512])

figure
plot(NMIN,'LineWidth',3)
legend('seasonal model')
title('N mineralization')
xlabel('timesteps')
ylabel('mg N /cm^3 soil/timestep')

figure
plot(DECOM_C,'LineWidth',3)
legend('seasonal model')
title('C decomposition')
xlabel('timesteps')
ylabel('mg C /cm^3 soil/timestep')

figure
plot(DECOM_N,'LineWidth',3)
legend('seasonal model')
title('N decomposition')
xlabel('timesteps')
ylabel('mg N /cm^3 soil/timestep')

figure
plot(UPT_C,'LineWidth',3)
legend('seasonal model')
title('C uptake')
xlabel('timesteps')
ylabel('mg C /cm^3 soil/timestep')

figure
plot(UPT_N,'LineWidth',3)
legend('seasonal model')
title('N uptake')
xlabel('timesteps')
ylabel('mg N /cm^3 soil/timestep')

figure
plot(DEATH_C,'LineWidth',3)
legend('seasonal model')
title('Microbial turnover C')
xlabel('timesteps')
ylabel('mg C /cm^3 soil/timestep')

figure
plot(DEATH_N,'LineWidth',3)
legend('seasonal model')
title('Microbial turnover N')
xlabel('timesteps')
ylabel('mg N /cm^3 soil/timestep')

figure
plot(Growth_C,'LineWidth',3)
legend('seasonal model')
title('C growth')
xlabel('timesteps')
ylabel('mg C /cm^3 soil/timestep')

figure
plot(overflow_C,'LineWidth',3)
legend('seasonal model')
title('overflow C')
xlabel('timesteps')
ylabel('mg C /cm^3 soil/timestep') 

figure
plot(soilM,'LineWidth',3)
legend('seasonal model')
title('SoilM')
xlabel('timesteps')
ylabel('mg C /cm^3 soil/timestep')

figure
plot(T,'LineWidth',3)
legend('seasonal model')
title('temperature')
xlabel('timesteps')
ylabel('mg C /cm^3 soil/timestep')

figure
plot(sol_SOC,'LineWidth',3)
legend('seasonal model')
title('sol_SOC')
xlabel('timesteps')
ylabel('mg C /cm^3 soil/timestep') 

figure
plot(O2,'LineWidth',3)
legend('seasonal model')
title('O2')
xlabel('timesteps')
ylabel('mg C /cm^3 soil/timestep')
