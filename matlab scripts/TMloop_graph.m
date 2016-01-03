x = [0.1 0.2 0.3 0.4 0.5 0.6];
y = [5 10 15 20 25];

for j=1:6
    for i= 1:5
filename = ['temp_',int2str(y(i)),'_soilM_',num2str(x(j)),'.mat'];
load(filename)

zc(i,j) = resp;
zn(i,j) = nmin;

    end
end

figure
surf(x,y,zc)
title('C mineralization')
xlabel('Soil Moisture (cm3 H20/cm3 soil)')
ylabel('Temperature (C)')

figure
surf(x,y,zn)
title('N mineralization')
xlabel('Soil Moisture (cm3 H20/cm3 soil)')
ylabel('Temperature (C)')
