
amp = 0.0005;
refin = 0.0005;

A1=amp/2;         %seasonal amplitude
A2=0;                %daily amplitude
w1=2.*pi/Nt;
w2=2.*pi;
ref=refin/2;
t=1:Nt;
DOC_input = ref+A1.*sin(w1.*t-pi/2)+A2.*sin(w2.*t);%0.0005;%external C input into DOC pool(e.g. root exudates,root turnover) mg C cm-3 soil hour-1

figure
plot(DOC_input)
trapz(DOC_input)