
par = dmcPars(); %load parameter names and values
s2c = struct2cell(par); %convert par to a cell array
sensPar = cell2mat(s2c(6:39)); %convert selection of par to mat array
otherPar = s2c(1:5); %save rest of par

allNames = fieldnames(par); %save parameter names

hiPar = zeros(1,1); %initialize vars in loop
loPar = zeros(1,1);
sensi = zeros(1,1);

for i = 1:length(sensPar) %for each parameter value

   sensPar(i) = sensPar(i)*2; %double parameter value
   hiPar = sensPar(i);
   m2c = mat2cell(sensPar, (rep(1,34))); %convert mat back to cell
   allPar = vertcat(otherPar, m2c); %concatenate sensPars and otherPars
   par = cell2struct(allPar, allNames); %convert cell back to struct
   hiOut = dmc(par); %pass parameters to dmc
   
   sensPar(i) = sensPar(i)*0.5; %halve parameter value
   loPar = sensPar(i);
   m2c = mat2cell(sensPar, (rep(1,34))); %convert mat back to cell
   allPar = vertcat(otherPar, m2c); %concatenate sensPars and otherPars
   par = cell2struct(allPar, allNames); %convert cell back to struct
   loOut = dmc(par); %pass parameters to dmc
   
   for j = 1:length(hiOut)
   sensi(i,j) = abs(log10(abs(hiOut(j)))-log10(abs(loOut(j))))/abs(log10(abs(hiPar))-log10(abs(loPar)));
   end
   
end

