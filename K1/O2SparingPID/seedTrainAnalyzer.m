

figure(1);clf;plot(t_total,VCD_total);xlabel('Total Time (days)');ylabel('Viable Cell Density (million/mL)');
figure(3);clf;plot(t_total,C(2,:)*126/1000);xlabel('Total Time (days)');ylabel('MAB Density (g/mL)');
figure(2);clf;plot(t_total,VCD_total.*vesselSize_total);xlabel('Total Time (days)');ylabel('Total Cells (million)');
figure(4);clf;plot(t_total,C(2,:)*126/1000.*vesselSize_total/1e6);xlabel('Total Time (days)');ylabel('MAB (kg)');
figure(5);clf;plot(t_total,log(VCD_total.*vesselSize_total));xlabel('Total Time (days)');ylabel('MAB (kg)');

b = log(VCD_total.*vesselSize_total)';
A = t_total';
A = [0*t_total'+1,A];
beta = A\b;
r = beta(2);
tau = log(2)/r*24; %doubling time in days


if system == 1 && writeFile
    fileName1 = 'initialConditionsPerfusion.mat';
    C0 = C(:,end);
    save(fileName1,'C0');
    
    fileName2 = 'PerfusionSeedTrain.mat';
    save(fileName2,'tau','Tend','rxnContainer','cContainer', 'tContainer', 'Rcontainer', 'extentContainer');
    
  
elseif system == 0 && writeFile
    fileName1 = 'initialConditionsBatch.mat';
    C0 = C(:,end);
    save(fileName1,'C0');
    
    fileName2 = 'BatchSeedTrain.mat';
    save(fileName2,'tau','Tend','rxnContainer','cContainer', 'tContainer', 'Rcontainer', 'extentContainer');
    
end


