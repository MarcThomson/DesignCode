% plot the seed train, do some basic analysis, and save the files

% Plot viable cell density
figure();clf;plot(t_total,VCD_total,'LineWidth',2);
xlabel('Total Time (days)');ylabel('Viable Cell Density (million/mL)');
set(gca,'FontSize',20);

% plot density of MAB
figure();clf;plot(t_total,C(2,:)*126/1000,'LineWidth',2);
xlabel('Total Time (days)');ylabel('MAB Density (g/mL)');
set(gca,'FontSize',20);

% plot total cell count
figure();clf;plot(t_total,VCD_total.*vesselSize_total,'LineWidth',2);
xlabel('Total Time (days)');ylabel('Total Cells (million)');
set(gca,'FontSize',20);

% plot total MAB
figure();clf;plot(t_total,C(2,:)*126/1000.*vesselSize_total/1e6,'LineWidth',2);
xlabel('Total Time (days)');ylabel('MAB (kg)');
set(gca,'FontSize',20);

% Plot log(total viable cells) vs time (should be linear)
figure();clf;plot(t_total,log(VCD_total.*vesselSize_total),'LineWidth',2);
xlabel('Total Time (days)');ylabel('ln(Viable Cells)');
set(gca,'FontSize',20);


% calculate the doubling time
b = log(VCD_total.*vesselSize_total)';
A = t_total';
A = [0*t_total'+1,A];
beta = A\b;
r = beta(2);
tau = log(2)/r*24; %doubling time in days

% save the outputs as initial conditions for the reactors
% save the containers with all data for future reference
if system == 1 && writeFile
    fileName1 = 'initialConditionsPerfusion.mat';
    C0 = C(:,end);
    save(fileName1,'C0');
    
    fileName2 = 'PerfusionSeedTrain.mat';
    save(fileName2,'tau','Tend','rxnContainer','cContainer',...
        'tContainer', 'Rcontainer', 'extentContainer');
    
  
elseif system == 0 && writeFile
    fileName1 = 'initialConditionsBatch.mat';
    C0 = C(:,end);
    save(fileName1,'C0');
    
    fileName2 = 'BatchSeedTrain.mat';
    save(fileName2,'tau','Tend','rxnContainer','cContainer',...
        'tContainer', 'Rcontainer', 'extentContainer');
    
end


