%Plots output of a simulation of CHO-K1
%close all;

% initialize figure
figure();set(gca,'FontSize',20);

% names of components as displayed on the graphs
Titles=[{'BIOM'},{'ANTI'},{'GLUC'},{'LAC '},{'ALA '},{'ASN '},{'ASP '},...
    {'C-C '},{'GLN '},{'GLY '},{'SER '},{'NH3 '},{'GLU '},{'VCD '},{'R  '}]

% load experimental data
load('expData.mat')

% if batch reactor, scale the BIOM to get the VCD
if Perfusion
    VCD = @(t)1;
else
    VCD = @(t)1-0.5259./(1+353.3*exp(-0.9381*t)); %emperical function
end

% indices for components
L = [5,2,8,12,1,3,4,6,9,11,15,13,10];     % in our framework
L_exp = [5,2,7,11,1,3,4,6,8,10,13,12,9];  % in the experimental framework

for j=1:length(L)
      subplot(3,5,j);
      % only plot experimental data if specified
      if plotExp
        plot(t,C(L(j),:),'r-','LineWidth',1.5);
        hold on;
        plot(t_exp,C_exp(L_exp(j),:),'b.','MarkerSize',18);
      else
          plot(t,C(L(j),:),'LineWidth',1.5);
      end
      hold off
      title(Titles(j))
      xlim([0,t(end)]);
      ylim([0,inf]);
      xlabel('Time (days)');
      ylabel('Concentration (mM)');
end
 
    
%plot VCD
subplot(3,5,14); 
if plotExp
    plot(t,(C(5,:)/2.31).*VCD(t),'r-','LineWidth',1.5);
    hold on;
    plot(t_exp,VCD_exp,'b.','MarkerSize',18);
else
    plot(t,(C(5,:)/2.31).*VCD(t),'LineWidth',1.5);
end
hold off
title(Titles(14))
xlim([0,t(end)]);
xlabel('Time (days)');
ylabel('10^6 cells/mL');

%Plot R
subplot(3,5,15); 

if plotExp
    plot(t,R,'r-','LineWidth',1.5);
    hold on
    plot(t_exp,R_exp,'b.','MarkerSize',18);
else
    plot(t,R,'LineWidth',1.5);
end
hold off
title(Titles(15))
xlim([0,t(end)]);
xlabel('Time (days)');
ylabel('R (unitless)');

