%Plots output of a simulation of CHO-K1
%close all;

figure();

% names of components as displayed on the graphs
Titles=[{'BIOM'},{'ANTI'},{'GLUC'},{'LAC '},{'ALA '},{'ASN '},{'ASP '},{'C-C '},{'GLN '},...
    {'GLY '},{'SER '},{'NH3 '},{'GLU '},{'VCD '},{'R  '}];

% load experimental data
load('expData.mat')

% if batch reactor, scale the BIOM to get the VCD
if Perfusion
    VCD = @(t)1;
else
    VCD = @(t)1-0.5259/(1+353.3*exp(-0.9381*t));
end

% indices for components
L = [5,2,8,12,1,3,4,6,9,11,15,13,10];
L_exp = [5,2,7,11,1,3,4,6,8,10,13,12,9];

for j=1:length(L)
      subplot(3,5,j);
      plot(t,C(L(j),:),'r-');
      hold on;
      % don't plot experimental data in the perfusion case
      if ~Perfusion
        plot(t_exp,C_exp(L_exp(j),:),'bs');
      end
      hold off
      title(Titles(j))
      xlim([0,t(end)]);
      ylim([0,inf]);
end
 
    
%plot VCD
subplot(3,5,14); 
plot(t,(C(5,:)/2.31).*VCD(t),'r-');
hold on;
if ~Perfusion
    plot(t_exp,VCD_exp,'bs');
end
hold off
title(Titles(14))
xlim([0,t(end)]);

%Plot R
subplot(3,5,15); 
plot(t,R,'r-');
hold on
if ~Perfusion
    plot(t_exp,R_exp,'bs');
end
hold off
title(Titles(15))
xlim([0,t(end)]);
