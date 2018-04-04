
%close all;
figure();
Titles=[{'BIOM'},{'ANTI'},{'GLUC'},{'LAC '},{'ALA '},{'ASN '},{'ASP '},{'C-C '},{'GLN '},...
    {'GLY '},{'SER '},{'NH3 '},{'GLU '},{'VCD '},{'R  '}];

load('expData.mat')
if Perfusion
    VCD = @(t)0.8;
else
    VCD = @(t)(1-1./(1+exp(-0.5*(t-9))));
end
L = [5,2,8,12,1,3,4,6,9,11,15,13,10];
L_exp = [5,2,7,11,1,3,4,6,8,10,13,12,9];
 
    for j=1:length(L)
      subplot(3,5,j);
      plot(t,C(L(j),:),'r-');
      hold on;
      plot(t_exp,C_exp(L_exp(j),:),'bs');
      hold off
      title(Titles(j))
      xlim([0,t(end)]);
      ylim([0,inf]);
    end
 
    
    
   subplot(3,5,14); 
   plot(t,(C(5,:)/2.31).*VCD(t),'r-');
   hold on;
   plot(t_exp,VCD_exp,'bs');
   hold off
   title(Titles(14))
   xlim([0,t(end)]);
   
   subplot(3,5,15); 
   plot(t,R,'r-');
   hold on
   plot(t_exp,R_exp,'bs');
   hold off
   title(Titles(15))
   xlim([0,t(end)]);
   
    hold off
    
