
%close all;
figure(5);clf;
Titles=[{'BIOM'},{'ANTI'},{'GLUC'},{'LAC '},{'ALA '},{'ASN '},{'ASP '},{'C-C '},{'GLN '},...
    {'GLY '},{'SER '},{'NH3 '},{'GLU '},{'VCD '},{'R  '}];

if Perfusion
    VCD = @(t)0.8;
else
    VCD = @(t)(1-1./(1+exp(-0.5*(t-9))));
end
L = [5,2,8,12,1,3,4,6,9,11,15,13,10];
 
    for j=1:length(L)
      subplot(3,5,j);
      plot(t,C(L(j),:));
      title(Titles(j))
      xlim([0,t(end)]);
      ylim([0,inf]);
    end
 
    
    
   subplot(3,5,14); 
   plot(t,(C(5,:)/2.31).*VCD(t));
   title(Titles(14))
   xlim([0,t(end)]);
   
   subplot(3,5,15); 
   plot(t,R);
   title(Titles(15))
   xlim([0,t(end)]);
   
    hold off
    
