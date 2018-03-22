
close all;
Titles=[{'BIOM'},{'ANTI'},{'GLUC'},{'LAC '},{'ALA '},{'ASN '},{'ASP '},{'C-C '},{'GLN '},...
    {'GLY '},{'SER '},{'NH3 '},{'GLU '},{'R  '},{'VCD '}];

L = [11,4,20,27,2,6,8,13,21,25,38,32,23,41];
 
    for j=1:length(L)
      subplot(3,5,j);
      plot(t,C(L(j),:));
      title(Titles(j))
      xlim([0,t(end)]);
      ylim([0,inf]);
    end
    
  pVCD=VCD(t);
  pVCD=pVCD.*C(11,:)/2.31;
  subplot(3,5,15);
  %plot(t,VCDtotal);
  title(Titles(15))
  xlim([0,t(end)]);
  hold off
    
