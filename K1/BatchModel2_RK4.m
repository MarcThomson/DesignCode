clear;close all;
%% Inputs %%
tend = 10;
cellType = 1;
shift = 0;
Perfusion = 0;
h = 0.05;
%% Load variables from definitions
reversibleLogicals;
internalLogicals;
initialConditions;
load('stoichMatrix.mat')
t = 0:h:tend;

I = find(~internal);
C0 = initialConditions_vec(I);

C = zeros(length(C0),length(t));
R = zeros(1,length(t));
objVec = zeros(length(t),1);
exitVec = zeros(length(t),1);

R(1) = initialConditions_vec(end);
C(:,1) = C0;
exitVec(1) = 1;
objVec(1) = -12;
v = zeros(34,length(t));
for i = 1:length(t)-1

    
   [K1, R1,~,~,~] =  instantRatesV2(C(:,i),t,Perfusion,cellType, shift, R(i));
   K1 = K1*h; 
   [K2, R2,~,~,~] = instantRatesV2(C(:,i)+K1/2, t + h/2 ,Perfusion,cellType, shift, R1);
   K2 = K2*h; 
   [K3, R3,~,~,~] = instantRatesV2(C(:,i)+K2/2, t + h/2 ,Perfusion,cellType, shift, R2);
   K3 = K3*h;    
   [K4, R4,~,~,~] = instantRatesV2(C(:,i)+K3, t + h ,Perfusion,cellType, shift, R3);
   K4 = K4*h;        
   
   rates = (K1+2*K2+2*K3+K4)/6;
   rates(end-1) = 0;
   

   C(:,i+1) = C(:,i) + rates;
   [~, R(i+1),v(:,i+1),objVec(i+1),exitVec(i+1)] =...
       instantRatesV2(C(:,i+1), t + h ,Perfusion,cellType, shift, R1);
   
   %Gluc
   if  C(8,i+1)<40 && t(i)>=3 &&t(i)-floor(t(i))==0
       C(8,i+1) = C(8,i+1) +15;
   end
   
   %CC
   if  C(6,i+1)<1 && t(i)>=3 &&t(i)-floor(t(i))==0
       C(6,i+1) = C(6,i+1) + 0.1;
   end
   
   %ASN
   if  C(3,i+1)<15 && t(i)>=3 &&t(i)-floor(t(i))==0
       C(3,i+1) = C(3,i+1) + 2.5;
   end
   
   %GLN
   if  C(9,i+1)<1.5 && t(i)>=3 &&t(i)-floor(t(i))==0
       C(9,i+1) = C(9,i+1) +0.5;
   end
   
   %SER
   if  C(15,i+1)<8 && t(i)>=3 &&t(i)-floor(t(i))==0
       C(15,i+1) = C(15,i+1) + 1.4;
   end
   
   Ineg = find(C(:,i+1)<0);
   C(Ineg,i+1) = 0;
  
end


plotToolV2


