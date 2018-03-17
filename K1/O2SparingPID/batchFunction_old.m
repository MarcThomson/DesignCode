function [f_obj] = batchFunction(parameters) 
%% Inputs %%
tend = 10;
cellType = 1;
shift = 0;
Perfusion = 0;
h = 0.05;
shiftDay = 3;
%% Load variables from definitions

parameters = exp(parameters);



reversibleLogicals;
internalLogicals;
initialConditions;
load('stoichMatrix.mat')
t = 0:h:tend;

I = find(~internal);
%C0 = C0(1:end-1);
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
failed = 0;
for i = 1:length(t)-1
    if t(i)>shiftDay
        shift = 1;
        load('stoichMatrix31.mat');
    end
    
   [K1, R1,~,~,~] =  instantRatesV2(C(:,i),t(i),Perfusion,cellType, shift, R(i), parameters);
   K1 = K1*h; K1(14) = 0; 
   [K2, R2,~,~,~] = instantRatesV2(C(:,i)+K1/2, t(i) + h/2 ,Perfusion,cellType, shift, R1, parameters);
   K2 = K2*h; K1(14) = 0;
   [K3, R3,~,~,~] = instantRatesV2(C(:,i)+K2/2, t(i) + h/2 ,Perfusion,cellType, shift, R2, parameters);
   K3 = K3*h;K1(14) = 0;
   [K4, R4,~,~,~] = instantRatesV2(C(:,i)+K3, t(i) + h ,Perfusion,cellType, shift, R3, parameters);
   K4 = K4*h;  K1(14) = 0;   
   
   rates = (K1+2*K2+2*K3+K4)/6;
   rates(end-1) = 0;
   

   Rtest = (R1+2*R2+2*R3+R4)/6;
   C(:,i+1) = C(:,i) + rates;
   [~, R(i+1),v(:,i+1),objVec(i+1),exitVec(i+1)] =...
       instantRatesV2(C(:,i+1), t(i) + h ,Perfusion,cellType, shift, Rtest, parameters);
   
   if  t(i)-floor(t(i))==0 && t(i)>=3
       
       C(1,i+1) = C(1,i+1)*0.97 + 0.03*6.00; %ALA
       C(4,i+1) = C(4,i+1)*0.97 + 0.03*15.04; %ASP
       C(8,i+1) = C(8,i+1)*0.97 + 0.03*444.4; %GLC
       C(6,i+1) = C(6,i+1)*0.97 + 0.03*4.68; %CC 
       C(3,i+1) = C(3,i+1)*0.97 + 0.03*54.01; %ASN 
       C(11,i+1) = C(11,i+1)*0.97+0.03*6.00; %GLY 
       C(10,i+1) = C(10,i+1)*0.97+0.03*6.00; %GLU
       C(15,i+1) = C(15,i+1)*0.97+0.03*45.16; %SER 
       %ALA
   end
  
   
   Ineg = find(C(:,i+1)<0);
   C(Ineg,i+1) = 0;
  
   if abs(R(i+1))>10 || C(8,i+1)<10
       R(i+1) = 0;
       failed = 1;
       break
   end
end

if failed ==1
    C = 0*C;
    C = C - 1000;
end

load('expData.mat');
I = [1;2;3;4;5;6;8;9;10;11;12;13;15];
f_obj = 0;
for i = 1:length(I)
    C_pred = interp1(t,C(I(i),:),t_exp);
    f_obj = f_obj + sum(((C_pred-C_exp(i,:))./C_exp(i,:)).^2);
end


