clear;close all;clc;                   % Clear out current variable list

% Parameters guess
load('finalParameters_v1.mat')
%% Inputs %%
tend = 10;
cellType = 1;
shift = 0;
Perfusion = 0;
h = 0.005;
shiftDay = 3;
Dt = 0.947; %m
A = Dt^2/4*pi;
writeFile = 1;
%% Load variables from definitions
reversibleLogicals;
internalLogicals;
initialConditionsSeedTrain;
load('initialConditionsBatch.mat')
load('stoichMatrix.mat')
t = 0:h:tend;

I = find(~internal);

%C0 = C0(1:end-1);
C0 = ((A*3*Dt * 0.79-0.4)*initialConditions_vec(I) + 0.4*C0)/(A*3*Dt * 0.79);

C = zeros(length(C0),length(t));
R = zeros(1,length(t));
%Nrad = zeros(1,length(t));
objVec = zeros(length(t),1);
exitVec = zeros(length(t),1);
%Qtotal=zeros(length(t),1);
%QO2=zeros(length(t),1);
%QCO2=zeros(length(t),1);
%O2SatVec = QO2;
%CO2SatVec = QO2;
Volume = [];
v = zeros(34,length(t));
%dBub = QO2;


R(1) = 6.2;
C(:,1) = C0;
exitVec(1) = 1;
objVec(1) = -12;
failed = 0;
%Nrad(1) = 1*2;
%Qtotal(1)=10e-3*A;
%QCO2(1)=6E-2*Qtotal(1);
%QO2(1)=0.4*Qtotal(1);
Volume(1) = A*3*Dt * 0.79;
totalVolume =  A*3*Dt;

totalFeed = zeros(15,7); %mmol

tChanged = -100;
for i = 1:length(t)-1
    
%    yO2 = QO2(i)/Qtotal(i);
%    yCO2 = QCO2(i)/Qtotal(i);
    
%    CO2Frac = (yCO2-400E-6)/(1-400E-6);
%    airFrac = (1-yO2-CO2Frac)/(1-0.2095);
    
    if t(i)>shiftDay
        shift = 1;
        load('stoichMatrix31.mat');
    end
   
%       [FO2,FCO2,Nrad(i+1), O2Prop, CO2Prop, dBub(i)] = spargerV2(C(14,i),C(7,i),airFrac,CO2Frac,Qtotal(i),Dt,Nrad(i),t(i),shift);
%   [Qtotal(i+1),QO2(i+1),QCO2(i+1), boolChanged] = ...
%       spargerControl(Qtotal(i), O2Prop, CO2Prop, QO2(i), QCO2(i), A, t(i), tChanged);
   
%   if boolChanged
%       tChanged = t(i);
%    end
%   O2SatVec(i) = O2Prop;
%   CO2SatVec(i) = CO2Prop;
  
   
    [K1, R1,~,~,~] =  instantRatesV2(C(:,i),t(i),Perfusion,cellType, shift, R(i), parameters,shiftDay);
%    [FO2,FCO2,~,~,~,~] = spargerV2(C(14,i),C(7,i),airFrac,CO2Frac,Qtotal(i),Dt,Nrad(i),t(i),shift);
%    FO2 = FO2*3600*24;
%    FCO2 = FCO2*3600*24;
%    K1(14) = K1(14) + FO2;
%    K1(7) = K1(7) + FCO2;
    K1 = K1*h;

    [K2, R2,~,~,~] = instantRatesV2(C(:,i)+K1/2, t(i) + h/2 ,Perfusion,cellType, shift, R1, parameters,shiftDay);
%     [FO2,FCO2,~,~,~,~] = spargerV2(C(14,i)+K1(14)/2,C(7,i)+K1(7)/2,airFrac,CO2Frac,Qtotal(i),Dt,Nrad(i),t(i),shift);
%     FO2 = FO2*3600*24;
%     FCO2 = FCO2*3600*24;
%     K2(14) = K2(14) + FO2;
%     K2(7) = K2(7) + FCO2;
     K2 = K2*h;
   
    [K3, R3,~,~,~] = instantRatesV2(C(:,i)+K2/2, t(i) + h/2 ,Perfusion,cellType, shift, R2, parameters,shiftDay);
%     [FO2,FCO2,~,~,~,~] = spargerV2(C(14,i)+K2(14)/2,C(7,i)+K2(7)/2,airFrac,CO2Frac,Qtotal(i),Dt,Nrad(i),t(i),shift);
%     FO2 = FO2*3600*24;
%     FCO2 = FCO2*3600*24;
%     K3(14) = K3(14) + FO2;
%     K3(7) = K3(7) + FCO2;
    K3 = K3*h;
   
   [K4, R4,~,~,~] = instantRatesV2(C(:,i)+K3, t(i) + h ,Perfusion,cellType, shift, R3, parameters, shiftDay);
%    [FO2,FCO2,~,~,~,~] = spargerV2(C(14,i)+K3(14),C(7,i)+K3(7),airFrac,CO2Frac,Qtotal(i),Dt,Nrad(i),t(i),shift);
%    FO2 = FO2*3600*24;
%    FCO2 = FCO2*3600*24;
%    K4(14) = K4(14) + FO2;
%    K4(7) = K4(7) + FCO2;
   K4 = K4*h;
   
   
   
   rates = (K1+2*K2+2*K3+K4)/6;
   Rtest = (R1+2*R2+2*R3+R4)/6;
   C(:,i+1) = C(:,i) + rates;
   
   
   [~, R(i+1),v(:,i+1),objVec(i+1),exitVec(i+1)] =...
       instantRatesV2(C(:,i+1), t(i) + h ,Perfusion,cellType, shift, Rtest, parameters, shiftDay);
   
   if  t(i)-floor(t(i))==0 && t(i)>=3 && t(i)<10
       Volume = [Volume;Volume(end) + 0.03*A*3*Dt];
       I_not_gas = [1:6,8:13,15];
       C(I_not_gas,i+1) = C(I_not_gas,i+1)*Volume(end-1)/Volume(end);
       Feed = zeros(15,1);
       Feed(1) = 6.00;%ALA
       Feed(3) = 54.01; %ASN
       Feed(4) = 15.04; %ASP
       Feed(6) = 4.68; %CC
       Feed(8) = 444.4; %GLc
       Feed(9) = 6; %GLN
       Feed(10) = 6;% GLU
       Feed(11) = 6; %GLY
       Feed(15) = 45.16; %SER
       
       C(:,i+1) = C(:,i+1) + 0.03*totalVolume/Volume(end)*Feed;
       totalFeed(:,floor(t(i))-2) = 0.03*totalVolume*Feed;
      
   end
   
  
%     Ineg = find(C(:,i+1)<0);
%     C(Ineg,i+1) = 0;
  
   if abs(R(i+1))>10 %|| C(8,i+1)<10
       R(i+1) = 0;
       failed = 1;
       break
   end
end

if failed ==1
    C = 0*C;
    C = C - 1000;
end

% figure(1);clf;subplot(1,2,1);plot(t(1:end-1),O2SatVec(1:end-1)*100,'LineWidth',2);xlabel('time (days)');ylabel('Oxygen Concentration (% sat.)');set(gca,'FontSize',20);
% subplot(1,2,2);plot(t(1:end-1),CO2SatVec(1:end-1)*100,'LineWidth',2);xlabel('time (days)');ylabel('CO2 Concentration (% sat.)');set(gca,'FontSize',20);
% subplot(1,2,1);
% axis([0,10,0,60])
% subplot(1,2,2);
% axis([0,10,0,10])


load('expData.mat');
plotToolV2;

I = [1;2;3;4;5;6;8;9;10;11;12;13;15];
f_obj = 0;
for i = 1:length(I)
    C_pred = interp1(t,C(I(i),:),t_exp);
    f_obj = f_obj + sum(((C_pred-C_exp(i,:))./C_exp(i,:)).^2);
end

reactionScraperBatch

if writeFile                                
    C_O2_vec = smooth(C(14,:),15);
    C_CO2_vec = smooth(C(7,:),15);
    fileName = 'Batch_Output.mat';
    save(fileName,'t','shiftDay','C_CO2_vec','C_O2_vec','A','Dt','C','MAB_Produced','rxn','extent');
end
