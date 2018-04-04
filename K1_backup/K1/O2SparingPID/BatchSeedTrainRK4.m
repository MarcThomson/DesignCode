clear;close all;clc;
%load('ObjFun.mat');
load('R2.mat')
%guess
parameters = x;
parameters = exp(parameters);
%parameters  = log(x0);

critDensity = 10;

system = 1; % 1 = perfusion, 0 = batch

if system == 0
    vesselSize =  [1,4,15,60,250,1000,4E3,20E3,80E3,400E3];
elseif system ==1
    vesselSize =  [1,4,15,60,250,1000,4E3,20E3,100E3];
end
%% Inputs %%
tend = 10;
cellType = 1;
shift = 0;
Perfusion = 0;
shiftDay = 100;
h = 0.005;
%% Load variables from definitions
reversibleLogicals;
internalLogicals;
initialConditionsSeedTrain;
initialConditionsCryo;
load('stoichMatrix.mat')
t = 0:h:tend;

I = find(~internal);
C0 = initialConditionsCryo_vec(I);

C_total = [C0];
t_total = [0];
R_total = [initialConditions_vec(end)];
obj_total = [-12];
exit_total = [1];
v_total = zeros(34,1);
VCD_total = [C_total(5)/2.31];
vesselSize_total = [vesselSize(1)];
Tend = [];

for vessel = 2:length(vesselSize)
    C0 =   (C_total(:,end) * vesselSize(vessel-1) + initialConditions_vec(I) *...
                    (vesselSize(vessel) - vesselSize(vessel-1))) / vesselSize(vessel);
           
    C0(5) = VCD_total(end)*vesselSize(vessel-1)/vesselSize(vessel)*2.31;            
   
    C = zeros(length(C0),length(t));
    R = zeros(1,length(t));
    objVec = zeros(1,length(t));
    exitVec = zeros(1,length(t));
    VCD = zeros(1,length(t));
    v = zeros(34,length(t));  
    vesselSize_vec = vesselSize(vessel)*ones(1,length(t));
    
    C(:,1) = C0;
    R(1) = R_total(end);
    exitVec(1) = 1;
    objVec(1) = -12;
    VCD(1) = C0(5)/2.31; 
    
    i = 1;
    while VCD(i) < 10
        [K1, R1,~,~,~] =  instantRatesV2(C(:,i),t(i),Perfusion,cellType, shift, R(i), parameters,shiftDay);
        K1 = K1*h;

        [K2, R2,~,~,~] = instantRatesV2(C(:,i)+K1/2, t(i) + h/2 ,Perfusion,cellType, shift, R1, parameters,shiftDay);
        K2 = K2*h;

        [K3, R3,~,~,~] = instantRatesV2(C(:,i)+K2/2, t(i) + h/2 ,Perfusion,cellType, shift, R2, parameters,shiftDay);
        K3 = K3*h;

        [K4, R4,~,~,~] = instantRatesV2(C(:,i)+K3, t(i) + h ,Perfusion,cellType, shift, R3, parameters, shiftDay);
        K4 = K4*h;

        rates = (K1+2*K2+2*K3+K4)/6;
        Rtest = (R1+2*R2+2*R3+R4)/6;
        C(:,i+1) = C(:,i) + rates;


        C(:,i+1) = C(:,i) + rates;
   [~, R(i+1),v(:,i+1),objVec(i+1),exitVec(i+1)] =...
       instantRatesV2(C(:,i+1), t(i) + h ,Perfusion,cellType, shift, Rtest, parameters, shiftDay);


        VCD(i+1) = C(5,i+1)/2.31 * (1-0.5259/(1+353.3*exp(-0.9381*t(i+1))));

        i = i + 1;
    end
C_total = [C_total,C(:,2:i)];
t_total = [t_total, t_total(end) + t(2:i)];
R_total = [R_total, R(2:i)];
obj_total = [obj_total,objVec(2:i)];
exit_total = [exit_total, exitVec(2:i)];
v_total = [v_total, v(:,2:i)];
VCD_total = [VCD_total,VCD(2:i)]; 
vesselSize_total = [vesselSize_total, vesselSize_vec(2:i)];
Tend = [Tend;t(i)];   


if vessel>=8 
    if system ==0
        systemString = 'Batch';
    elseif system ==1
        systemString = 'Perfusion';
    end
    
    Volume = num2str(vesselSize(vessel)/1000);
    fileName = [systemString,Volume,'.mat'];
    Dt = (4*vesselSize(vessel)/(3*pi))^(1/3)/100; %m
    A = pi/4*Dt^2;
    C_O2_vec = smooth(C(14,2:i),15);
    C_CO2_vec = smooth(C(7,2:i),15);
    t2 = t;
    t = t(2:i);
    save(fileName,'t','shiftDay','C_CO2_vec','C_O2_vec','Dt','A');
    t = t2;
end   

end



C = C_total;
t = t_total;
R = R_total;
obj = obj_total;
exit_vec = exit_total;
v = t_total;
VCD = VCD_total;

seedTrainAnalyzer
%plotToolV2