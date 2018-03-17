clear;close all;
critDensity = 10;

vesselSize =  [1,40,180,850,4E3,20E3,80E3,400E3];%,2000E3];

%% Inputs %%
tend = 10;
cellType = 1;
shift = 0;
Perfusion = 0;
h = 0.05;
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
       [K1, R1,~,~,~] =  instantRatesV2(C(:,i),t,Perfusion,cellType, shift, R(i));
       K1 = K1*h; K1(14) = 0;
       [K2, R2,~,~,~] = instantRatesV2(C(:,i)+K1/2, t + h/2 ,Perfusion,cellType, shift, R(i));
       K2 = K2*h; K2(14) = 0;
       [K3, R3,~,~,~] = instantRatesV2(C(:,i)+K2/2, t + h/2 ,Perfusion,cellType, shift, R(i));
       K3 = K3*h; K3(14) = 0;
       [K4, R4,~,~,~] = instantRatesV2(C(:,i)+K3, t + h , Perfusion,cellType, shift, R(i));
       K4 = K4*h; K4(14) = 0;   

       rates = (K1+2*K2+2*K3+K4)/6;
       rates(end-1) = 0;


       C(:,i+1) = C(:,i) + rates;
       [~, R(i+1),v(:,i+1),objVec(i+1),exitVec(i+1)] =...
           instantRatesV2(C(:,i+1), t + h ,Perfusion,cellType, shift, R1);

       Ineg = find(C(:,i+1)<0);
       C(Ineg,i+1) = 0;
       
       VCD(i+1) = C(5,i+1)/2.31 * (1-1./(1+exp(-0.5*(t(i+1)-9))));
       
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
    
    
end



C = C_total;
t = t_total;
R = R_total;
obj = obj_total;
exit_vec = exit_total;
v = t_total;
VCD = VCD_total;


plotToolV2