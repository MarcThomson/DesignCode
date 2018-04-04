clear;%close all;

%load('ObjFun.mat');
load('R2.mat')
%guess
parameters = x;
parameters = exp(parameters);
%parameters  = log(x0);

%% Inputs %%
tend = 30;
cellType = 1;
shift = 0;
Perfusion = 1;
h = 0.005;
alpha = (1+19.2703)*(1-0.99975);
tau = 1;
shiftDay = 31;
Dt = (4*500E3/(3*pi))^(1/3)/100; %m
A = pi/4*Dt^2;
%% Load variables from definitions
reversibleLogicals;
internalLogicals;
load('stoichMatrix.mat')
initialConditionsSeedTrain;

load('initialConditionsPerfusion.mat')

I = find(~internal);
C_i_in = initialConditions_vec(I);
%C0 = C0(1:end-1);
C0 = ( 400 * C_i_in + 100*C0)/(500);


t = 0:h:tend;

I = find(~internal);
C0 = [C0;C0(5)];

C = zeros(length(C0),length(t));
R = zeros(1,length(t));
objVec = zeros(length(t),1);
exitVec = zeros(length(t),1);

R(1) = 6.2;%%%% ADD LATER
C(:,1) = C0;
exitVec(1) = 1;
objVec(1) = -12;
v = zeros(34,length(t));

I_notBIOM = [1:4,6:length(C0)-1];
for i = 1:length(t)-1
    
    if t(i)>shiftDay
        shift = 1;
        load('stoichMatrix31.mat');
    end
    
    
   [K1, R1,~,~,~] =  instantRatesV2(C(1:end-1,i),t,Perfusion,cellType, shift, R(i), parameters, shiftDay);
   K1(I_notBIOM) = K1(I_notBIOM) - C(I_notBIOM,i)/tau + C_i_in(I_notBIOM)/tau;
   K1 = [K1; K1(5)];
   K1(5) = K1(5) - alpha/tau*C(5,i);
   K1 = K1*h; %K1(14) = 0;
   
   [K2, R2,~,~,~] = instantRatesV2(C(1:end-1,i)+K1(1:end-1)/2, t + h/2 ,Perfusion,cellType, shift, R1, parameters, shiftDay);
   K2(I_notBIOM) = K2(I_notBIOM) - C(I_notBIOM,i)/tau + C_i_in(I_notBIOM)/tau;
   K2 = [K2; K2(5)];
   K2(5) = K2(5) - alpha/tau*C(5,i);
   K2 = K2*h; %K2(14) = 0;
   
   
   [K3, R3,~,~,~] = instantRatesV2(C(1:end-1,i)+K2(1:end-1)/2, t + h/2 ,Perfusion,cellType, shift, R2, parameters, shiftDay);
   K3(I_notBIOM) = K3(I_notBIOM) - C(I_notBIOM,i)/tau + C_i_in(I_notBIOM)/tau;
   K3 = [K3; K3(5)];
   K3(5) = K3(5) - alpha/tau*C(5,i);
   K3 = K3*h;% K3(14) = 0;
   
   
   [K4, R4,~,~,~] = instantRatesV2(C(1:end-1,i)+K3(1:end-1), t + h ,Perfusion,cellType, shift, R3, parameters, shiftDay);
   K4(I_notBIOM) = K4(I_notBIOM) - C(I_notBIOM,i)/tau + C_i_in(I_notBIOM)/tau;
   K4= [K4; K4(5)];
   K4(5) = K4(5) - alpha/tau*C(5,i);
   K4 = K4*h; %K4(14) = 0;
   
   Rtest = (R1+2*R2+2*R3+R4)/6;
   rates = (K1+2*K2+2*K3+K4)/6;
   %rates(14) = 0;
   

   C(:,i+1) = C(:,i) + rates;
   
   
   
   [~, R(i+1),v(:,i+1),objVec(i+1),exitVec(i+1)] =...
       instantRatesV2(C(1:end-1,i+1), t + h ,Perfusion,cellType, shift, Rtest, parameters, shiftDay);
  
   
   
   
%    Ineg = find(C(:,i+1)<0);
%    C(Ineg,i+1) = 0;
  
   if abs(R(i))>10
       break
   end
end

figure(1);clf;
subplot(2,1,1);plot(t,C(5,:)/2.31); hold on;
subplot(2,1,1);plot(t,C(end,:)/2.31);
subplot(2,1,2);plot(t,C(5,:)./C(end,:))
plotToolV2
reactionScraperPerfusion

C_O2_vec = smooth(C(14,:),15);
C_CO2_vec = smooth(C(7,:),15);
fileName = 'Perfusion_Output.mat';
save(fileName,'t','shiftDay','C_CO2_vec','C_O2_vec','A','Dt');

