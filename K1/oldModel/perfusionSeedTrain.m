clear;close all;
%% Inputs %%
cellType = 1;
shift = 0;
Perfusion = 0;
R0 = 6.2;
a = 0.7968;
V0 = 50E-3; %L
tend = log(100/V0)/a;
fileName = 'pSeedTrainData.mat';
%% Load variables from definitions
reversibleLogicals;
internalLogicals;
perfusionInputs;
perfusionInitialConditions;
load('stoichMatrix.mat')
tspan = [0 tend];

Rvec = [R0];
objVec = [];
exitVec = [];
ttestVec = [0];
save(fileName, 'Rvec','objVec','exitVec','ttestVec');

dCdt_particular = @(t,C) dCdt(t,C,Perfusion, cellType, shift, fileName, a);

C0 = [C0; V0];
C0(5) = 10*2.31/5;
[t, C] = ode45(dCdt_particular, tspan,C0);
load(fileName)
function derivativeSystem = dCdt(t,C,Perfusion, cellType, shift, fileName, a)

    V = C(end);
    flow = a * V;
    reversibleLogicals;
    internalLogicals;
    perfusionInputs;
    perfusionInitialConditions;
    load('stoichMatrix.mat')
    rates = instantRates(C(1:end-1),t,Perfusion,cellType,shift, fileName);
    derivativeSystem = rates + C_i_in*flow/V - C(1:end - 1)*flow/V;
    derivativeSystem = [derivativeSystem; a*V];
end
