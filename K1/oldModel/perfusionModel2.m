clear;close all;
%% Inputs %%
tau = 1;
K = 0.3714;
tend = 30;
cellType = 1;
shift = 0;
Perfusion = 1;
R0 = 6.2;
fileName = 'pData.mat';
%% Load variables from definitions
reversibleLogicals;
internalLogicals;
perfusionAlpha;
perfusionInputs;
perfusionInitialConditions;
load('stoichMatrix.mat')
tspan = 0:1:tend;

Rvec = [R0];
objVec = [];
exitVec = [];
ttestVec = [0];
save(fileName, 'Rvec','objVec','exitVec','ttestVec');

dCdt_particular = @(t,C) dCdt(t,C,Perfusion, cellType, shift, fileName, tau, K);

[t, C] = ode45(dCdt_particular, tspan,C0);
load(fileName)
function derivativeSystem = dCdt(t,C,Perfusion, cellType, shift, fileName, tau, K)
    reversibleLogicals;
    internalLogicals;
    perfusionAlpha;
    perfusionInputs;
    perfusionInitialConditions;
    load('stoichMatrix.mat')
    rates = instantRates(C,t,Perfusion,cellType,shift, fileName);
    derivativeSystem = (C_i_in + rates*tau-(1-alphaComponents)*(1+K).*C)/tau;
end
