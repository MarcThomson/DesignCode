clear;close all;
%% Inputs %%
K = 0.10;
tau = 1; %Defined twice!
tend = 40;
%% Load variables from definitions
reversibleLogicals;
internalLogicals;
perfusionAlpha;
perfusionInputs;
perfusionInitialConditions;
load('stoichMatrix.mat')
tspan = [0 tend];

[t, C] = ode45(@dCdt,tspan,C0);

function derivativeSystem = dCdt(t,C)
    tau = 1;
    K = 0.10;
    
    reversibleLogicals;
    internalLogicals;
    perfusionAlpha;
    perfusionInputs;
    perfusionInitialConditions;
    load('stoichMatrix.mat')
    derivativeSystem = (C_i_in + perfusionRates(C)*tau-(1-alphaComponents)*(1+K).*C)/tau;
end
