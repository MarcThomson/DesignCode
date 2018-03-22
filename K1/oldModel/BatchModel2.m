clear;close all;
%% Inputs %%
tend = 10;
cellType = 1;
shift = 0;
Perfusion = 0;
R0 = 6.2;
fileName = 'bData.mat';
%% Load variables from definitions
reversibleLogicals;
internalLogicals;
initialConditions;
load('stoichMatrix.mat')
tspan = 0:1:tend;

Rvec = [R0];
objVec = [];
exitVec = [];
ttestVec = [0];
save(fileName, 'Rvec','objVec','exitVec','ttestVec');

dCdt_particular = @(t,C) dCdt(t,C,Perfusion, cellType, shift, fileName);

I = find(~internal);
C0 = initialConditions_vec(I);
[t, C] = ode45(dCdt_particular, tspan,C0);
load(fileName)
function derivativeSystem = dCdt(t,C,Perfusion, cellType, shift, fileName)
    reversibleLogicals;
    internalLogicals;
    initialConditions;
    load('stoichMatrix.mat')
    rates = instantRatesV2(C,t,Perfusion,cellType,shift, fileName);
    rates(end-1) = 0;
    derivativeSystem = rates;
end
