% run the batch function to compare the outputs for different parameters

% plot the outputs
plotOn = 1;

% choose the parameter set
%load('finalParameters_v1.mat')
load('NolanLeeParameters.mat')

% find the objective function
f_obj = batchFunction(log(parameters'),plotOn)





