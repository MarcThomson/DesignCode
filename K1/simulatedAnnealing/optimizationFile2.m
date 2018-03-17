%guess
clear;close all;clc
load('ClusterTestPatOut')
%x0 = log(x0);
lb = x0 - log(3);
ub = x0 + log(3);
%x0 = log(x0);
%lb = x0 - log(3);
%ub = x0 + log(3);

%   options = optimoptions('patternsearch','Display','iter','PlotFcn',@psplotbestf,'UseCompletePoll',true,'UseParallel',true,'Cache','on','MaxIterations',100);

options = optimoptions('simulannealbnd','PlotFcns',{@saplotbestx,...
          @saplotbestf,@saplotx,@saplotf},...
          'MaxTime',3600*1.5,'Display','iter','ReannealInterval',100,...
          'TemperatureFcn',@temperaturefast...
          );

%parpool(2)
%lb = zeros(47,1);
%ub = 10000*ones(47,1);
fun = @batchFunction;
tic
%[x,fval,exitflag,output]  = patternsearch(fun,x0,[],[],[],[],lb,ub,[],options)
[x,fval,exitflag,output] = simulannealbnd(fun,x0,lb,ub,options);
toc