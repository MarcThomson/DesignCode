%guess
load('paramOutput_v3')
% lb = x0/3;
% ub = x0*3;

%x0 = log(x0);
lb = x0 - log(3);
ub = x0 + log(3);   

options = optimoptions('patternsearch','Display','iter',...
    'PlotFcn',@psplotbestf,'UseCompletePoll',true,...
    'UseParallel',true,'Cache','on','MaxTime',3000);

parpool(2)
%lb = zeros(47,1);
%ub = 10000*ones(47,1);
fun = @batchFunction;
tic
[x,fval,exitflag,output]  = patternsearch(fun,x0,[],[],[],[],lb,ub,[],options)
toc