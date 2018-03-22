maxTime = 24; %hr

maxTime = (maxTime - 2/3) * 3600; %give extra time


%% 
%guess

%x0 = log(x0);
load('ParticleSwarm.mat');
x0 = x;
lb = x0 - log(3);
ub = x0 + log(3);   

options = optimoptions('patternsearch','UseCompletePoll',true,...
    'UseParallel',true,'Cache','on','CacheSize',5e4,'MaxTime',maxTime,'MaxFunctionEvaluations',188000);

parpool(22)

fun = @batchFunction;
[x,fval,exitflag,output]  = patternsearch(fun,x0,[],[],[],[],lb,ub,[],options);
fileName = 'PatternSearch.mat';
save(fileName,'x','fval','exitflag','output');
