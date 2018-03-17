maxTime = 1; %hr

maxTime = (maxTime - 1/3) * 3600; %give extra time


%% 
guess

x0 = log(x0);
lb = x0 - log(3);
ub = x0 + log(3);   

options = optimoptions('patternsearch','UseCompletePoll',true,...
    'UseParallel',true,'Cache','on','CacheSize',5e4,'MaxTime',maxTime,'MaxFunctionEvaluations',188000);

%parpool(2)

fun = @batchFunction;
[x,fval,exitflag,output]  = patternsearch(fun,x0,[],[],[],[],lb,ub,[],options);
fileName = 'PatternSearchTest.mat';
save(fileName,'x','fval','exitflag','output');