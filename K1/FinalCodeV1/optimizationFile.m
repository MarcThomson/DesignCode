% configure the time on the computer
maxTime = 24; %hr
maxTime = (maxTime - 1/3) * 3600; %give extra time for saving

% load a guess
load('NolanLeeParameters.mat');
x0 = parameters;
x0 = log(x0);

% give appropriate bounds
lb = x0 - log(3); 
ub = x0 + log(3);   

% use multiple cores on the supercomputer
parpool(22)

% call the modified batch script
fun = @(x)batchFunction(x,0);

% particle swarm options
options = optimoptions(@particleswarm,'Display','iter','SwarmSize',470,...
    'HybridFcn',@patternsearch,'UseParallel', true,'MaxTime',maxTime);
[x,fval,exitflag,output]  = particleswarm(fun,47,lb,ub,options);

% pattern search  
% options = optimoptions('patternsearch','UseCompletePoll',true,...
%    'UseParallel',true,'Cache','on','CacheSize',5e4,'MaxTime',maxTime);
% [x,fval,exitflag,output]  = patternsearch(fun,x0,[],[],[],[],lb,ub,[],options);

% save the parameters
fileName = 'ParametersOut.mat';
save(fileName,'x','fval','exitflag','output');