maxTime = 24; %hr

maxTime = (maxTime - 1/3) * 3600; %give extra time


%% 
guess

x0 = log(x0);
lb = x0 - log(3);
ub = x0 + log(3);   
options = optimoptions(@particleswarm,'Display','iter','SwarmSize',470,'HybridFcn',@patternsearch,'UseParallel', true);


%parpool(22)

fun = @batchFunction;
[x,fval,exitflag,output]  = particleswarm(fun,47,lb,ub,options);
fileName = 'ParticleSwarm.mat';
save(fileName,'x','fval','exitflag','output');