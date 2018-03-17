maxTime = 24; %hr

maxTime = (maxTime - 2/3) * 3600; %give extra time


%% 
guess

x0 = log(x0);
lb = x0 - log(3);
ub = x0 + log(3);   
options = optimoptions(@particleswarm,'Display','iter','SwarmSize',470,'UseParallel', true,'MaxTime',maxTime);


parpool(22)

fun = @batchFunction;
[x,fval,exitflag,output]  = particleswarm(fun,47,lb,ub,options);
fileName = 'ParticleSwarm.mat';
save(fileName,'x','fval','exitflag','output');
