maxTime = 24; %hr

maxTime = (maxTime - 1/3) * 3600; %give extra time


%% 
x0 =     [2.0000
    2.0000
    0.3000
  100.0000
   75.0000
    0.2000
   0
    0.0500
         0];

lb = 0*x0;
ub = lb+200; 
options = optimoptions(@particleswarm,'Display','iter','SwarmSize',100,'UseParallel', true,'MaxTime',maxTime);



parpool(22)

fun = @ PIDSparging;
[x,fval,exitflag,output]  = particleswarm(fun,9,lb,ub,options);
fileName = 'ParticleSwarm_PID.mat';
save(fileName,'x','fval','exitflag','output');