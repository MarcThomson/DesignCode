x1 = linspace(0,10);
P = test2Fxn(x1);

fileName = 'test2.mat';

save(fileName, 'P','x1');