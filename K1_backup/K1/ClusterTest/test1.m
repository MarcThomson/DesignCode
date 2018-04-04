y= @(x)x.^2;

x1 = linspace(0,10);
P = y(x1);

fileName = 'test1.mat';

save(fileName, 'P','x1');