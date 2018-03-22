function componentRates = instantRates(C,t,Perfusion,cellType,shift, fileName)
%% Extra inputs
Rtol = 0.01;

if Perfusion
    VCD = 0.8;
else
    VCD = (1-1./(1+exp(-0.5*(t-9))));
end
%shift = 0;
%% Call other parameters
parameterDefinitions;
reversibleLogicals;
internalLogicals;
load('stoichMatrix.mat')
load(fileName);

numRxns = length(stoichMatrix(1,:));
numComponents = length(stoichMatrix(:,1));


%%Modify concentrations
I_external = find(internal == 0);
C_new = zeros(numComponents,1);
C_new(I_external) = C;

%% Optimization Initialization
tol = 1e-10;
options = optimoptions('quadprog','OptimalityTolerance',tol,...
    'MaxIterations',30000,'Algorithm','interior-point-convex',...
    'StepTolerance',1e-12,'Display','off');

Aeq = stoichMatrix;
for j=flip(1:length(I_external))
    Aeq=[Aeq(1:I_external(j)-1,:);Aeq(I_external(j)+1:end,:)];
end
beq = zeros(length(Aeq(:,1)),1);

A = zeros(numRxns,numRxns);
for j=flip(1:numRxns)
    if reversible(j)==1
        A=[A(1:j-1,:);A(j+1:end,:)];
    else
        A(j,j) = -1;
    end
end
b = zeros(length(A(:,1)),1);


[~,IR] = min(abs(ttestVec-t));
Rold = Rvec(IR);
R = Rold;
y = [C_new; shift;  R; cellType];
inputs = mat2cell(y,ones(1,length(y)),1);
x0 = rateLaws(inputs{:});
I = find(x0);
H = zeros(numRxns,numRxns);
f = zeros(numRxns,1);    
for j = 1:length(I)
       H(I(j),I(j)) = 2/x0(I(j))^2;
       f(I(j)) = -2/x0(I(j)); 
end
[v, obj, exit] = quadprog(H,f,A,b,Aeq,beq,[],[],x0,options);
Rnew = (2*v(1)+0.6391*v(16))/v(34);

it = 0;
maxIt = 300;
damp = 0.05;
while abs(Rnew-R)/abs(Rnew)>Rtol & it<maxIt
    it = it + 1;
    R = R + damp*(Rnew-R);
    y = [C_new; shift;  R; cellType];
    inputs = mat2cell(y,ones(1,length(y)),1);
    x0 = rateLaws(inputs{:});
    I = find(x0);
    H = zeros(numRxns,numRxns);
    f = zeros(numRxns,1);    
    for j = 1:length(I)
           H(I(j),I(j)) = 2/x0(I(j))^2;
           f(I(j)) = -2/x0(I(j)); 
    end
    [v, obj, exit] = quadprog(H,f,A,b,Aeq,beq,[],[],x0,options);
    Rnew = (2*v(1)+0.64*v(16))/v(34).
end    
if it == maxIt
    R = Rold;
    y = [C_new; shift;  R; cellType];
    inputs = mat2cell(y,ones(1,length(y)),1);
    x0 = rateLaws(inputs{:});
    I = find(x0);
    H = zeros(numRxns,numRxns);
    f = zeros(numRxns,1);    
    for j = 1:length(I)
           H(I(j),I(j)) = 2/x0(I(j))^2;
           f(I(j)) = -2/x0(I(j)); 
    end
    [v, obj, exit] = quadprog(H,f,A,b,Aeq,beq,[],[],x0,options);
end
   
Rvec = [Rvec; Rnew];
objVec = [objVec; obj];
exitVec = [exitVec;exit];
ttestVec = [ttestVec;t];
save(fileName, 'Rvec','objVec','exitVec','ttestVec');

componentRates = stoichMatrix*v  * (C_new(11)*VCD/2.31)/1000; 

componentRates = componentRates(I_external);
end