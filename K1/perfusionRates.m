function componentRates = perfusionRates(C)
%% Extra inputs
cellType = 1;
Rguess = 6.2;
shift = 0;
R_tol = 0.001;
VCD = 0.8;

%% Call other parameters
parameterDefinitions;
reversibleLogicals;
internalLogicals;
load('stoichMatrix.mat')

numRxns=length(stoichMatrix(1,:));
numComponents = length(stoichMatrix(:,1));


%%Modify concentrations
I_external = find(internal ==0);
C_new = zeros(numComponents,1);
C_new(I_external)= C;

%% Optimization Initialization
tol = 1e-10;
options = optimoptions('quadprog','OptimalityTolerance',tol,'MaxIterations',30000,'Algorithm','interior-point-convex','StepTolerance',1e-12,'Display','off');

Aeq = stoichMatrix;
for j=flip(1:length(I_external))
    %Aeq(I_external(j),:) = 0;
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

%% First Run
y = [C_new; shift;  Rguess; cellType];
inputs = mat2cell(y,ones(1,length(y)),1);
x0 = rateLaws(inputs{:});
I = find(x0);
H = zeros(numRxns,numRxns);
f = zeros(numRxns,1);    
for j = 1:length(I)
       H(I(j),I(j)) = 2/x0(I(j))^2;
       f(I(j)) = -2/x0(I(j)); 
end
v = quadprog(H,f,A,b,Aeq,beq,[],[],x0,options);
Rnew = (2*v(1)+0.64*v(16))/v(34);
Rold = Rguess;


%% Minimize R difference
it = 1;
while abs((Rnew-Rold)/Rold) > R_tol
    it = it+1;
    Rold = Rnew;
    y = [C_new; shift; Rold; cellType];
    inputs = mat2cell(y,ones(1,length(y)),1);
    x0 = rateLaws(inputs{:});
    I = find(x0);
    H = zeros(numRxns,numRxns);
    f = zeros(numRxns,1);    
    for j = 1:length(I)
           H(I(j),I(j)) = 2/x0(I(j))^2;
           f(I(j)) = -2/x0(I(j)); 
    end
    v = quadprog(H,f,A,b,Aeq,beq,[],[],x0,options);    
    Rnew = (2*v(1)+0.64*v(16))/v(34);
end
    
componentRates = stoichMatrix*v  * (C_new(11)*VCD/2.31)   /1000; 

componentRates = componentRates(I_external);
end