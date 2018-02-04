
%% Inputs %%
clear all;close all;
cellType = 1;   % Type 1=A, 0=B
dt = 0.25;      % days
tend = 10;      % days
shiftDay=4;     % days

%% Load variables from definitions
parameterDefinitions;
reversibleLogicals;
initialConditions;
internalLogicals;
load('stoichMatrix.mat')

%% Optimization Settings
t = 0:dt:tend;
VCD = @(t1)(1-1./(1+exp(-0.5*(t1-9))));
tol = 1e-10;
options = optimoptions('quadprog','OptimalityTolerance',tol,'MaxIterations',30000,'Algorithm','interior-point-convex','StepTolerance',1e-12);
    
%% Initializations
steps = length(t)-1;
numRxns=length(stoichMatrix(1,:));
numComponents = length(stoichMatrix(:,1));
v = zeros(numRxns,length(t)-1);
C = zeros(numComponents+2,length(t)); % Adds R and shift
C(:,1) = initialConditions_vec;
maxIndex=zeros(12,1);
objectives=zeros(steps,1);
exitFlagVec=zeros(steps,1);

%% Iterations for Feasible Solution
for i=1:length(t)-1
    
    if i==1
        y = [initialConditions_vec; cellType];
    else
        y = [C(:,i); cellType];
    end
    
    inputs = num2cell(y,length(y));
    rates = rateLaws(inputs{:});
    x0 = rates;
%% Nonlinear Optimization Initialization
    I_external=find(1-internal);       % We are only assuming pseudosteady-state for intracellular components.
    
    Aeq = stoichMatrix;
    for j=1:length(I_external)
        Aeq(I_external(j),:) = 0;
    end
    
    for j=flip(1:length(I_external))
        Aeq=[Aeq(1:I_external(j)-1,:);Aeq(I_external(j)+1:end,:)];
    end  
    
    A = zeros(length(x0),length(x0));  
    for j=flip(1:length(x0))
        if reversible(j)==1
            A=[A(1:j-1,:);A(j+1:end,:)];
        else
            A(j,j) = -1;
        end
    end

    I = find(rates);
    H = zeros(length(x0),length(x0));
    f = zeros(length(x0),1); 
    for j = 1:length(I)
       H(I(j),I(j)) = 2/x0(I(j))^2;
       f(I(j)) = -2/x0(I(j)); 
    end
    b = zeros(length(A(:,1)),1);
    beq = zeros(length(Aeq(:,1)),1);
    
%% Nonlinear Optimization Execution Initialization
    % Execution of QuadProg solver
    [v(:,i), objectives(i), exitFlagVec(i)] =quadprog(H,f,A,b,Aeq,beq)    ;
    R = (2*v(1,i)+0.64*v(16,i))/v(34,i);
    % Update Concentration Solution Matrix
    C(1:end-2,i+1) =  C(1:end-2,i) +  (stoichMatrix*v(:,i))*dt* VCD(t(i))/2.3*C(11,i)/1000; % Updating cell density
    % Set O2 concentration constant.
    C(34,i+1) = C(34,i);
    % Set variable 34 to be constant
    C(end-1,i+1) = C(end-1,i);
    % Log R value
    C(end, i+1) = R;  
    % Temperature Shift Handeling
    if t(i) > shiftDay
        C(end-1,i+1)=1;
    end
    
end
