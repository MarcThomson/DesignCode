% Modeling File

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

%tol = 1e-10;
%options = optimoptions('fmincon','ConstraintTolerance',tol,'StepTolerance',tol,'OptimalityTolerance',tol);
tol = 1e-10;
options = optimoptions('fmincon','ConstraintTolerance',tol,'OptimalityTolerance',tol,'MaxFunctionEvaluations',30000,'display','off','Algorithm','sqp');


steps = length(t)-1;
numRxns=length(stoichMatrix(1,:));
numComponents = length(stoichMatrix(:,1));

v = zeros(numRxns,length(t)-1);
C = zeros(numComponents+2,length(t)); % Adds R and shift
C(:,1) = initialConditions_vec;
maxIndex=zeros(12,1);
objectives=zeros(steps,1);
%%

for i=1:length(t)-1
    if i==1
        y = [initialConditions_vec; cellType];
    else
        y = [C(:,i); cellType];
    end
    
    inputs = mat2cell(y,ones(1,length(y)),1);
    rates = rateLaws(inputs{:});
    x0 = rates;
    %% Nonlinear Optimization Initialization
    I = find(rates);
    Aeq = stoichMatrix;
    beq = zeros(numComponents,1);
    I_external=find(1-internal);       % We are only assuming pseudosteady-state for intracellular components.
    
    for j=1:length(I_external)
        Aeq(I_external(j),:) = 0;
    end
    
    A = zeros(length(x0),length(x0));
    for j=1:length(x0)
        A(j,j) = reversible(j)-1;
    end
    
    b = reversible;
    
    % Nonlinear Optimization Execution
    f_obj=@(x) sum( ((x(I)-x0(I))./(x0(I))).^2);
    %options = optimoptions('fmincon','MaxFunctionEvaluations',30000);
    v(:,i) = fmincon(f_obj,x0,A,b,Aeq,beq,[],[],[],options);
    R = (2*v(1,i)+0.64*v(16,i))/v(34,i);
    objectives(i)=f_obj(v(:,i));
    % f_obj_2=@(x)max( ((x(I)-x0(I))./(x0(I))).^2);
    % [objectives(i),maxIndex(i)]=f_obj_2(v(:,i));
    % Update Solution Matrix
    C(1:end-2,i+1) =  C(1:end-2,i) +  (stoichMatrix*v(:,i))*dt* VCD(t(i))/2.3*C(11,i)/1000; % Updating cell density
    C(34,i+1) = C(34,i);
    C(end-1,i+1) = C(end-1,i);
    C(end, i+1) = R;
    
    % Temperature Shift
    if t(i) > shiftDay
        C(end-1,i+1)=1;
    end
    
end
