% Modeling File

%% Inputs %%
clear all;close all;
cellType = 1;   % Type 1=A, 0=B
%dt = 0.25;      % days
tend = 10;      % days


%% Load variables from definitions

load('stoichMatrix.mat')
initialConditions;
tspan = [0 tend];
C0 = [initialConditions_vec;cellType;1];



[t, C] = ode113(@dCdt,tspan,C0);
load gong
sound(y,Fs)
%%

function derivativeSystem = dCdt(t,C)
    shiftDay=4;     % days
    cellType = C(end-1);
    parameterDefinitions;
    reversibleLogicals;
    load('stoichMatrix.mat')
    internalLogicals;
    numRxns=length(stoichMatrix(1,:));
    numComponents = length(stoichMatrix(:,1));
    VCD = @(t1)(1-1./(1+exp(-0.5*(t1-9))));

    
    if t > shiftDay
        C(end-3) = 1;
    end
    

    %% Nonlinear Optimization Initialization
    inputs = mat2cell(C(1:end-1),ones(1,length(C)-1),1);
    rates = rateLaws(inputs{:});
    x0 = rates;
    
    tol = 1e-10;
    options = optimoptions('fmincon','ConstraintTolerance',tol,'OptimalityTolerance',tol,'MaxFunctionEvaluations',30000,'display','off','Algorithm','sqp');   
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
    v = fmincon(f_obj,x0,A,b,Aeq,beq,[],[],[],options);


    derivativeSystem = stoichMatrix*v* VCD(t)/2.3*C(11)/1000; % Updating cell density
    derivativeSystem(34) = 0;
    
    
    
    R = (2*v(1)+0.64*v(16))/v(34);
    obj=f_obj(v);
    derivativeSystem= [derivativeSystem;0;100*(R-C(end-2));0;100*(obj-C(end))];

    
end
