% Modeling File

%% Inputs %%
close all; clear
cellType = 1;                                                              % Cells of Type 1=A and 0=B can be selected.
dt = 0.01;                                                                 % days. Timestep for iteration.
tend = 10;                                                                 % days. Total days for calculation.
shiftDay =  3;                                                                % days. The day of the temperature shift where the temperature is slightly raised.

%% Load variables from definitions
parameterDefinitions;
reversibleLogicals;
initialConditions;
internalLogicals;
load('stoichMatrix.mat')

%% Optimization Settings

t = dt:dt:tend;                                                             % Defining the time interval for iteration
VCD = @(t1)(1-1./(1+exp(-0.5*(t1-9))));                                    % Defining the fractional viable cell density
tol = 1e-15;                                                               % Tolerance for nonlinear optimization

steps = length(t)-1;                                                       % Total number of iterative steps.
numRxns=length(stoichMatrix(1,:));                                         % Total number of reactions.
numComponents = length(stoichMatrix(:,1));                                 % Total number of components.

v = zeros(numRxns,length(t)-1);                                            % Initializes the rates of each reaction for each timestep.
C = zeros(numComponents+2,length(t));                                      % Adds R and shift
C(:,1) = initialConditions_vec;                                            % Adds initial conditions to the first column of the solution matrix.
objectives=zeros(steps,1);                                                 % Initializes objectives. These are minimized in the nonlinear optimization.
exitFlagVec=zeros(steps,1);                                                % Initialilzes exitflags. These signify if the optimization ran into an error.
worst = [];
%% Begin Iteration
for i=1:length(t)-1
    %if t(i)>shiftDay
        %load('StoichMatrix31.mat');
    %end
    
    if i==1                                                                % Load inputs
        y = [initialConditions_vec; cellType];                             % If it is the first step, inputs are loaded from the initial conditions and respective celltype.
    else
        y = [C(:,i); cellType];                                            % If it is any other step, inputs are loaded from the previous step.
    end
    
    inputs = num2cell(y,length(y));                                        % Converts inputs to a usable format.
    rates = rateLaws(inputs{:});                                           % Calculates rates for the input conditions.
    x0 = rates;                                                            %
    
    %% Nonlinear Optimization Initializations
    
%     %%%%% test
%     initialFlux;
%      I = [1,2,3,8,9,10,11,12,13,16,17,33];
%     %%%%% end test
    
    
    A = zeros(length(x0),length(x0));
    for j=flip(1:length(x0))
        if reversible(j)==1
            A=[A(1:j-1,:);A(j+1:end,:)];
        else
            A(j,j) = -1;
        end
    end
    
    Aeq = stoichMatrix;                                                    % Constraint for optimization
    Aeq = Aeq(find(internal),:);
    %Aeq(18,:) = Aeq(18,:)+Aeq(19,:);
    %Aeq = Aeq([1:17,20:end],:);
    %Aeq = Aeq([1:8,10:17,20:end],:);
    %Aeq(end+1,end) = 1; %%%%% test
    H = zeros(length(x0),length(x0));
    f = zeros(length(x0),1);
    I = [1,2,3,8,9,10,11,12,13,16,17,33];
    for j = 1:length(I)
        H(I(j),I(j)) = 2/x0(I(j))^2;
        f(I(j)) = -2/x0(I(j));
    end
    
    beq = zeros(length(Aeq(:,1)),1);
    %beq(end) = 1000;
    
    b = zeros(length(A(:,1)),1);
    %testing!!!!
%     H = Aeq'*Aeq*2;
%     Aeq = zeros(12,34);
%     Aeq([I,I]) = 1;
%     beq = x0(I);
%     f = 0*f;
    %% Nonlinear Optimization Execution
    options = optimoptions('quadprog','OptimalityTolerance',tol,...        % Define options for quadratic programming
        'MaxIterations',90000,'Algorithm','interior-point-convex',...
        'StepTolerance',1e-12,'Display','off');                     
    initialFlux;
    [v(:,i),objectives(i),exitFlagVec(i)]=...                              % Execute the quadratic optimization
        quadprog(H,f,A,b,Aeq,beq,[],[],initialFluxVec,options);
    
    %% some diagnostic testing to get the right fluxes
    Itest = find(x0);
    p=((v(Itest,i)-x0(Itest))./x0(Itest)).^2;
    [Mtest, I2] = max(p);
    worst(i) = I2;
    initialFlux;
    diff = ((initialFluxVec-v(:,1))./initialFluxVec).^2;
    if v(34,1)<800
        denom = 800;
    else
        denom = v(34,1);
    end
    
    if t(i)>6
        2+2;
    end
    R = (2*v(1,i)+0.64*v(16,i))/v(34,i);      %5.889*exp(-0.6695*t(i))+0.1072*t(i);%                             % Define the R term in terms of current rates

    
    C(1:end-2,i+1) =  C(1:end-2,i) +  ...                                  % Updating cell density
    (stoichMatrix*v(:,i))*dt* VCD(t(i))/2.31*C(11,i)/1000;                  % 
    C(end-1,i+1) = C(end-1,i);
    C(end, i+1) = R;                                                       % Store R in solution matrix
    
    % Temperature Shift
    if t(i) > shiftDay
        C(end-1,i+1)=1;                                                    % If the current time is past the designated shift day, shift the temperature.
    end
         
    if t(i)-floor(t(i))==0 && C(20,i)<=40  && t(i)>=3                                % If the sugar concentration drops below the threshold of 40, more sugar is added. 
        C(20,i+1)=C(20,i)+10;
    end
end
plotTool                                                                   % Generates a figure for the species with designated importance.

