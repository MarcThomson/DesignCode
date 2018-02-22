clear;close all;
critDensity = 10; % VCD, 10^6 cell/mL
%vesselSize=[10,50,250,1E3,20E3,80E3,400E3,2000E3]; %mL
vesselSize =  [40,180,850,4E3,20E3,80E3,400E3,2000E3];%[10,57.1860367967821,327.024280452292,1870.12225353860,10694.4880005339,61157.5384321278,349735.724318027,2000000.00000000]; %mL
vesselSize = vesselSize(2:end);
vessel = 1; %inital vessel

    
%% Inputs %%
cellType = 1;                                                              % Cells of Type 1=A and 0=B can be selected.
dt = 0.05;                                                                 % days. Timestep for iteration.
tend = 10;                                                                 % days. Total days for calculation.
shiftDay=11;                                                                % days. The day of the temperature shift where the temperature is slightly raised.


totalStructure=containers.Map('KeyType','double','ValueType','any');

%% Load variables from definitions
initialConditions;
reversibleLogicals;
parameterDefinitions;
internalLogicals;
load('stoichMatrix.mat')

t = zeros(1,tend/dt+1);
steps = length(t)-1;                                                       % Total number of iterative steps.
numRxns=length(stoichMatrix(1,:));                                         % Total number of reactions.
numComponents = length(stoichMatrix(:,1));                                 % Total number of components.

C = zeros(numComponents+2,length(t));                                      % Adds R and shift
initialConditions_vec2 = initialConditions_vec;
initialConditions_vec2(11)=4.5/vesselSize(1)*25*2.31*0.8;                                      % Defining the time interval for iteration
for vessel=1:length(vesselSize)
t = 0:dt:tend;
if vessel>1
       initialConditions_vec(11) = 0;
       Cold = totalStructure(vessel-1);
       told = totalStructure(-(vessel-1));
       initialConditions_vec2 =   (Cold(:,end) * vesselSize(vessel-1) + initialConditions_vec *...
                            (vesselSize(vessel) - vesselSize(vessel-1))) / vesselSize(vessel);
       initialConditions_vec2(11) = VCD(told(end))*initialConditions_vec2(11);
       
       initialConditions_vec2(end-1:end) =  Cold(end-1:end,end);
       
end
%% Optimization Settings


VCD = @(t1)(1-1./(1+exp(-0.5*(t1-9))));                                    % Defining the fractional viable cell density
tol = 1e-10;                                                               % Tolerance for nonlinear optimization



v = zeros(numRxns,length(t)-1);                                            % Initializes the rates of each reaction for each timestep.

C(:,1) = initialConditions_vec2;                                            % Adds initial conditions to the first column of the solution matrix.
objectives=zeros(steps,1);                                                 % Initializes objectives. These are minimized in the nonlinear optimization.
exitFlagVec=zeros(steps,1);                                                % Initialilzes exitflags. These signify if the optimization ran into an error.

 

pVCD = zeros(length(t),length(vesselSize)); %???
i = 1;
%% Begin Iteration
while pVCD(i,vessel) < critDensity
   
    if i==1                                                                % Load inputs
        y = [initialConditions_vec2; cellType];                             % If it is the first step, inputs are loaded from the initial conditions and respective celltype.
    else
        y = [C(:,i); cellType];                                            % If it is any other step, inputs are loaded from the previous step.
    end
    
    inputs = num2cell(y,length(y));                                        % Converts inputs to a usable format.
    rates = rateLaws(inputs{:});                                           % Calculates rates for the input conditions.
    x0 = rates;                                                            %
    
    %% Nonlinear Optimization Initializations
    
    A = zeros(length(x0),length(x0));
    for j=flip(1:length(x0))
        if reversible(j)==1
            A=[A(1:j-1,:);A(j+1:end,:)];
        else
            A(j,j) = -1;
        end
    end
    
    Aeq = stoichMatrix;                                                    % Constraint for optimization
    I_external=find(1-internal);                                           % We are only assuming pseudosteady-state for intracellular components.
    for j=flip(1:length(I_external))
        Aeq=[Aeq(1:I_external(j)-1,:);Aeq(I_external(j)+1:end,:)];
    end
    
    H = zeros(length(x0),length(x0));
    f = zeros(length(x0),1);
    I = find(rates);
    for j = 1:length(I)
        H(I(j),I(j)) = 2/x0(I(j))^2;
        f(I(j)) = -2/x0(I(j));
    end
    
    beq = zeros(length(Aeq(:,1)),1);
    b = zeros(length(A(:,1)),1);
    
    %% Nonlinear Optimization Execution
    options = optimoptions('quadprog','OptimalityTolerance',tol,...        % Define options for quadratic programming
        'MaxIterations',30000,'Algorithm','interior-point-convex',...
        'StepTolerance',1e-12,'Display','off');                     
    
    [v(:,i),objectives(i),exitFlagVec(i)]=...                              % Execute the quadratic optimization
        quadprog(H,f,A,b,Aeq,beq,[],[],[],options);
    
    R = (2*v(1,i)+0.64*v(16,i))/v(34,i);                                   % Define the R term in terms of current rates
    C(1:end-2,i+1) =  C(1:end-2,i) +  ...                                  % Updating cell density
    (stoichMatrix*v(:,i))*dt* VCD(t(vessel))/2.3*C(11,i)/1000;                  % 
    C(end-1,i+1) = C(end-1,i);
    C(end, i+1) = R;                                                       % Store R in solution matrix

    
    
    % Glucose monitoring 
    
    if t(i)-floor(t(i))==0 && C(20,i)<=40                                % If the sugar concentration drops below the threshold of 40, more sugar is added. 
        C(20,i+1)=C(20,i)+5;
    end
    
    pVCD(i,vessel)=VCD(t(i))*C(11,i)/2.31;
    
    if pVCD(i,vessel) > critDensity && vessel~=length(vesselSize)
       totalStructure(vessel)=C(:,1:i+1);
       totalStructure(-vessel)=t(1:i+1);
       break
    elseif t(i+1)==tend
       totalStructure(vessel)=C(:,1:i+1);
       totalStructure(-vessel)=t(1:i+1);
       %vessel=vessel+1; 
       break
    end
    i = i + 1;
end
end


