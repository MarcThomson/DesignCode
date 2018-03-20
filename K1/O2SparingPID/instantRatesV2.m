function [componentRates, R, v, obj, exit] = instantRatesV2(C,t,Perfusion,cellType,shift,R, parameters, shiftDay)
% outputs the reaction rate as a functon of current system parameters
% inputs:
% C =  vector of concentrations of all external components (15x1) in standard order
% t = current time (days)
% Perfusion = 1 if perfusion system, 0 otherwise
% cellType = cell line. A = 1, 0 = B. A is standard
% shift = 1 if the shift has occured, 0 otherwise
% R = current redox parameter
% paremeters = 47x1 vector of the non-logged kinetic parameters
% shiftDay = the day on which the shift would occur (days)
% outputs: 
% componentRates = 15x1 vector of the rate of appearance/disapearance for
% each component
% R = new redox parameter
% v = 34x1 vector of reaction fluxes 
% obj = objective function for MFB
% exitflag = flag produced by quadprog indicating how the program finished.
% -1 means a successful minimization

% set the VCD function
if Perfusion
    VCD = 1;
    carryingCapacity = 90; %10 ^6/mL
else
    VCD = 1-0.5259/(1+353.3*exp(-0.9381*t));
end

%% Call other parameters
parameterDefinitions;
reversibleLogicals;
internalLogicals;

% if the shift has occured, take a weighted average of the stoichiometric
% matrices
% these multipliers assume a roughly 3 hour transition time between the hot
% and cold temperatures
if shift == 1
    preMultiplier = exp(-33.9400974951428*(t-shiftDay));
    postMultiplier = 1-exp(-33.9400974951428*(t-shiftDay));
else
    preMultiplier = 1;
    postMultiplier = 0;
end

load('stoichMatrix.mat')
preStoich = stoichMatrix;
load('stoichMatrix31.mat')
postStoich = stoichMatrix;
stoichMatrix =  preMultiplier * preStoich + postMultiplier*postStoich;

numRxns = length(stoichMatrix(1,:));
numComponents = length(stoichMatrix(:,1));

% define a new vector including all components for MFB
I_external = find(internal == 0);
C_new = zeros(numComponents,1);
C_new(I_external) = C;

%% Optimization Initialization
tol = 1e-10;
options = optimoptions('quadprog','OptimalityTolerance',tol,...
    'MaxIterations',30000,'Algorithm','interior-point-convex',...
    'StepTolerance',1e-12,'Display','off');

% condition 1: Aeq * x = b
% this enforces pseudo-steady state for the internal components
Aeq = stoichMatrix; %initialize Aeq as the stoichiometric matrix, then trim
                    % off external components
for j=flip(1:length(I_external))
 Aeq=[Aeq(1:I_external(j)-1,:);Aeq(I_external(j)+1:end,:)];
end
% these rates should be 0
beq = zeros(length(Aeq(:,1)),1);

% condition 2: x(i)>=0 if the ith reaction is irreversible
% this just puts -1 along the diagonal, and enforces A*<=0
A = zeros(numRxns, numRxns);
for j=flip(1:numRxns)
    if reversible(j)==1
        A=[A(1:j-1,:);A(j+1:end,:)];
    else
        A(j,j) = -1;
    end
end
b = zeros(length(A(:,1)),1);


% inputs for initial rate calculations
y = [C_new; shift;  R; cellType];
inputs = mat2cell(y,ones(1,length(y)),1); %this is done to streamline the inputs
inputs = [inputs;{parameters};t;shiftDay]; 
x0 = rateLaws(inputs{:}); % kinetically estimated rates, nmol/ (10^6 cells)/day

% in the perfusion case, scale down the rate of cell growth by comparing
% the density to the carrying capacity
if Perfusion
    x0(16) = x0(16)*(1 - C(5)*VCD/2.31  / carryingCapacity );
end

% initialize the quadratic programming
% obj = 1/2 * x' * H * x + f'*x
% in this case, H(i,i) = 2/x0(i)^2 for i calculated
% f(i) = -2/x(i) for i calculated
I = find(x0); % nonzero rates 
H = zeros(numRxns,numRxns); 
f = zeros(numRxns,1);    
for j = 1:length(I)
       H(I(j),I(j)) = 2/x0(I(j))^2;
       f(I(j)) = -2/x0(I(j)); 
end

% if there is a numerical error in setting up H, don't run quadprog
if isequal(H,H') && sum(sum(imag(H)))==0
    % minimize the objective function
    [v, obj, exit] = quadprog(H,f,A,b,Aeq,beq,[],[],[],options);
else
    % output if it failed
    v = ones(34,1);
    obj = 1000;
    exit = -1;
end

% calculate the new value of R
R = (2*v(1) + 0.64*v(16))/v(34);

% calculate the total rates in mmol/day by multiplying by cell density and
% performing some unit conversions
componentRates = stoichMatrix* v  * (C_new(11)*VCD/2.31)/1000; 

% cut down the vector to the 15x1 external components
componentRates = componentRates(I_external);
end