load('finalParameters_v1.mat')
f_obj = batchFunction(parameters') 


function [f_obj] = batchFunction(parameters) 
%% Inputs %%
tend = 10;              % Length of batch reactor (days)
cellType = 1;           % Identifies cell type (1=A)
shift = 0;              % Initializes shift as false (bool)
Perfusion = 0;          % Initializes perfusion as false (bool)
h = 0.05;               % Stepsize in days
shiftDay = 3    ;       % Identifies the day of the temperature shift
Dt = (4*2000E3/(3*pi))^(1/3)/100;          % Diameter of tank used in batch (m). Calculated based on volume and 3:1 D:h ratio.
A = Dt^2/4*pi;          % Tank cross sectional area in m^2
plotExp = 1;            % Boolean defining  whether or not the experimental data is to be plotted
%% Load variables from definitions

load('reversibleLogicals.mat');      % Loads reversible/irreversible logicals for metabolic reactions
load('internalLogicals.mat');        % Loads internal/external logicals for components
load('initialConditions_NolanLee.mat'); %Load the initial conditions measured by Nolan and Lee
load('stoichMatrix.mat')             % Loads the Stoichiometric Matrix
t = 0:h:tend;

I = find(~internal);                 % Identifies external components

C0 = initialConditions_vec(I);

C = zeros(length(C0),length(t));    % Intializes concentration matrix
R = zeros(1,length(t));             % Initializes redox variable vector
objVec = zeros(length(t),1);        % Initial objective function vector
exitVec = zeros(length(t),1);       % Initializes exitflag vector
v = zeros(34,length(t));            % Initializes mass flux vector 

R(1) = 6.2;         % Initial condition on redox variable
% Initial conditions on concentration matrix based on media composition
C(:,1) = C0;
exitVec(1) = 1;     % Initial exitflag set to 1 for sizing issues
objVec(1) = -12;    % Initial objective function set for sizing issues

failed = 0;  % initialize failed as 0
for i = 1:length(t)-1
    if t(i)>shiftDay
        shift = 1;
        load('stoichMatrix31.mat');
    end
    % The next few lines of code perform RK4 to find reaction rates   
   [K1, R1,~,~,~] =  instantRatesV2(C(:,i),t(i),Perfusion,cellType, shift, R(i), parameters, shiftDay);
   K1 = K1*h; 
   % Second step of RK4
   [K2, R2,~,~,~] = instantRatesV2(C(:,i)+K1/2, t(i) + h/2 ,Perfusion,cellType, shift, R1, parameters, shiftDay);
   K2 = K2*h; 
   % Third step of RK4
   [K3, R3,~,~,~] = instantRatesV2(C(:,i)+K2/2, t(i) + h/2 ,Perfusion,cellType, shift, R2, parameters, shiftDay);
   K3 = K3*h;
    % Fourth step of RK4  
   [K4, R4,~,~,~] = instantRatesV2(C(:,i)+K3, t(i) + h ,Perfusion,cellType, shift, R3, parameters, shiftDay);
   K4 = K4*h;   
    % Using the terms from RK4, the reaction rates for each component are
     % calculated. 
   rates = (K1+2*K2+2*K3+K4)/6;
   Rtest = (R1+2*R2+2*R3+R4)/6;
   % These rates are then added to the previous step's
   % concentration matrix and define this step's concentrations.
   C(:,i+1) = C(:,i) + rates;
   
   % the next value of R is determined by running another trial
   [~, R(i+1),v(:,i+1),objVec(i+1),exitVec(i+1)] =...
       instantRatesV2(C(:,i+1), t(i) + h ,Perfusion,cellType, shift, Rtest, parameters, shiftDay);
   
   
   % if it is a feeding day, then feed the cells
   if  t(i)-floor(t(i))==0 && t(i)>=3
       C(1,i+1) = C(1,i+1)*0.97 + 0.03*6.00; %ALA
       C(4,i+1) = C(4,i+1)*0.97 + 0.03*15.04; %ASP
       C(8,i+1) = C(8,i+1)*0.97 + 0.03*444.4; %GLC
       C(6,i+1) = C(6,i+1)*0.97 + 0.03*4.68; %CC 
       C(3,i+1) = C(3,i+1)*0.97 + 0.03*54.01; %ASN 
       C(11,i+1) = C(11,i+1)*0.97+0.03*6.00; %GLY 
       C(10,i+1) = C(10,i+1)*0.97+0.03*6.00; %GLU
       C(15,i+1) = C(15,i+1)*0.97+0.03*45.16; %SER 
       %ALA
   end
  
   I_not_gas = [1:6,8:13,15]; %components that are not O2/CO2
   Ineg = find(C(I_not_gas,i+1)<0);
   
   % fail the system if R is too large or any component is negative
    if abs(R(i+1))>10 || length(Ineg)>0
       R(i+1) = 0;
       failed = 1;
       break
   end
end

% if failed, give a high objective function
if failed ==1
    C = 0*C;
    C = C - 1000;
end

plotToolV2;

% compute the objective function (R^2, adj)
load('expData.mat');
I = [1;2;3;4;5;6;8;9;10;11;12;13;15]; % components to compare

% compute ss total
meanC = mean(mean(C_exp));
SStotal = sum(sum((C_exp-meanC).^2));

% compute SSR
SSR = 0;
for i = 1:length(I)
    C_pred = interp1(t,C(I(i),:),t_exp);
    SSR = SSR + sum((C_pred-C_exp(i,:)).^2);
end

% compute R2
R2 = 1-SSR/SStotal;

% compute R2, adj
n = length(C_exp(:));
R2adj = 1- (((1-R2)*(n-1))/(n-47-1));
f_obj  = 1 - R2adj;


end
