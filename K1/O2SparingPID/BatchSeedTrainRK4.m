clear;close all;clc;                   % Clear out current variable list
% Parameters guess
load('finalParameters_v1.mat')

critDensity = 10;
system = 0;                         % 1 = perfusion, 0 = batch


if system == 0                      % Vessel sizes for the batch seed train
    vesselSize =  [1,4,15,60,250,1000,4E3,20E3,80E3,400E3];
elseif system ==1                   % Vessel sizes for the perfusion seed train
    vesselSize =  [1,4,15,60,250,1000,4E3,20E3,100E3];
end

rxnContainer = containers.Map('KeyType','double','ValueType','any');
cContainer = containers.Map('KeyType','double','ValueType','any');
tContainer = containers.Map('KeyType','double','ValueType','any');
Rcontainer = containers.Map('KeyType','double','ValueType','any');
extentContainer = containers.Map('KeyType','double','ValueType','any');
%% Inputs %%
tend = 10;                          % Defines the length of the seed train
cellType = 1;                       % Identifies the cell type (1=A)
shift = 0;                          % Initializes the temperature shift bool
Perfusion = 0;                      % Identifies the type of reactor (batch)
shiftDay = 100;                     % Day of the temperature shift
h = 0.005;                          % Step size for Runge Kutta 4th order in days
writeFile = 1;
%% Load variables from definitions
reversibleLogicals;                 % Loads logic for reversible/irreversible reactions
internalLogicals;                   % Loads logic for internal/external components
initialConditionsSeedTrain;         % Loads initial media conditions
initialConditionsCryo;              % Loads initial cryovial conditions
load('stoichMatrix.mat')            % Loads the stoichiometry matrix
t = 0:h:tend;                       % Defines timescale for seed train

I = find(~internal);                % Finds external components
C0 = initialConditionsCryo_vec(I);  % Finds initial conditions for the cryovial     

C_total = [C0];                             %
t_total = [0];                              %
R_total = [initialConditions_vec(end)];     %
obj_total = [-12];                          %
exit_total = [1];                           %
v_total = zeros(34,1);                      %
VCD_total = [C_total(5)/2.31];              %
vesselSize_total = [vesselSize(1)];         %
Tend = [];                                  %


% Begins iteration on each vessel size. Each iteration corresponds to a
% single vessel in the seed train.
for vessel = 2:length(vesselSize) 
    % The next line of code calculations the dilution of transitioning from
    % one vessel to the next. This is done in every step starting from the
    % cryovial.
    C0 =   (C_total(:,end) * vesselSize(vessel-1) + initialConditions_vec(I) *...
                    (vesselSize(vessel) - vesselSize(vessel-1))) / vesselSize(vessel);
    %        
    C0(5) = VCD_total(end)*vesselSize(vessel-1)/vesselSize(vessel)*2.31;            
   
    C = zeros(length(C0),length(t)); % Initializes concentration matrix.
    R = zeros(1,length(t));          % Initializes the redox variable vector
    objVec = R;                      % Initializes objective function vector 
    exitVec = R;                     % Initializes the exitflags vector
    VCD = R;                         % Initializes viable cell density vector
    v = zeros(34,length(t));         % Initializes the mass flux vector for each component
    % Initializes the vessel size vector. 
    vesselSize_vec = vesselSize(vessel)*ones(1,length(t));
    
    C(:,1) = C0;                    % Defines initial media concentrations in the solution matrix
    R(1) = R_total(end);            % Defines the intial value of the redox variable
    exitVec(1) = 1;                 % Sets the initial exit flag as success. This is done to help sizing issues.
    objVec(1) = -12;                % Sets the initial objective to a value. This is done to help sizing issues.
    VCD(1) = C0(5)/2.31;            % Sets the intial viable cell density (cells/mL)
    
    
    %% Runge Kutta (4th Order)
    i = 1;                          % Sets iterant for following while loop
    while VCD(i) < 10               % While the VCD is below the transfer value of 10, growth is allowed to proceed.
    % The line below performs the first step of Runge Kutta (4th order).
    % This calls the function instantRatesV2, which calculates reaction rates
    % at the current timestep.
        [K1, R1,~,~,~] =  instantRatesV2(C(:,i),t(i),Perfusion,cellType, shift, R(i), parameters,shiftDay);
        K1 = K1*h;
    % Second step of Runge Kutta
        [K2, R2,~,~,~] = instantRatesV2(C(:,i)+K1/2, t(i) + h/2 ,Perfusion,cellType, shift, R1, parameters,shiftDay);
        K2 = K2*h;
    % Third step of Runge Kutta
        [K3, R3,~,~,~] = instantRatesV2(C(:,i)+K2/2, t(i) + h/2 ,Perfusion,cellType, shift, R2, parameters,shiftDay);
        K3 = K3*h;
    % Fourth step of Runge Kutta
        [K4, R4,~,~,~] = instantRatesV2(C(:,i)+K3, t(i) + h ,Perfusion,cellType, shift, R3, parameters, shiftDay);
        K4 = K4*h;
    % The next line calculates the final reaction rate for each component
        rates = (K1+2*K2+2*K3+K4)/6;
        Rtest = (R1+2*R2+2*R3+R4)/6;
    % The next line calculates the concentration at the next step
        C(:,i+1) = C(:,i) + rates; 
    
   [~, R(i+1),v(:,i+1),objVec(i+1),exitVec(i+1)] =...
       instantRatesV2(C(:,i+1), t(i) + h ,Perfusion,cellType, shift, Rtest, parameters, shiftDay);

    %
        VCD(i+1) = C(5,i+1)/2.31 * (1-0.5259/(1+353.3*exp(-0.9381*t(i+1))));

        i = i + 1;                  % Itera
    end
C_total = [C_total,C(:,2:i)];
t_total = [t_total, t_total(end) + t(2:i)];
R_total = [R_total, R(2:i)];
obj_total = [obj_total,objVec(2:i)];
exit_total = [exit_total, exitVec(2:i)];
v_total = [v_total, v(:,2:i)];
VCD_total = [VCD_total,VCD(2:i)]; 
vesselSize_total = [vesselSize_total, vesselSize_vec(2:i)];
Tend = [Tend;t(i)];   

t2 = t;
t = t(2:i);
C = C(:,2:i);
R = R(2:i);
reactionScraperBatch;

rxnContainer(vessel) = rxn;
cContainer(vessel) = C;
tContainer(vessel) = t;
Rcontainer(vessel) = R;
extentContainer(vessel) = extent;

t = t2;
if (vessel >= 8 && system ==0) || (vessel == 9 && system ==1);
    t2 = t;
    t = t(2:i);
    % Identifies the reactor type for saving the file properly.
    if system ==0
        systemString = 'Batch';
    elseif system ==1
        systemString = 'Perfusion';
    end
    %
    Volume = num2str(vesselSize(vessel)/1000);
    fileName = [systemString,Volume,'.mat'];  
    Dt = (4*vesselSize(vessel)/(3*pi))^(1/3)/100; %m
    A = pi/4*Dt^2;
    C_O2_vec = smooth(C(14,:),15);
    C_CO2_vec = smooth(C(7,:),15);
    save(fileName,'t','shiftDay','C_CO2_vec','C_O2_vec','A','Dt','C','MAB_Produced','rxn','extent');
    t = t2;
end   
end


% The next 7 lines define variables in terms of variables found in
% seedTrainAnalyzer and plotToolV2. 
C = C_total;
t = t_total;
R = R_total;
obj = obj_total;
exit_vec = exit_total;
v = t_total;
VCD = VCD_total;

seedTrainAnalyzer
plotToolV2