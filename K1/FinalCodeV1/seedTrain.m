clear;close all;clc;                   % Clear out current variable list
% load the optimized parameters
load('finalParameters_v1.mat')
%% Inputs %%
critDensity = 10;
system = 0;                         % 1 = perfusion, 0 = batch
tend = 10;                          % Defines the length of the seed train
cellType = 1;                       % Identifies the cell type (1=A)
shift = 0;                          % Initializes the temperature shift bool
Perfusion = 0;                      % Identifies the type of reactor (batch)
shiftDay = 100;                     % Day of the temperature shift
h = 0.05;                           % Step size for Runge Kutta 4th order in days
writeFile = 0;                      % Boolean to determine if an output file is written
plotExp = 0;                        % Boolean to determine if experimental data should be added to the plot
%% initialize the set of vessels 
% all vessel sizes in mL
if system == 0                      % Vessel sizes for the batch seed train
    vesselSize =  [1,4,15,60,250,1000,4E3,20E3,80E3,400E3];
elseif system ==1                   % Vessel sizes for the perfusion seed train
    vesselSize =  [1,4,15,60,250,1000,4E3,20E3,100E3];
end

%% Load variables from definitions
load('reversibleLogicals.mat');                 % Loads logic for reversible/irreversible reactions
load('internalLogicals.mat');                   % Loads logic for internal/external components
load('initialConditionsMedia.mat');             % Loads initial media conditions
load('initialConditionsCryovial.mat');          % Loads initial cryovial conditions
load('stoichMatrix.mat')            % Loads the stoichiometry matrix
t = 0:h:tend;                       % Defines timescale for seed train

I = find(~internal);                % Finds external components
C0 = initialConditionsCryo_vec(I);  % Finds initial conditions for the cryovial     

C_total = [C0];                             % Concentrations across all seed train vessels
t_total = [0];                              % time across all seed train vessels
R_total = [initialConditions_vec(end)];     % R across all seed train vessels
obj_total = [-12];                          % objevtice function across all seed train vessels
exit_total = [1];                           % exitflag across all seed train vessels
v_total = zeros(34,1);                      % fluxes across all seed train vessels
VCD_total = [C_total(5)/2.31];              % VCD across all seed train vessels
vesselSize_total = [vesselSize(1)];         % Total vessel sizes across all seed train vessels
Tend = [];                                  % Time for each seed train vessel


% initialize containers to store output data that will be written to a file
rxnContainer = containers.Map('KeyType','double','ValueType','any'); % rxn String
cContainer = containers.Map('KeyType','double','ValueType','any'); % concentrations 
tContainer = containers.Map('KeyType','double','ValueType','any'); % times
Rcontainer = containers.Map('KeyType','double','ValueType','any'); % R
extentContainer = containers.Map('KeyType','double','ValueType','any'); % extent of reaction

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
    
   % R is calculated at the new point with an additional iteration
   [~, R(i+1),v(:,i+1),objVec(i+1),exitVec(i+1)] =...
       instantRatesV2(C(:,i+1), t(i) + h ,Perfusion,cellType, shift, Rtest, parameters, shiftDay);

    % Viable cell density is updates based on cell concentration and an emperical functtion    
        VCD(i+1) = C(5,i+1)/2.31 * (1-0.5259/(1+353.3*exp(-0.9381*t(i+1))));

        i = i + 1;                  % Iteration
    end
C_total = [C_total,C(:,2:i)]; %total concentrations expanded
t_total = [t_total, t_total(end) + t(2:i)]; %total times expanded 
R_total = [R_total, R(2:i)]; %total R expanded 
obj_total = [obj_total,objVec(2:i)]; %total objectives expanded
exit_total = [exit_total, exitVec(2:i)]; %total exit vector expanded 
v_total = [v_total, v(:,2:i)]; % total fluxes expanded 
VCD_total = [VCD_total,VCD(2:i)]; % total VCD expanded 
vesselSize_total = [vesselSize_total, vesselSize_vec(2:i)]; % total vessel size expanded
Tend = [Tend;t(i)];  % end time expanded

t2 = t; %backup t
t = t(2:i); %shorten t, C, R to truncate later times
C = C(:,2:i);
R = R(2:i);
% output the reaction from the current vessel
reactionAnalyzerBatch;

% save the current vessel properties to the containers
rxnContainer(vessel) = rxn;
cContainer(vessel) = C;
tContainer(vessel) = t;
Rcontainer(vessel) = R;
extentContainer(vessel) = extent;

% restore t
t = t2;

% if a vessel that requires a spargers, save the data
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
    A = pi/4*Dt^2; %m^2
    C_O2_vec = smooth(C(14,:),15); %smooth to erase numerical error
    C_CO2_vec = smooth(C(7,:),15);
    % save relevent data
    save(fileName,'t','shiftDay','C_CO2_vec','C_O2_vec','A','Dt','C',...
        'MAB_Produced','rxn','extent');
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

% produce some plots and output data
seedTrainAnalyzer

% plot the data
plotToolV2