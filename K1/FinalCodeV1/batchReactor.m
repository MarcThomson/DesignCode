% Batch Reactor

% This code simulates the dynamic behavior of CHO-K1 metabolism in a batch
% reactor. The model for CHO-K1 metabolism model was created with the help
% of Lee et al. (2011).

% Written by:
% Edward DeCrescenzo
% Marc Thomson


clear;close all;clc;                   % Clear out current variable list

% Parameters 
load('finalParameters_v1.mat')
%% Inputs %%
tend = 10;              % Length of batch reactor (days)
cellType = 1;           % Identifies cell type (1=A)
shift = 0;              % Initializes shift as false (bool)
Perfusion = 0;          % Initializes perfusion as false (bool)
h = 0.005;              % Stepsize in days
shiftDay = 4    ;       % Identifies the day of the temperature shift
Dt = (4*2000E3/(3*pi))^(1/3)/100;          % Diameter of tank used in batch (m). Calculated based on volume and 3:1 D:h ratio.
A = Dt^2/4*pi;          % Tank cross sectional area in m^2
writeFile =0;          % Boolean to determine if the data is written to a file
plotExp = 0;            % Boolean defining  whether or not the experimental data is to be plotted
%% Load variables from definitions
load('reversibleLogicals.mat');      % Loads reversible/irreversible logicals for metabolic reactions
load('internalLogicals.mat');        % Loads internal/external logicals for components
load('initialConditionsMedia.mat');  % Loads initial condtions of the media
load('initialConditionsBatch.mat');  % Loads the initial conditions from the seed train
load('stoichMatrix.mat')             % Loads the Stoichiometric Matrix
t = 0:h:tend;

I = find(~internal);                 % Identifies external components

% get the initial conditions in the reactor taking into account mixing of
% media and seed train output
C0 = ((A*3*Dt * 0.79-0.4)*initialConditions_vec(I) + 0.4*C0)/(A*3*Dt * 0.79);

C = zeros(length(C0),length(t));    % Intializes concentration matrix
R = zeros(1,length(t));             % Initializes redox variable vector
objVec = zeros(length(t),1);        % Initial objective function vector
exitVec = zeros(length(t),1);       % Initializes exitflag vector
Volume = [];                        % Initializes a vector containing the variable reactor volume over time
v = zeros(34,length(t));            % Initializes mass flux vector 

R(1) = 6.2;         % Initial condition on redox variable
% Initial conditions on concentration matrix based on media composition
C(:,1) = C0;
exitVec(1) = 1;     % Initial exitflag set to 1 for sizing issues
objVec(1) = -12;    % Initial objective function set for sizing issues
Volume(1) = A*3*Dt * 0.79;  % Initial fluid volume in the reactor
totalVolume =  A*3*Dt; % Total reactor volume used for dilution

totalFeed = zeros(15,7); % Initializes the feed vector for periodic feeding, mmol

for i = 1:length(t)-1
    if t(i)>shiftDay % If the shift day has passed,
        shift = 1;   % then set shift to TRUE
        load('stoichMatrix31.mat'); % and load the new stoichmatrix
    end
    
    % The next few lines of code perform RK4 to find reaction rates   
    [K1, R1,~,~,~] =  instantRatesV2(C(:,i),t(i),Perfusion,cellType, shift, R(i), parameters,shiftDay);
    K1 = K1*h;
    % Second step of RK4
    [K2, R2,~,~,~] = instantRatesV2(C(:,i)+K1/2, t(i) + h/2 ,Perfusion,cellType, shift, R1, parameters,shiftDay);
     K2 = K2*h;
   % Third step of RK4
    [K3, R3,~,~,~] = instantRatesV2(C(:,i)+K2/2, t(i) + h/2 ,Perfusion,cellType, shift, R2, parameters,shiftDay);
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
   if  t(i)-floor(t(i))==0 && t(i)>=3 && t(i)<10
       Volume = [Volume;Volume(end) + 0.03*A*3*Dt]; % increase the volume by 3%
       I_not_gas = [1:6,8:13,15]; %components that are diluted; O2/CO2 are assumed to remaind constant during feeding
       C(I_not_gas,i+1) = C(I_not_gas,i+1)*Volume(end-1)/Volume(end); %dilute components
       Feed = zeros(15,1); %concentration in the 3% by volume solution
       Feed(1) = 6.00;%ALA
       Feed(3) = 60; %ASN
       Feed(4) = 15.04; %ASP
       Feed(6) = 4.68; %CC
       Feed(8) = 444.4; %GLc
       Feed(9) = 14; %GLN
       Feed(10) = 6;% GLU
       Feed(11) = 6; %GLY
       Feed(15) = 45.16; %SER
       
       % Add the feeds
       C(:,i+1) = C(:,i+1) + 0.03*totalVolume/Volume(end)*Feed;
       
       % Calculated the total amount fed
       totalFeed(:,floor(t(i))-2) = 0.03*totalVolume*Feed;
   end 
end

plotToolV2;   % Plots component concentrations ove time

% analyzes the output to get the reaction, and other information
reactionAnalyzerBatch

% if desired, this data is saved. O2/CO2 are smoothed to erase numerical
% error
if writeFile                                
    C_O2_vec = smooth(C(14,:),15);
    C_CO2_vec = smooth(C(7,:),15);
    fileName = 'Batch_Output.mat';
    save(fileName,'t','shiftDay','C_CO2_vec','C_O2_vec','A','Dt','C','MAB_Produced','rxn','extent');
end
