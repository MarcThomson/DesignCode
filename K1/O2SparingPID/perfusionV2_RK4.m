clear;                   % Clear out current variable list
%load('ObjFun.mat');     % Fits for all 47 parameters
load('R2.mat')           % Fits (based on R2) for all 47 parameters

% Parameters guess
parameters = x;
parameters = exp(parameters);
%parameters  = log(x0);

%% Inputs %%
tend = 30;                                % Ending day of perfusion (days)
cellType = 1;                             % Cell type
shift = 0;                                % Initialize shift variable (bool)
Perfusion = 1;                            % Type of reactor initialized
h = 0.01;                                 % Step size for Runge Kutta 4 (days)
alpha = (1+19.2703)*(1-0.99975);          % Cell death parameter. A measure of cells killed in the filter.
tau = 1;                                  % Reactor residence time (1 day)
shiftDay = 6;                             % Day of temperature shift (days)
Dt = (4*500E3/(3*pi))^(1/3)/100;          % Diameter of tank used in perf (m). Calculated based on volume and 3:1 D:h ratio.
A = pi/4*Dt^2;                            % Area of tank base based on Dt.
writeFile=0;                              % Should the solution be written to a file?
%% Load variables from definitions
reversibleLogicals;                       % Define reversible/irreversible reactions
internalLogicals;                         % Define internal/external reactions
load('stoichMatrix.mat')                  % Define reaction stoichiometry
initialConditionsSeedTrain;               % Consider renaming, media concentrations.

load('initialConditionsPerfusion.mat')    % Output from perf seed train

I = find(~internal);                      % Identify external components
C_i_in = initialConditions_vec(I);        % Isolate concentrations of external components
C0 = ( 400 * C_i_in + 100*C0)/(500);      % Find new concentrations from dilution of seed train 
t = 0:h:tend;

C0 = [C0;C0(5)];                          %?????

C = zeros(length(C0),length(t));          % Initialize solution matrix
C(:,1) = C0;                              % Initial conditions on media
R = zeros(1,length(t));                   % Initialize redox variable vector
R(1) = 6.2;                               % Initial condition on R
objVec = R;                               % Initialize objectives vector for optimization
objVec(1) = -12;                          % For uniform vector sizing an initial objVec is defined
exitVec = R;                              % Initialize exitFlags for optimization
exitVec(1) = 1;                           % For uniform vector sizing an initial exitFlag is defined 
v = zeros(34,length(t));                  % Initialize mass flux vector
I_notBIOM = [1:4,6:length(C0)-1];         % 

for i = 1:length(t)-1                     % Start loop on time (days)
    
    if t(i)>shiftDay                      % Check if shiftDay>=currentDay
        shift = 1;                        % Set shift true
        load('stoichMatrix31.mat');       % Load the shifted stoich matrix    
    end
    
  % Calls the instant rates file which generates production/use rates for
  % the current timestep. 
   [K1, R1,~,~,~] =  instantRatesV2(C(1:end-1,i),t(i),Perfusion,cellType, shift, R(i), parameters, shiftDay);
  % K1 is rate of each component. Below is a mass balance on each component
  % in the perfusion system. This excludes biomass as biomass is handled
  % separately in the next line. 
   K1(5) = K1(5) - alpha/tau*C(5,i);
   K1 = K1*h;   % Generates the actual component changes based on timestep.
  % This is the first term generated for Runge Kutta 4
   
   [K2, R2,~,~,~] = instantRatesV2(C(1:end-1,i)+K1(1:end-1)/2, t(i) + h/2 ,Perfusion,cellType, shift, R1, parameters, shiftDay);
   K2(I_notBIOM) = K2(I_notBIOM) - C(I_notBIOM,i)/tau + C_i_in(I_notBIOM)/tau;
   K2 = [K2; K2(5)];
   K2(5) = K2(5) - alpha/tau*C(5,i);
   K2 = K2*h; 
  % This is the second term generated for Runge Kutta 4
   
   [K3, R3,~,~,~] = instantRatesV2(C(1:end-1,i)+K2(1:end-1)/2, t(i) + h/2 ,Perfusion,cellType, shift, R2, parameters, shiftDay);
   K3(I_notBIOM) = K3(I_notBIOM) - C(I_notBIOM,i)/tau + C_i_in(I_notBIOM)/tau;
   K3 = [K3; K3(5)];
   K3(5) = K3(5) - alpha/tau*C(5,i);
   K3 = K3*h;
  % This is the third term generated for Runge Kutta
   
   [K4, R4,~,~,~] = instantRatesV2(C(1:end-1,i)+K3(1:end-1), t(i) + h ,Perfusion,cellType, shift, R3, parameters, shiftDay);
   K4(I_notBIOM) = K4(I_notBIOM) - C(I_notBIOM,i)/tau + C_i_in(I_notBIOM)/tau;
   K4= [K4; K4(5)];
   K4(5) = K4(5) - alpha/tau*C(5,i);
   K4 = K4*h; 
  % This is the fourth and final term generated for Runge Kutta 4
   
   Rtest = (R1+2*R2+2*R3+R4)/6;
   rates = (K1+2*K2+2*K3+K4)/6;
  
   
   C(:,i+1) = C(:,i) + rates;
   
   [~, R(i+1),v(:,i+1),objVec(i+1),exitVec(i+1)] =...
       instantRatesV2(C(1:end-1,i+1), t(i) + h ,Perfusion,cellType, shift, Rtest, parameters, shiftDay); 
end

figure(1);clf;
subplot(2,1,1);plot(t,C(5,:)/2.31); hold on; % Plots biomass
subplot(2,1,1);plot(t,C(end,:)/2.31);        % Plots VCD on same graph as BIO
subplot(2,1,2);plot(t,C(5,:)./C(end,:))      
plotToolV2
reactionScraperPerfusion


if writeFile
C_O2_vec = smooth(C(14,:),15);
C_CO2_vec = smooth(C(7,:),15);
fileName = 'Perfusion_Output.mat';
save(fileName,'t','shiftDay','C_CO2_vec','C_O2_vec','A','Dt');
end
