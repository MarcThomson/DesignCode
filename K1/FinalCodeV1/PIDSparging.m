clear;close all;clc
% Calculates the necesary O2, CO2, and air flow rates for a particular
% process


%% Input
% Dataset:
% load('Batch_Output.mat');load('Batch_Output_PIDParameters.mat');
 load('Perfusion_Output.mat');load('Perfusion_Output_PIDParameters.mat');
% load('Perfusion20.mat');load('Perfusion20_PIDParameters.mat');
% load('Perfusion100.mat');load('Perfusion100_PIDParameters.mat');
% load('Batch20.mat');load('Batch20_PIDParameters.mat'); 
% load('Batch80.mat'); load('Batch80_PIDParameters.mat');
%  load('Batch400.mat');load('Batch400_PIDParameters.mat');

% Other inputs for manual tuning. Keep commented unless tuning
% batchReactor  = 0; %1 for batch reactor, 0 otherwise
%     CO2_0 = 8E-2; %percentage of initial feed that is CO2
%     O2_0 = 45E-2; %percentage of initial feed that is O2
%   ImpVesRatio = 0.33;   % ratio of impeller diameter to vessel diameter
%  MinBubble = 0.69e-2;  % minimum bubble diameter for deadband control, m
%  MaxBubble = 0.7e-2;  % maximum bubble diameter for deadband control, m
% NImpeller = 1;
% 
% % PID Parameters
%  K = [  1     50       0; %O2
%         10     100       0; %CO2   
%           0     0         0];%Flow
% %       K1    K2        K3


%% calculations
% reset the initial concentrations to a reasonable level (otherwise they
% might start negative due to the nature of the date files)
C_O2_vec = C_O2_vec - C_O2_vec(1) + 0.4189; 
C_CO2_vec = C_CO2_vec - C_CO2_vec(1) + 1.7481;


h1 = 10/24/60/10; %mass transfer time step
loopTime = 10/24/60; %PID loop step time 

% initially, the temperature is 37 dC
shift = 0;


if batchReactor
    % in the batch case, define t as with more precision near the feeding
    % days
    t_new = t(1):h1:2.9;
    for i = 3:9
        t_new = [t_new,i-0.1+h1:h1/3:i+0.1];
        t_new = [t_new,i+0.1+h1:h1:i+1-0.1];
    end
    t_new = [t_new,t_new(end)+h1:h1:t(end)];
else
    % otherwise, define t with constant steps
    t_new = t(1):h1:t(end);
end

% interpolate the values from the input. These will represent the O2/CO2
% concentration over time if there were no external source/sink
C_O2_rxn = interp1(t,C_O2_vec,t_new,'linear');
C_CO2_rxn = interp1(t,C_CO2_vec,t_new,'linear');
C_VCD_rxn = interp1(t,C_VCD_vec,t_new,'linear');

% these vectors represent the total O2/CO2 concentrations only from the
% sparging system
C_O2_sparge = 0 * C_O2_rxn;
C_CO2_sparge = 0 * C_CO2_rxn;

% initialize vectors to house the saturation proportion over time
CO2prop_vec = 0 * C_O2_rxn;
O2prop_vec = 0 * C_CO2_rxn;

% initialize some other important vectors
flowTotal = zeros(length(t_new),1); %total flow rate, m^3/s
flowO2 = zeros(length(t_new),1); % flow rate of O2, m^3/s
flowCO2 = zeros(length(t_new),1); %flow ratw of CO2, m^3/s
Nrad = zeros(length(t_new),1); %rotation speed of impellers, rad/s
PV = zeros(length(t_new),1);   % EDR in W/m^3
Qtherm = zeros(length(t_new),1);   % Q (heat from cells + ht from air) in J
T_tank = zeros(length(t_new),1);   % Tank temperature over time
T_jacket = zeros(length(t_new),1);   % Jacket temperature over time

% Select some reasonable initial values 
Nrad(1) = 1*2*pi; %60 rpm
flowTotal(1) = 10e-3*A; % aprroximately the middle allowable flow rate
flowCO2(1) = CO2_0 *flowTotal(1); % 6% of the initial flow is 
flowO2(1) = O2_0*flowTotal(1);
PV(1) = 0;         
Qtherm(1) = 0;

for i = 2 : length(t_new)
    % define h more easily within this loop
    h = t_new(i) - t_new(i-1);
    
    % shift temperature if past the shift day
    if t_new(i)>shiftDay
        shift = 1;
    end
    
   % compute the mole fractions of O2,CO2 in the gas feed
   yO2 = flowO2(i-1)/flowTotal(i-1);
   yCO2 = flowCO2(i-1)/flowTotal(i-1);
   
   % Use lever rule to calculate the proportion of the feed that must be
   % CO2 supplement, and the amount that must be air
   CO2Frac = (yCO2-400E-6)/(1-400E-6);
   airFrac = (1-yO2-CO2Frac)/(1-0.2095);
        
   % Use Runga-Kutta 4 to integrate the mass transfer
   
    %calculate the approximate saturation percentage andd initial flux
    [FO2_1, FCO2_1, ~, O2prop_vec(i-1), CO2prop_vec(i-1), ~,~, Qbub1] = ...
    massTransferModel( C_O2_sparge(i-1) + C_O2_rxn(i-1),...
               C_CO2_sparge(i-1) + C_CO2_rxn(i-1),...
               airFrac,CO2Frac,flowTotal(i-1), Dt, Nrad(i-1), t_new(i-1),...
               shift, MinBubble, MaxBubble, ImpVesRatio, NImpeller);
    FO2_1 = 24*3600*FO2_1; %convert from /s basis to /day
    FCO2_1 = 24*3600*FCO2_1; %convert from /s basis to /day
    Qbub1 = Qbub1 - 38*H*A*C_VCD_rxn(i-1);
    Qbub1 =  24*3600*Qbub1; %convert from /s basis to /day
    
   %complete step 2 of RK4 
   [FO2_2, FCO2_2, ~, ~, ~, ~,~, Qbub2] = ...
    massTransferModel( C_O2_sparge(i-1) + FO2_1*h/2 + 1/2*(C_O2_rxn(i-1)+ C_O2_rxn(i))   ,...
               C_CO2_sparge(i-1) + FCO2_1*h/2 + 1/2*(C_CO2_rxn(i-1) + C_CO2_rxn(i)),...
               airFrac, CO2Frac, flowTotal(i-1), Dt, Nrad(i-1),...
               t_new(i-1) + h/2, shift, MinBubble, MaxBubble, ImpVesRatio, NImpeller);
    FO2_2 = 24*3600*FO2_2; %convert from /s basis to /day
    FCO2_2 = 24*3600*FCO2_2; %convert from /s basis to /day
    Qbub2 = Qbub2 - 38*H*A*(C_VCD_rxn(i-1)+C_VCD_rxn(i))/2;
    Qbub2 =  24*3600*Qbub2; %convert from /s basis to /day
    
    %complete step 3 of RK4 
    [FO2_3, FCO2_3, ~, ~, ~, ~, ~, Qbub3] = ...
    massTransferModel( C_O2_sparge(i-1) + FO2_2*h/2 + 1/2*(C_O2_rxn(i-1) + C_O2_rxn(i)),...
               C_CO2_sparge(i-1)+ FCO2_2*h/2 + 1/2*(C_CO2_rxn(i-1) + C_CO2_rxn(i)),...
                    airFrac, CO2Frac, flowTotal(i-1), Dt, Nrad(i-1),...
                    t_new(i-1) + h/2, shift, MinBubble, MaxBubble, ImpVesRatio, NImpeller);
    FO2_3 = 24*3600*FO2_3; %convert from /s basis to /day
    FCO2_3 = 24*3600*FCO2_3; %convert from /s basis to /day
    Qbub3 = Qbub3 - 38*H*A*(C_VCD_rxn(i-1)+C_VCD_rxn(i))/2;
    Qbub3 =  24*3600*Qbub3; %convert from /s basis to /day
    
    %complete step 4 of RK4
    [FO2_4, FCO2_4, ~, ~, ~, ~, ~, Qbub4] = ...
    massTransferModel(C_O2_sparge(i-1) + FO2_3*h + C_O2_rxn(i),...
              C_CO2_sparge(i-1) + FCO2_3*h + C_CO2_rxn(i)  ,...
                    airFrac, CO2Frac, flowTotal(i-1), Dt, Nrad(i-1),...
                    t_new(i-1) + h, shift, MinBubble, MaxBubble, ImpVesRatio, NImpeller);
    FO2_4 = 24*3600*FO2_4; %convert from /s basis to /day
    FCO2_4 = 24*3600*FCO2_4; %convert from /s basis to /day
    Qbub4 = Qbub4 - 38*H*A*C_VCD_rxn(i);
    Qbub4 =  24*3600*Qbub4; %convert from /s basis to /day
    
    
    % compute the overall derivative from a weighted average
    FO2 = (FO2_1 + 2*FO2_2 + 2*FO2_3 + FO2_4)/6;
    FCO2 = (FCO2_1 + 2*FCO2_2 + 2*FCO2_3 + FCO2_4)/6;            
    Qbub  = (Qbub1 + 2*Qbub2 + 2*Qbub3 + Qbub4)/6;     
    
    
    % update for ith iteration
    C_O2_sparge(i) = C_O2_sparge(i-1) + h * FO2;
    C_CO2_sparge(i) = C_CO2_sparge(i-1) + h * FCO2;
    Qtherm(i) = Qtherm(i-1) + h * Qbub;
    
    % calculate one final time to update the O2prop, CO2prop, impeller
    % speed and bubble size
    [~, ~, Nrad(i), O2prop, CO2prop, dB(i), PV(i),~] = ...
    massTransferModel(C_O2_sparge(i) + C_O2_rxn(i),...
              C_CO2_sparge(i) + C_CO2_rxn(i),...
              airFrac, CO2Frac, flowTotal(i-1), Dt, Nrad(i-1), t_new(i), ...
              shift, MinBubble, MaxBubble, ImpVesRatio, NImpeller);
    
    % add the saturation proportions to the vectors
    CO2prop_vec(i) = CO2prop;
    O2prop_vec(i) = O2prop;
   
    % PID 
    % compute errors
    errorO2 = 0.36 - O2prop_vec(1:i);
    errorCO2 = 0.07 - CO2prop_vec(1:i);
    % after the 2nd iteration, begin implementing PI
    if i > 2
        % compute derivatives with backward differences 
        derO2 = (-3*errorO2(i) - errorO2(i-2) + 4*errorO2(i-1))/(2*h);
        derCO2 = (-3*errorCO2(i) - errorCO2(i-2) + 4*errorCO2(i-1))/(2*h);
        
        % use derivatives and the error function to update the output
        % function
        uO2 = K(1,1)*errorO2(i) + K(1,2) * trapz(t_new(1:i), errorO2(1:i)) + K(1,3) * derO2;
        uCO2 = K(2,1)*errorCO2(i) + K(2,2) * trapz(t_new(1:i), errorCO2(1:i)) + K(2,3) * derCO2;
        uFlow = K(3,1)*errorCO2(i) + K(3,2) * trapz(t_new(1:i), errorCO2(1:i)) + K(3,3) * derCO2;
    else
       uO2 = 0; 
       uCO2 = 0;
       uFlow = 0;
    end
    
    % Only update the parameter if on a loop time
    if mod(t_new(i),loopTime) < h/2       
        
 
        flowO2(i) = flowO2(1)*(1 + uO2);    
        flowCO2(i) = flowCO2(1)*(1 + uCO2);
        flowTotal(i) = flowTotal(1)*(1 - uFlow);
        
        % Reset flow if flow is outside of acceptable bounds
        if flowTotal(i)/A > 18.1e-3
            flowTotal(i) = A*18.1e-3;
        elseif flowTotal(i)/A < 3.7e-3
            flowTotal(i) = A*3.7e-3;
        end
       
        % if O2 flow is higher than the total flow minus CO2 flow, reset it
        if flowO2(i)>flowTotal(i)*0.95 
            % if the O2 flow is too high, adjust the total flow to
            % compensate, but don't let it go beyond acceptable bounds
            flowTotal(i) = flowO2(i)/(1-CO2Frac)/0.95;
            if flowTotal(i)/A > 18.1e-3
                flowTotal(i) = A*18.1e-3;
            end
            flowO2(i) = flowTotal(i)*(1-CO2Frac)*0.95;
            %QCO2(i) = CO2Frac * Qtotal(i);
        end
        % if O2 concentration is lower than 0.2095, reset it
        if flowO2(i) < 0.2095*(flowTotal(i)-flowCO2(i))
            % if the O2 flow is too low, adjust the total flow to
            % compensate, but don't let it go beyond acceptable bounds
            flowO2(i) = 0.2095*flowTotal(i);
            flowTotal(i) = flowO2(i)/(1-CO2Frac);
            if flowTotal(i)/A < 3.7e-3
                flowTotal(i) = A*3.7e-3;
            end
            flowO2(i) = 0.2095*flowTotal(i)*(1-CO2Frac);
            %QCO2(i) = CO2Frac * Qtotal(i);
        end
        
        % if CO2 goes negative, force it up
        if flowCO2(i) < 0
            flowCO2(i) = 0;
        end
        
        % if CO2 gets too big, force it down
        if flowCO2(i) > flowTotal(i)
            flowCO2(i) = flowTotal(i)*0.95;
        end
        
       
        
    % if not on a loop timestep, don't change it    
    else
        flowO2(i) = flowO2(i-1);
        flowCO2(i) = flowCO2(i-1);
        flowTotal(i) = flowTotal(i-1);
    end
    
    
     % plot on each day for PID tuning. Delete later
     if floor(t_new(i))-t_new(i) ==  0
        figure(1);plot(t_new(1:i-1),100 * O2prop_vec(1:i-1),'LineWidth',2);
        xlabel('Time (days)');ylabel('O_2 Saturation Percentage');set(gca,'FontSize',20);
        drawnow
        
        figure(2);plot(t_new(1:i-1),100*CO2prop_vec(1:i-1),'LineWidth',2);
        xlabel('Time (days)');ylabel('CO_2 Saturation Percentage');set(gca,'FontSize',20);
        drawnow
     end
    
    
end

% plot at the end
figure(1);plot(t_new(1:i-1),100 * O2prop_vec(1:i-1),'LineWidth',2);
xlabel('Time (days)');ylabel('O_2 Saturation Percentage');set(gca,'FontSize',20);


figure(2);plot(t_new(1:i-1),100*CO2prop_vec(1:i-1),'LineWidth',2);
xlabel('Time (days)');ylabel('CO_2 Saturation Percentage');set(gca,'FontSize',20);
drawnow
% total flow 
yO2 = flowO2./flowTotal;
yCO2 = flowCO2./flowTotal;
CO2Frac = (yCO2-400E-6)/(1-400E-6);
airFrac = (1-yO2-CO2Frac)/(1-0.2095);

%all in m^3

%rate of flow, m^3/s
CO2_flow = flowTotal.*CO2Frac;
air_flow = flowTotal.*airFrac;
O2_flow = flowTotal.*(1-airFrac-CO2Frac);

% total flows, m^3
CO2_total = trapz(t_new,CO2_flow )*3600*24;
air_total = trapz(t_new,air_flow)*3600*24;
O2_total = trapz(t_new,O2_flow)*3600*24;

visc = (t_new>=shiftDay).*7.842E-7 + 6.969E-7*(t_new<shiftDay); %kinematic viscosity of water m^2/s
visc = visc';
phi = 15; % ratio of max EDR to mean EDR
lambda = (visc.^3./(phi * (PV/1000))).^(1/4)*1E6; % Kolmogorov eddy size micro m

% derivative of Q
Qdot = gradient(Qtherm)./gradient(t_new');

% analyze T in reactor
heatTransferAnalysis


% put the outputs into the format of the excel doc
excelOutput = [batchReactor, MinBubble,MaxBubble,ImpVesRatio,CO2_0, ...
              O2_0,NImpeller, K(1,:),K(2,:),K(3,:), O2_total,CO2_total,...
              air_total,max(O2_flow),max(CO2_flow), max(air_flow),...
              max(Nrad)/(2*pi)*60, min(lambda)]
          
excelOutputHT = [T_beforeShift, T_afterShift,flowJacket,volumeJacket,...
                 minBeforeShift, maxBeforeShift, minAfterShift,maxAfterShift,...
                 totalDuty]
                 
