clear;close all;clc
% Calculates the necesary O2, CO2, and air flow rates for a particular
% process


%% Input
% Dataset:
% load('Batch_Output.mat');load('Batch_Output_PIDParameters.mat');
% load('Perfusion_Output.mat');load('Perfusion_Output_PIDParameters.mat');
% load('Perfusion100.mat');load('Perfusion100_PIDParameters.mat');
 load('Batch20.mat');load('Batch20_PIDParameters.mat');
% load('Batch80.mat'); load('Batch80_PIDParameters.mat');
% load('Batch400.mat');load('Batch400_PIDParameters.mat');

% Other inputs for manual tuning. Keep commented unless tuning
% batchReactor  = 0; %1 for batch reactor, 0 otherwise
%    CO2_0 = 8E-2; %percentage of initial feed that is CO2
%    O2_0 = 45E-2; %percentage of initial feed that is O2
% ImpVesRatio = 0.5;   % ratio of impeller diameter to vessel diameter
% MinBubble = 0.6e-2;  % minimum bubble diameter for deadband control, m
 MaxBubble = 0.68e-2;  % maximum bubble diameter for deadband control, m
% 
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

% these vectors represent the total O2/CO2 concentrations only from the
% sparging system
C_O2_sparge = 0 * C_O2_rxn;
C_CO2_sparge = 0 * C_CO2_rxn;

% initialize vectors to house the saturation proportion over time
CO2prop_vec = 0 * C_O2_rxn;
O2prop_vec = 0 * C_CO2_rxn;

% initialize some other important vectors
Qtotal = zeros(length(t_new),1); %total flow rate, m^3/s
QO2 = zeros(length(t_new),1); % flow rate of O2, m^3/s
QCO2 = zeros(length(t_new),1); %flow ratw of CO2, m^3/s
Nrad = zeros(length(t_new),1); %rotation speed of impellers, rad/s

% Select some reasonable initial values 
Nrad(1) = 1*2*pi; %60 rpm
Qtotal(1) = 10e-3*A; % aprroximately the middle allowable flow rate
QCO2(1) = CO2_0 *Qtotal(1); % 6% of the initial flow is 
QO2(1) = O2_0*Qtotal(1);


for i = 2 : length(t_new)
    % define h more easily within this loop
    h = t_new(i) - t_new(i-1);
    
    % shift temperature if past the shift day
    if t_new(i)>shiftDay
        shift = 1;
    end
    
   % compute the mole fractions of O2,CO2 in the gas feed
   yO2 = QO2(i-1)/Qtotal(i-1);
   yCO2 = QCO2(i-1)/Qtotal(i-1);
   
   % Use lever rule to calculate the proportion of the feed that must be
   % CO2 supplement, and the amount that must be air
   CO2Frac = (yCO2-400E-6)/(1-400E-6);
   airFrac = (1-yO2-CO2Frac)/(1-0.2095);
        
   % Use Runga-Kutta 4 to integrate the mass transfer
   
    %calculate the approximate saturation percentage andd initial flux
    [FO2_1, FCO2_1, ~, O2prop_vec(i-1), CO2prop_vec(i-1), ~] = ...
    spargerV2( C_O2_sparge(i-1) + C_O2_rxn(i-1),...
               C_CO2_sparge(i-1) + C_CO2_rxn(i-1),...
               airFrac,CO2Frac,Qtotal(i-1), Dt, Nrad(i-1), t_new(i-1),...
               shift, MinBubble, MaxBubble, ImpVesRatio);
    FO2_1 = 24*3600*FO2_1; %convert from /s basis to /day
    FCO2_1 = 24*3600*FCO2_1; %convert from /s basis to /day
    
   %complete step 2 of RK4 
   [FO2_2, FCO2_2, ~, ~, ~, ~] = ...
    spargerV2( C_O2_sparge(i-1) + FO2_1*h/2 + 1/2*(C_O2_rxn(i-1)+ C_O2_rxn(i))   ,...
               C_CO2_sparge(i-1) + FCO2_1*h/2 + 1/2*(C_CO2_rxn(i-1) + C_CO2_rxn(i)),...
               airFrac, CO2Frac, Qtotal(i-1), Dt, Nrad(i-1),...
               t_new(i-1) + h/2, shift, MinBubble, MaxBubble, ImpVesRatio);
    FO2_2 = 24*3600*FO2_2; %convert from /s basis to /day
    FCO2_2 = 24*3600*FCO2_2; %convert from /s basis to /day
    
    %complete step 3 of RK4 
    [FO2_3, FCO2_3, ~, ~, ~, ~] = ...
    spargerV2( C_O2_sparge(i-1) + FO2_2*h/2 + 1/2*(C_O2_rxn(i-1) + C_O2_rxn(i)),...
               C_CO2_sparge(i-1)+ FCO2_2*h/2 + 1/2*(C_CO2_rxn(i-1) + C_CO2_rxn(i)),...
                    airFrac, CO2Frac, Qtotal(i-1), Dt, Nrad(i-1),...
                    t_new(i-1) + h/2, shift, MinBubble, MaxBubble, ImpVesRatio);
    FO2_3 = 24*3600*FO2_3; %convert from /s basis to /day
    FCO2_3 = 24*3600*FCO2_3; %convert from /s basis to /day
    
    %complete step 4 of RK4
    [FO2_4, FCO2_4, ~, ~, ~, ~] = ...
    spargerV2(C_O2_sparge(i-1) + FO2_3*h + C_O2_rxn(i),...
              C_CO2_sparge(i-1) + FCO2_3*h + C_CO2_rxn(i)  ,...
                    airFrac, CO2Frac, Qtotal(i-1), Dt, Nrad(i-1),...
                    t_new(i-1) + h, shift, MinBubble, MaxBubble, ImpVesRatio);
    FO2_4 = 24*3600*FO2_4; %convert from /s basis to /day
    FCO2_4 = 24*3600*FCO2_4; %convert from /s basis to /day
    
    % compute the overall derivative from a weighted average
    FO2 = (FO2_1 + 2*FO2_2 + 2*FO2_3 + FO2_4)/6;
    FCO2 = (FCO2_1 + 2*FCO2_2 + 2*FCO2_3 + FCO2_4)/6;            
    
    % update for ith iteration
    C_O2_sparge(i) = C_O2_sparge(i-1) + h * FO2;
    C_CO2_sparge(i) = C_CO2_sparge(i-1) + h * FCO2;
    
    % calculate one final time to update the O2prop, CO2prop, impeller
    % speed and bubble size
    [~, ~, Nrad(i), O2prop, CO2prop, dB(i)] = ...
    spargerV2(C_O2_sparge(i) + C_O2_rxn(i),...
              C_CO2_sparge(i) + C_CO2_rxn(i),...
              airFrac, CO2Frac, Qtotal(i-1), Dt, Nrad(i-1), t_new(i), ...
              shift, MinBubble, MaxBubble, ImpVesRatio);
    
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
        
 
        QO2(i) = QO2(1)*(1 + uO2);    
        QCO2(i) = QCO2(1)*(1 + uCO2);
        Qtotal(i) = Qtotal(1)*(1 - uFlow);
        
        % Reset flow if flow is outside of acceptable bounds
        if Qtotal(i)/A > 18.1e-3
            Qtotal(i) = A*18.1e-3;
        elseif Qtotal(i)/A < 3.7e-3
            Qtotal(i) = A*3.7e-3;
        end
       
        % if O2 flow is higher than the total flow minus CO2 flow, reset it
        if QO2(i)>Qtotal(i)*0.95 
            % if the O2 flow is too high, adjust the total flow to
            % compensate, but don't let it go beyond acceptable bounds
            Qtotal(i) = QO2(i)/(1-CO2Frac)/0.95;
            if Qtotal(i)/A > 18.1e-3
                Qtotal(i) = A*18.1e-3;
            end
            QO2(i) = Qtotal(i)*(1-CO2Frac)*0.95;
            %QCO2(i) = CO2Frac * Qtotal(i);
        end
        % if O2 concentration is lower than 0.2095, reset it
        if QO2(i) < 0.2095*(Qtotal(i)-QCO2(i))
            % if the O2 flow is too low, adjust the total flow to
            % compensate, but don't let it go beyond acceptable bounds
            QO2(i) = 0.2095*Qtotal(i);
            Qtotal(i) = QO2(i)/(1-CO2Frac);
            if Qtotal(i)/A < 3.7e-3
                Qtotal(i) = A*3.7e-3;
            end
            QO2(i) = 0.2095*Qtotal(i)*(1-CO2Frac);
            %QCO2(i) = CO2Frac * Qtotal(i);
        end
        
        % if CO2 goes negative, force it up
        if QCO2(i) < 0
            QCO2(i) = 0;
        end
        
        % if CO2 gets too big, force it down
        if QCO2(i) > Qtotal(i)
            QCO2(i) = Qtotal(i)*0.95;
        end
        
       
        
    % if not on a loop timestep, don't change it    
    else
        QO2(i) = QO2(i-1);
        QCO2(i) = QCO2(i-1);
        Qtotal(i) = Qtotal(i-1);
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
yO2 = QO2./Qtotal;
yCO2 = QCO2./Qtotal;
CO2Frac = (yCO2-400E-6)/(1-400E-6);
airFrac = (1-yO2-CO2Frac)/(1-0.2095);

%all in m^3
CO2_total = trapz(t_new,Qtotal.*CO2Frac)*3600*24;
air_total = trapz(t_new,Qtotal.*airFrac)*3600*24;
O2_total = trapz(t_new,Qtotal.*(1-airFrac-CO2Frac))*3600*24;

% put the outputs into the format of the excel doc
%excelOutput = [batchReactor, MinBubble,MaxBubble,ImpVesRatio,CO2_0, O2_0,...
%               K(1,:),K(2,:),K(3,:), O2_total,CO2_total,...
%               air_total,max(Nrad)/(2*pi)*60]

