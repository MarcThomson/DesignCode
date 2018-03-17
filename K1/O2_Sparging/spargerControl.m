function [Qgas,QO2,QCO2,boolChanged] = spargerControl(Qgas0,O2Prop,CO2Prop,QO20,QCO20,A,t,tChanged)
%controls the flow rate through the sparger
% Inputs
% Qgas0 = current flow rate of gas into the sparger, m^3/s
% O2Prop = current saturation percentage of the O2
% CO2Prop = current saturation percentage of the CO2
% QO20 = current flow rate of supplemental pure O2 in the feed, m^3/s
% QCO20 = current flow rate of supplemental pure CO2 in the feed, m^3/s
% A = area of the sparger, m^2
% t = current time, day
% tChanged = last time the system changed, days
%
%Outputs:
%Qgas = new total gas flow rate, m^3/s
%QO2 = new pure O2 flow rate, m^3/s
%QCO2 = new pure CO2 flow rate, m^3/s
%boolChanged = whether or not the system has changed on this iteration



%Control Parameters
percentChangeO2 = 0.15; %how much the system changes when a change is requested, for O2
percentChangeCO2 = 0.04; %how much the system changes when a change is requested, for CO2/Flow
minCO2=0.065; %min sat % CO2 for deadband
maxCO2=0.075; %max say % CO2 for deadband
minO2=0.4; %min sat % O2 for deadband
maxO2=0.5; %max sat % O2 for deadband
minFlow = 3.7e-3; %min superficial velocty, m/s
maxFlow = 18.1e-3; %max superficial velocity, m/s
tFeed = 1/24; %how close a time must be to a feeding to add extra O2, decrease CO2 (hours)
dtCrit = 1/24/2; % how often the control system can change the parameters

% if the parameters have changed in the last hours, do not change them
% again
if (t-tChanged < dtCrit) && ~ (t > 2.5 && abs(round(t)-t) < tFeed && round(t)>t && abs(round(t)-tChanged) > tFeed  ) 
    boolChanged=0;
    Qgas = Qgas0;
    QO2 = QO20;
    QCO2 = QCO20;
else
    %if O2 too low or nearing a feeding, increase O2 flow
    %if O2 too high, reduce O2 flow
    if O2Prop < minO2 || (t > 2.5 && abs(round(t)-t) < tFeed && round(t)>t)
        QO2 = (1 + percentChangeO2) * QO20;
    elseif O2Prop > maxO2   
        QO2 = (1 - percentChangeO2) * QO20;
    else %otherwise stay the same
        QO2 = QO20;
    end   

    %if CO2 too high or nearing feed, increase total flow and decrease CO2
    %in feed
    %if CO2 too low, decrease total flor and increase CO2
    if CO2Prop > maxCO2 || (t > 2.5 && abs(t-round(t)) < tFeed && round(t)>t)
        Qgas = (1+percentChangeCO2) * Qgas0;       
        QCO2 = (1-percentChangeCO2) * QCO20;
    elseif CO2Prop < minCO2
        Qgas = (1-percentChangeCO2) * Qgas0;
        QCO2 = (1+percentChangeCO2) * QCO20;
    else %otherwise stay the same
        Qgas = Qgas0;
        QCO2 = QCO20;
    end    
   
    Qgas = spargerResetFlow(A,Qgas,minFlow,maxFlow);
    
    %If the oxygen desired oxygen concentration in the feed is infeasible,
    %reset it. Also, change flow rate slightly, but do not go beyond the bounds
    if QO2/Qgas < .2095
       Qgas = (1 - percentChangeCO2)*Qgas;
       Qgas = spargerResetFlow(A,Qgas,minFlow,maxFlow);
       QO2 = Qgas * .2095; 
    end
   if QO2/Qgas > 1
       Qgas = (1 + percentChangeCO2)*Qgas;
       Qgas = spargerResetFlow(A,Qgas,minFlow,maxFlow);
       QO2 = Qgas; 
   end
    
   %if the desired CO2 concentration in the feed is infeasible, reset it.
   %Also, slightly modify the flow, but do not go beyond the bounds
    if QCO2/Qgas < 400E-6
       Qgas = (1 - percentChangeCO2)*Qgas;
       Qgas = spargerResetFlow(A,Qgas,minFlow,maxFlow);
       QCO2 = Qgas*400E-6; 
    end
   if QCO2/Qgas > (1-QO2/Qgas)
       Qgas = (1 + percentChangeCO2)*Qgas; 
       Qgas = spargerResetFlow(A,Qgas,minFlow,maxFlow);
       QCO2 = Qgas*(1-QO2/Qgas); 
    end
   
   
    % if the parameters have changed record it and prevent another change
    % for some amount of time
    if QO2~= QO20 || Qgas ~= Qgas0 || QCO2~= QCO20
        boolChanged = 1; 
    else
        boolChanged = 0;
    end
    
end

end

