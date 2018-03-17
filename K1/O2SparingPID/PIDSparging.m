function fobj = PIDSparging(parameters)
%inputs
batchReactor  = 0; %1 for batch reactor, 0 otherwise
%load('Batch_Output.mat');
%load('Perfusion_Output.mat');
 load('Perfusion20.mat');
% load('Perfusion100.mat');
% load('Batch20.mat');
% load('Batch80.mat');
% load('Batch400.mat');

% K = [2 100 0.01 ; %O2
%      2  75    0.01; %CO2   
%      0.3   .2   0.00];%Flow
%   %K1  K2   K3

%%

parameters = reshape(parameters,[3,3]);
K = parameters;

C_O2_vec = C_O2_vec - C_O2_vec(1) + 0.4189*1.2; 
C_CO2_vec = C_CO2_vec - C_CO2_vec(1) + 1.7481;

h1 = 10/24/60/60; %mass transfer
loopTime = 10/24/60; %loop step time

% Dt = 0.95; %m
% A=Dt^2/4*pi;
shift = 0;
%

% only applies to batch!!
if batchReactor
    t_new = t(1):h1:2.9;
    for i = 3:9
        t_new = [t_new,i-0.1+h1:h1/2:i+0.1];
        t_new = [t_new,i+0.1+h1:h1:i+1-0.1];
    end
    t_new = [t_new,t_new(end)+h1:h1:t(end)];
else
    t_new = t(1):h1:t(end);
end

C_O2_rxn = interp1(t,C_O2_vec,t_new,'linear');
C_CO2_rxn = interp1(t,C_CO2_vec,t_new,'linear');

Qtotal=zeros(length(t_new),1);
QO2=zeros(length(t_new),1);
QCO2=zeros(length(t_new),1);
Nrad=zeros(length(t_new),1);

Nrad(1) = 1*2*pi;
Qtotal(1)=10e-3*A;
QCO2(1)=6E-2*Qtotal(1);
QO2(1)=0.4*Qtotal(1);

C_O2_sparge = 0 * C_O2_rxn;
C_CO2_sparge = 0 * C_CO2_rxn;

CO2prop_vec = 0 * C_O2_rxn;
O2prop_vec = 0 * C_CO2_rxn;



for i = 2 : length(t_new)
    h = t_new(i) - t_new(i-1);
    if t_new(i)>shiftDay
        shift = 1;
    end
    
   yO2 = QO2(i-1)/Qtotal(i-1);
   yCO2 = QCO2(i-1)/Qtotal(i-1);
    
   CO2Frac = (yCO2-400E-6)/(1-400E-6);
   airFrac = (1-yO2-CO2Frac)/(1-0.2095);
        
    [FO2_1, FCO2_1, ~, O2prop_vec(i-1), CO2prop_vec(i-1), ~] = ...
    spargerV2( C_O2_sparge(i-1) + C_O2_rxn(i-1),...
               C_CO2_sparge(i-1) + C_CO2_rxn(i-1),...
                    airFrac,CO2Frac,Qtotal(i-1),Dt,Nrad(i-1),t_new(i-1),shift);
    FO2_1 = 24*3600*FO2_1;
    FCO2_1 = 24*3600*FCO2_1;
    
   [FO2_2, FCO2_2, ~, ~, ~, ~] = ...
    spargerV2( C_O2_sparge(i-1) + C_O2_rxn(i-1)*0.5+  C_O2_rxn(i)*0.5   + FO2_1*h/2,...
               C_CO2_sparge(i-1)+ C_CO2_rxn(i-1)*0.5+ C_CO2_rxn(i)*0.5  + FCO2_1*h/2,...
                    airFrac,CO2Frac,Qtotal(i-1),Dt,Nrad(i-1),t_new(i-1)+h/2,shift);
    FO2_2 = 24*3600*FO2_2;
    FCO2_2 = 24*3600*FCO2_2;

    [FO2_3, FCO2_3, ~, ~, ~, ~] = ...
    spargerV2( C_O2_sparge(i-1) + C_O2_rxn(i-1)*  0.5+  C_O2_rxn(i)*0.5 + FO2_2*h/2,...
               C_CO2_sparge(i-1) + C_CO2_rxn(i-1)*0.5 + C_CO2_rxn(i)*0.5+ FCO2_2*h/2,...
                    airFrac,CO2Frac,Qtotal(i-1),Dt,Nrad(i-1),t_new(i-1)+h/2,shift);
    FO2_3 = 24*3600*FO2_3;
    FCO2_3 = 24*3600*FCO2_3;
    
    [FO2_4, FCO2_4, ~, ~, ~, ~] = ...
    spargerV2( C_O2_sparge(i-1) + C_O2_rxn(i) + FO2_3*h,  C_CO2_sparge(i-1) + C_CO2_rxn(i)  + FCO2_3*h,...
                    airFrac,CO2Frac,Qtotal(i-1),Dt,Nrad(i-1),t_new(i-1)+h,shift);
    FO2_4 = 24*3600*FO2_4;
    FCO2_4 = 24*3600*FCO2_4;    
    
    FO2 = (FO2_1 + 2*FO2_2 + 2*FO2_3 + FO2_4)/6;
    FCO2 = (FCO2_1 + 2*FCO2_2 + 2*FCO2_3 + FCO2_4)/6;            
    
     if floor(t_new(i))-t_new(i) ==      0
        figure(1);plot(t_new(1:i-1),100*O2prop_vec(1:i-1),'LineWidth',2);xlabel('time (days');ylabel('O_2 Saturation Percentage');drawnow
        figure(2);plot(t_new(1:i-1),100*CO2prop_vec(1:i-1),'LineWidth',2);;xlabel('time (days');ylabel('CO_2 Saturation Percentage');drawnow
        
     end
    %update for ith iteration
    C_O2_sparge(i) = C_O2_sparge(i-1) + h * FO2;
    C_CO2_sparge(i) = C_CO2_sparge(i-1) + h * FCO2;
    
    %g
    [~, ~, Nrad(i), O2prop, CO2prop, dB(i)] = ...
    spargerV2(C_O2_sparge(i) + C_O2_rxn(i),  C_CO2_sparge(i) + C_CO2_rxn(i)  + FCO2_3*h,...
                    airFrac,CO2Frac,Qtotal(i-1),Dt,Nrad(i-1),t_new(i),shift);
    
    CO2prop_vec(i) = CO2prop;
    O2prop_vec(i) = O2prop;
              
    errorO2 = 0.5 - O2prop_vec(1:i);
    errorCO2 = 0.07 - CO2prop_vec(1:i);
    if i > 2
        derO2 = (-3*errorO2(i) - errorO2(i-2) + 4*errorO2(i-1))/(2*h);
        uO2 = K(1,1)*errorO2(i) + K(1,2) * trapz(t_new(1:i), errorO2(1:i)) + K(1,3) * derO2;
        
        derCO2 = (-3*errorCO2(i) - errorCO2(i-2) + 4*errorCO2(i-1))/(2*h);
        uCO2 = K(2,1)*errorCO2(i) + K(2,2) * trapz(t_new(1:i), errorCO2(1:i)) + K(2,3) * derCO2;
        uFlow = K(3,1)*errorCO2(i) + K(3,2) * trapz(t_new(1:i), errorCO2(1:i)) + K(3,3) * derCO2;
    else
       uO2 = 0; 
       uCO2 = 0;
       uFlow = 0;
    end
    
    if mod(t_new(i),loopTime) < h/2        
        QO2(i) = QO2(1)*(1 + uO2);% QO2_SS +  ;   
        QCO2(i) = QCO2(1)*(1 + uCO2);
        Qtotal(i) = Qtotal(i-1)*(1 + uFlow);
        
        if Qtotal(i)/A>18.1e-3
            Qtotal(i) = A*18.1e-3;
        elseif Qtotal(i)/A<3.7e-3
            Qtotal(i) = A*3.7e-3;
        end
        
        

        if QO2(i)>Qtotal(i)*(1-CO2Frac)
            Qtotal(i) = QO2(i)/(1-CO2Frac);
            if Qtotal(i)/A>18.1e-3
                Qtotal(i) = A*18.1e-3;
            end
            QO2(i) = Qtotal(i)*(1-CO2Frac);
            QCO2(i) = CO2Frac*Qtotal(i);
        end
        if QO2(i) < 0
            QO2(i) = 0;
        end
        if QCO2(i) < 0
            QCO2(i) = 0;
        end
    else
        QO2(i) = QO2(i-1);
        QCO2(i) = QCO2(i-1);
        Qtotal(i) = Qtotal(i-1);
    end
end

n = 4;
ErrorO2 = ((O2prop_vec-0.5)/0.5).^n;
ErrorO2 = trapz(t_new,ErrorO2)^(1/n);
ErrorCO2 = ((CO2prop_vec-0.07)/0.07).^n;
ErrorCO2 = trapz(t_new,ErrorCO2)^(1/n);

fobj = ErrorO2 + ErrorCO2;
        figure(1);plot(t_new(1:i-1),100*O2prop_vec(1:i-1),'LineWidth',2);xlabel('time (days');ylabel('O_2 Saturation Percentage');drawnow
        figure(2);plot(t_new(1:i-1),100*CO2prop_vec(1:i-1),'LineWidth',2);;xlabel('time (days');ylabel('CO_2 Saturation Percentage');drawnow
end