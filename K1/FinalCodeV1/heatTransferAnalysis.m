% Analyzes the tank/jacket temperature for a reactor

U = 350; % heat transfer coefficient (themopedia) W/mK
V = H*A; %Tank volume, m^3
AShell = (A + Dt*pi*H)/1.5612; %area of contact  with the jacket, m^2


deriv_ode = @(t,y)dydt(t,y,t_new,Qdot,shiftDay, flowJacket,volumeJacket, T_beforeShift,T_afterShift, U, AShell, V,H,Dt, batchReactor);
tHTRange = [0 t_new(end)]; 
% Initialize both at 37 dC
y0 = [37;37];

% Calculate the temperatures over time
[tHT, y] = ode23s(deriv_ode,tHTRange,y0);

% plot them
figure(1);
plot(tHT,y,'LineWidth',2);xlabel('Time (days)');ylabel('Temperature (Celsius)');legend('Tank Temperature','Jacket Temperature','Location','NorthEast')
title(['Feed = ',num2str(flowJacket),'kg/s'])


% Find the min/max T before/after the shift
[~,Ibefore] = min(abs(tHT-(shiftDay)));
[~,Iafter] = min(abs(tHT-(shiftDay+1)));
T_tank = y(:,1);
maxBeforeShift = max(T_tank(1:Ibefore));
minBeforeShift = min(T_tank(1:Ibefore));
maxAfterShift =  max(T_tank(Iafter:end));
minAfterShift = min(T_tank(Iafter:end));


% Find the variable Ashell, Volume for the batch reactor
AShellVec = 0*tHT;
VVec = 0*tHT;
for i = 1:length(VVec)
    if batchReactor
    Nfeeds = 0;
        if tHT(i)>=2
            Nfeeds = floor(tHT(i))-2;
            if Nfeeds == 8
                Nfeeds = 7;
            end
        end
        VVec(i) = V*(0.79+0.03*Nfeeds);
        AShellVec(i) = AShell*(0.79+0.03*Nfeeds);
    else
         VVec(i) = V;
         AShellVec(i) = AShell;
    end
end
    
% Find the total rate of heat transfer and total duty
Qtransfer = U*AShellVec.*(y(:,1)-y(:,2))*3600*24; %J/day
totalDuty = trapz(tHT,Qtransfer); % J


% derivative for ODE45
function derivatives = dydt(t,y,t_new,Qdot,shiftDay, flowJacket,volumeJacket, T_beforeShift,T_afterShift, U, AShell, V,H,Dt, batchReactor)
    % deal with the variable Volume in the batch reactor
    if batchReactor
        Nfeeds = 0;
        if t>=2
            Nfeeds = floor(t)-2;
            if Nfeeds == 8
                Nfeeds = 7;
            end
        end
        V = V*(0.79+0.03*Nfeeds);
        H = H*(0.79+0.03*Nfeeds);
        A = pi*Dt^2/4;
        AShell = (A + Dt*pi*H)/2.7622;
    end
    
    % call water properties, scale by the size of the reactor
    Cp = 4184; %J/kg K
    rhoL = 1000; %kg/m^3
    MCpShell = Cp*volumeJacket*rhoL; %J/K
    MCpTank = Cp*V*rhoL;   %J/K
    
    if t > shiftDay
        TFeed = T_afterShift;
    else
        TFeed = T_beforeShift;
    end
    
    % scale to 1/day
    flowJacket = flowJacket*3600*24;
    
    % estimate the cell/air heat transfer by interpolation
    Qtherm_t = interp1(t_new,Qdot',t,'Linear',"extrap");
    
    % compute the jacket/tank heat transfer
    Qtransfer = U*AShell*(y(1)-y(2))*3600*24;
    
    % put them into derivatives
    derivatives = zeros(2,1);
    if t_new(end) == 30 %perfusion derivative
        V_L = 1000*V;
        derivatives(1) = 1/MCpTank * (-Qtherm_t - Qtransfer + V_L*Cp*(37 - y(1))   );
    else
        derivatives(1) = 1/MCpTank * (-Qtherm_t - Qtransfer);
    end
    derivatives(2) = 1/MCpShell*(Qtransfer + flowJacket*Cp*TFeed - flowJacket*Cp*y(2));
end