function Qgas = spargerResetFlow(A,Qgas,minFlow,maxFlow)
% resets the flow rate of gas through the sparger if the superficial
% velocity is outside of acceptible bounds
%inputs:
% A = area of the sparger, m^2
% Qgas = current total gas flow rate, m^3/2
% minFlow is minimum superficial velocity, m/s
% maxFlow is the maximum supercial velocity, m/s
% outputs:
% Qgas is the adjusted gas flow rate, if necesary

    if Qgas/A < minFlow
        Qgas = minFlow*A;
    elseif Qgas/A > maxFlow
        Qgas = maxFlow*A;
    end
end