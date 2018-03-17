function [FO2, FCO2, Nrad, O2prop, CO2prop, dB] = ...
    spargerV2(CL_O2,CL_CO2,airFrac,CO2Frac,Qgas,Dt,Nrad,t,shift)
% give the oxygen flux for the current sparger system
% inputs:
% CL_O2 = concentration of O2 in liquid phase, mM
% CL_CO2 = concentration of CO2 in liquid phase, mM
% airFrac = proportion of feed coming from air
% CO2Frac = proportion of feed coming from pure CO2
% Qgas = inlet flow rate of gas, m^3/s
% Dt = diameter of tank, m
% Nrad = initial speed of impellers, rad/s
% t = current time, day
% outputs:
% FO2 = O2 flux due to mass transfer, mM/s
% FCO2 = CO2 flux due to mass transfer, mM/s
% Nrad = new impeller speed, rad/s
% O2prop = current O2 saturation percentage
% CO2prop = current CO2 saturation percentage
% dB = current mean bubble diameter, m

% other inputs
ImpVesRatio = 3/9;   % ratio of impeller diameter to vessel diameter
MinBubble = .5e-2;  % minimum bubble diameter for deadband control, m
MaxBubble = .6e-2;  % maximum bubble diameter for deadband control, m
if shift 
    T = 31; % temperature, celsius
else
    T = 37; % temperature, celsius
end            

% vessel geometry
hDRat = 3;          % height to diameter ratio
height= hDRat*Dt;   % m
A = pi*Dt^2/4;      % m^2
VR = A*height;      % m^3
impDiam = ImpVesRatio*Dt;     % diameter of impeller fin, m

% misc properties
Nrev = Nrad/(2*pi);             % rotation in Hz
g = 9.81;                       % gravity, m^2/s 
Pmean = height*g*1000/2 + 1E5;  % mean vessel pressure, Pa
rhoL = 1000;                    % fluid density, kg/m^3


% Henry's law
PO2 = Pmean*(0.2095*airFrac+(1-airFrac-CO2Frac)); %Total Pressure of O2 in feed, Pa
PCO2 = Pmean*(400e-6 * airFrac +   CO2Frac);      %Total Pressure of CO2 in feed, Pa
hO2 = 1/(1.385-2.635*10^-2*(T-20)+4.288*10^-4*(T-20)^2)*1.013*10^5; % Henry's law coefficient, J/mol            SOURCE??
hCO2 = 101330*exp(11.549-2440.4/(T+273)) * 1E-3;  % Henry's law coefficient, J/mol          SOURCE??
CstarO2 = PO2/hO2;      % saturation concentation O2 in liquid, mol/m^3
CstarCO2 = PCO2/hCO2;   % saturation concentration CO2 in liquid, mol/m^3

% Impeller power
NP1 = 0.75; %impeller numbers, assuming A315
NP2 = 0.75; %impeller numbers, assuming A315
P01 = NP1*rhoL*Nrev^3*(impDiam)^5; %power input for nonaerated system, W
P02 = NP2*rhoL*Nrev^3*(impDiam)^5; %power input for nonaerated system, W
P0 = P01 + P02; %total power input for nonaerated system, W

Pg = P0*(1-38.2*Qgas/A/sqrt(g*impDiam)); %actual power input, W
PV = (Pg)/VR; %EDR, W/m^3

dB=(PV^-.2)*10.1*10^-3; % Sauter mean diameter

% perform deadband control on bubble size
while dB<MinBubble || dB>MaxBubble
    if dB > MaxBubble
        Nrev = 1.01*Nrev;
    elseif dB<MinBubble
        Nrev = 0.99*Nrev;
    end
    
    P01 = NP1*rhoL*Nrev^3*(impDiam)^5; %power input for nonaerated system, W
    P02 = NP2*rhoL*Nrev^3*(impDiam)^5; %power input for nonaerated system, W
    P0 = P01 + P02; %total power input for nonaerated system, W

    Pg = P0*(1-38.2*Qgas/A/sqrt(g*impDiam)); %actual power input, W
    PV = (Pg)/VR; %EDR, W/m^3

    dB=(PV^-.2)*10.1*10^-3; % Sauter mean diameter
end


kLaO2 = 0.0218*PV^.5*Qgas^.6; % kLa for O2, Hz
kLaCO2 = 0.92*kLaO2; % kLa for CO2, Hz

FO2 = kLaO2*(CstarO2-CL_O2);  %O2 flux, mM/s
FCO2 = kLaCO2*(CstarCO2-CL_CO2);  %CO2 flux, mM/s


CstarO2 = 1E5/hO2;    % O2 sat conc under pure O2 atmosphere mol/m^3 = mM
CstarCO2 = 1E5/hCO2;  % CO2 sat conc under pure CO2 atmosphere mol/m^3 = mM

CO2prop = CL_CO2/CstarCO2; % CO2 saturation concentration
O2prop = CL_O2/CstarO2;    % O2 saturation concentration

Nrad = Nrev*2*pi; % new impeller speed, rad/s
end
