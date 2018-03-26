function [FO2, FCO2, Nrad, O2prop, CO2prop, dB, PV, Q] = ...
    massTransferModel(CL_O2, CL_CO2, airFrac, CO2Frac, Qgas, Dt, Nrad, t,...
              shift, MinBubble, MaxBubble, ImpVesRatio, NImpeller)
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
% Minbubble = minimum desired bubble size, m
% Maxbubble = maximum desired bubble size, m
% ImpVesRatio = ratio of impeller diameter to vessel diameter
% NImpeller = number of impellers
% outputs:
% FO2 = O2 flux due to mass transfer, mM/s
% FCO2 = CO2 flux due to mass transfer, mM/s
% Nrad = new impeller speed, rad/s
% O2prop = current O2 saturation percentage
% CO2prop = current CO2 saturation percentage
% dB = current mean bubble diameter, m
% PV = EDR in W/m^3
% Q = heat transferred to gas, W

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
if NImpeller == 2
    P0 = P01 + P02; %total power input for nonaerated system, W
elseif NImpeller == 1
    P0 = P01;
end

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
    if NImpeller == 2
        P0 = P01 + P02; %total power input for nonaerated system, W
    elseif NImpeller == 1
        P0 = P01;
    end

    Pg = P0*(1-38.2*Qgas/A/sqrt(g*impDiam)); %actual power input, W
    PV = (Pg)/VR; %EDR, W/m^3

    dB=(PV^-.2)*10.1*10^-3; % Sauter mean diameter
end


kLaO2 = 0.0218*PV^.5*Qgas^.6; % kLa for O2, Hz

if NImpeller == 1
    kLaO2 = 0.7 * kLaO2;
end

kLaCO2 = 0.92*kLaO2; % kLa for CO2, Hz

FO2 = kLaO2*(CstarO2-CL_O2);  %O2 flux, mM/s
FCO2 = kLaCO2*(CstarCO2-CL_CO2);  %CO2 flux, mM/s

eps = 22.4*(PV^0.24)*(Qgas^0.65)/100; % void fraction. CHECK
a = 6*eps/(dB*(1-eps)); % bubble specific area m^2/m^3 
ahead = A/VR; % head specific area m^2/m^3 

kLaO2_head = ahead/a*kLaO2;     %find kLa for headspace, Hz
kLaCO2_head = ahead/a*kLaCO2;   %find kLa for headspace, Hz


% Recompute Henry's law for headspace
PO2 = 1E5*0.2095; %Total Pressure of O2 in air, Pa
PCO2 = 1E5*400e-6;%Total Pressure of CO2 in air, Pa
CstarO2 = PO2/hO2;      % saturation concentation O2 in liquid, mol/m^3
CstarCO2 = PCO2/hCO2;   % saturation concentration CO2 in liquid, mol/m^3

FO2 = FO2 + kLaO2_head*(CstarO2-CL_O2);  % total O2 flux, mM/s
FCO2 = FCO2 + kLaCO2_head*(CstarCO2-CL_CO2);  % total CO2 flux, mM/s


CstarO2 = 1E5/hO2;    % O2 sat conc under pure O2 atmosphere mol/m^3 = mM
CstarCO2 = 1E5/hCO2;  % CO2 sat conc under pure CO2 atmosphere mol/m^3 = mM

CO2prop = CL_CO2/CstarCO2; % CO2 saturation concentration
O2prop = CL_O2/CstarO2;    % O2 saturation concentration

Nrad = Nrev*2*pi; % new impeller speed, rad/s


% Account for gas-liquid heat transfer
Cpgas = 1.005E3; %J/kgK
Pentry = height*g*1000 + 1E5; %Pa
rhogas = 1.2 * Pentry/1E5; %kg/m^3

Q = (T-25)*Cpgas*rhogas*Qgas; %W

end
     