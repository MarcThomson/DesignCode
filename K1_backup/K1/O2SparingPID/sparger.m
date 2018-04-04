function [FO2,FCO2]=sparger(CL_O2,CL_CO2,Qgas0,A,VR)



% Parameter Initialization
T=37;
g = 9.81;
height = VR/A;
Dt=2*sqrt(A/pi);
Pbar = height*g*1000/2 + 1E5; %Pa
R = 8.314;                    % ideal gas constant
rhoL = 1000;                  % kg/m^3
rhoG = 1.225;                 % kg/m^3
sigmaL = 6.96E-2;             % N/m
N  = 0.01*2*pi;               % revolve rate
muL =  6.9635e-04;            % viscosity of W
muG = 18.03E-6;               % viscosity of air
Tr=.9*Dt;
w=pi*N*Tr;
Vs=Qgas0/A;                   %superficial gas vel
L=1/6*height;                 % Length of mixer shaft
Re=rhoL*N*Dt^2/muL;
difG=2.9E-5/100^2;             % bad value


%  Spline interp
y1=[2.62,2.33,2.44];
x1=[60900,73600,84000];
if Re<max(x1) && Re>min(x1)
eps=interp1(x1,y1,Re)*N^3*Dt^2*rhoL;
else
   fprintf('Help me') 
   return 
end

lengthScale = 2*(sigmaL/(.4*rhoL))^(3/5)*(L^(2/5)/w^(6/5))*(rhoL/rhoG)^0.1;
syms phix
phi=double(solve(phix/(1-phix) == 0.5*(Vs^2/(g*lengthScale))^(1/3)*rhoL/(rhoL-rhoG)));

dbS = 0.7*sigmaL^0.6/(eps^.4*rhoL^.2)*(muL/muG)^.1
a_int = 6*phi/dbS;


hO2 = 1/(1.385-2.635*10^-2*(T-20)+4.288*10^-4*(T-20)^2)*1.013*10^5; % J/mol;
hCO2 = 101330*exp(11.549-2440.4/(T+273)) * 1E-3;  % J/mol


kLaO2 =1.41E-3*(muL/(rhoL*difG))^0.5*(Tr^2*N/(2*pi)*rhoL/muL)^0.67*...
    (rhoL*(N/(2*pi))^2*Dt^3/sigmaL)^1.29; %1.41E-3*(mu/rhoL/DL)^0.5*
kLaCO2 =.91*kLaO2  ;

alphaO2 = VR/Qgas0*R*(Tr+273)/hO2 * kLaO2;
alphaCO2 = VR/Qgas0*R*(Tr+273)/hCO2 * kLaCO2;

CG_CO2_0 = Pbar/(R*(Tr+273)) * 400/1E6; %mol/m^3
CG_O2_0 = Pbar/(R*(Tr+273)) * 0.2095;  %mol/m^3

CG_CO2 = CG_CO2_0/CL_CO2 * hCO2/(R*(Tr+273)); %%%%%%%%%%%%%%%%%%%%%%%%
CG_O2 = CG_O2_0/ hO2/(R*(Tr+273)); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


cSatO2 = R*(273+Tr)*CG_O2/hO2;
cSatCO2 = R*(273+Tr)*CG_CO2/hCO2;


FO2 = (CG_O2_0 - CG_O2) *Qgas0;
FCO2 = (CG_CO2_0 - CG_CO2) *Qgas0;



end

