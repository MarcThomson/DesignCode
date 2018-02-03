function rates = rateLaws(...
    c_AKG, c_ALA_ext, c_ALA_int, c_ANTI_ext, c_ANTI_int, c_ASN_ext, ...
    c_ASN_int, c_ASP_ext, c_ASP_int, c_ATP, c_BIOM_ext, c_BIOM_int, ...
    c_C_C_ext, c_C_C_int, c_CO2_ext, c_CO2_int, c_CYS, c_FADH2, c_G6P, ...
    c_GLC_ext, c_GLN_ext, c_GLN_int, c_GLU_ext, c_GLU_int, c_GLY_ext, ...
    c_GLY_int, c_LAC_ext, c_LAC_int, c_MAL, c_NADH_cyto, c_NADH_mito, ...
    c_NH3_ext, c_NH3_int, c_O2_ext, c_O2_int, c_OXA, c_PYR, c_SER_ext, ...
    c_SER_int)
parameterDefinitions

% Edward DeCrescenzo
% Marc Thomson
% William Cordell
% Rate laws provided by: Lee et al.


% 1 %nmol/(10^6 cells)/day
V1 = v_max1f/TC1b*c_GLC; 
% CHECKED BY:
% DATE:

% 2 % nmol/(10^6 cells)/day
if  R>=1
V2 = (R-1)v_max2;  
else
V2=(R-1)*v_max2*(c_LAC/K_m1)/(1+c_LAC/K_m2);
end
% CHECKED BY:
% DATE:

% 3 %nmol/(10^6 cells)/day
if R>=1 
V3=(v_max3f/TC_3b*(c_GLC/Km_3c)-v_max3r*TC_3b*c_ALA/K_m3b))...
/(1+(c_GLC/K_m3c)+(c_ALA)/K_m3b));
else 
V3=(v_max3f*(c_LAC/K_m3a)-v_max3r*(c_ALA/K_m3b))...
/(1+(c_LAC/K_m3a)+(c_ALA/K_m3b)); 
end
% CHECKED BY:
% DATE:

% 8 %nmol/(10^6 cells)/day
V8 = (v_max8f/TC_8*(c_GLN/K_m8a)-v_max8r*(c_GLU/K_m8b)*(c_NH3/K_m8c))... 
/(1+(c_GLN/K_m8a)+(c_GLU/K_m8b)+(c_NH3/K_m8c)+(c_GLU/K_m8b)*(c_NH3/K_m8c);
% CHECKED BY:
% DATE:

% 9 %nmol/(10^6 cells)/day
V9 = (v_max9f*v_3(c_NH3/K_m9)-v_max9r*v_11)...
/(1+(c_NH3/K_m9)); 
% CHECKED BY:
% DATE:

% 10 %nmol/(10^6 cells)/day
V10 = (v_max10f/TC_10*(c_ASN/k_m10a)-v_max10r/TC_10*(c_ASP/K_m10b)*(c_NH3/K_m10c)...
/(1+(c_ASN/K_m10a)+(c_ASP/K_m10b)+(c_NH3/K_m10c)+(c_ASP/K_m10b)*(c_NH3/K_m10c)); 
% CHECKED BY:
% DATE:

% 11 %nmol/(10^6 cells)/day
V11 = v_max11*v_10; 
% CHECKED BY:
% DATE:

% 12 %nmol/(10^6 cells)/day
V12 = (v_max12f*v_16*(c_SER/K_m12a)-v_max12r*(c_GLY/K_m12b)^2)...
/(1+(c_SER/K_m12a)+(c_GLY/K_m12b)+(c_GLY/k_m12b)^2; 
% CHECKED BY:
% DATE:

% 13 %nmol/(10^6 cells)/day
if shift==0                     % MUST DEFINE SHIFT BOOLEAN
r=R                % This could probably be initialized as R and this simplified
elseif shift ==1 
r=1
else
error(Non-boolean shift)
end
V13 = r*v_max13*(C-C)/K_m13...
/(1+(C-C)/K_m13);               % WHAT IS GOING ON HERE
% CHECKED BY:
% DATE:

% 16 %nmol/(10^6 cells)/day
if shift == 0 
V16 = v_max16/TC_16*(c_GLN/K_m16a)*(c_ASN/K_m16b)...
/(1+(c_GLN/K_m16)+(c_ASN/K_m16b)+(c_GLN/K_m16a)*(c_ASN/K_m16b);
elseif shift==1
V16 = v_max16/TC_16*(c_ASN/K_m16b)*(c_GLC/K_m16c)...
/(1+(c_ASN/K_m16b)+(c_GLC/K_m16c)+(c_ASN/K_m16b)*(c_GLC/K_m16c);
else
error(Non-boolean shift)
end
% CHECKED BY:
% DATE:

% 17 %nmol/(10^6 cells)/day
V17 = v_max17/(1+(c_LAC/K_i17)^exp_17); 
% CHECKED BY:
% DATE:

% 33 %nmol/(10^6 cells)/day
V33 = TC_33*v_max33*(c_GLN/K_m33a)...
/(1+(c_GLN/K_m33a)+0.01*v_16-0.05*v_17; 
% CHECKED BY:
% DATE:

rates = [
    V1; V2; V3; 0; 0; 0; 0; V8; V9; V10; V11;
    V12; V13; 0; 0; V16; V17; 0; 0; 0; 0; 0;
    0; 0; 0; 0; 0; 0; 0; 0; 0; 0; V33; 0
    ];
