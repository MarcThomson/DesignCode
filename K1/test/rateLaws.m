function rates = rateLaws(...
    c_AKG, c_ALA_ext, c_ALA_int, c_ANTI_ext, c_ANTI_int, c_ASN_ext, ...
    c_ASN_int, c_ASP_ext, c_ASP_int, c_ATP, c_BIOM_ext, c_BIOM_int, ...
    c_C_C_ext, c_C_C_int, c_CO2_ext, c_CO2_int, c_CYS, c_FADH2, c_G6P, ...
    c_GLC_ext, c_GLN_ext, c_GLN_int, c_GLU_ext, c_GLU_int, c_GLY_ext, ...
    c_GLY_int, c_LAC_ext, c_LAC_int, c_MAL, c_NADH_cyto, c_NADH_mito, ...
    c_NH3_ext, c_NH3_int, c_O2_ext, c_O2_int, c_OXA, c_PYR, c_SER_ext, ...
    c_SER_int, shift, R, cellType)
parameterDefinitions;

% Edward DeCrescenzo
% Marc Thomson
% William Cordell
% Rate laws provided by: Lee et al. 2012


% 1 %nmol/(10^6 cells)/day
V1 = v_max1f/TC_1b*(c_GLC_ext/Km_1)/ (1+(c_LAC_ext/Ki_1)^exp_1)/...
     (1+c_GLC_ext/Km_1); 
% CHECKED BY:
% DATE:

% 2 % nmol/(10^6 cells)/day
if  R>=1
    V2 = (R-1)*v_max2;  
else
    V2 = (R-1)*v_max2*(c_LAC_ext/Km_2)/(1+c_LAC_ext/Km_2);
end
% CHECKED BY:
% DATE:

% 3 %nmol/(10^6 cells)/day
if R>=1 
    V3 = (v_max3f/TC_3b*(c_GLC_ext/Km_3c)-v_max3r*TC_3b*c_ALA_ext/Km_3b)...
    /(1+(c_GLC_ext/Km_3c)+(c_ALA_ext/Km_3b));
else 
    V3 = (v_max3f*(c_LAC_ext/Km_3a)-v_max3r*(c_ALA_ext/Km_3b))...
    /(1+(c_LAC_ext/Km_3a)+(c_ALA_ext/Km_3b)); 
end
% CHECKED BY:
% DATE:

% 8 %nmol/(10^6 cells)/day
V8 = (v_max8f/TC_8b*(c_GLN_ext/Km_8a)-v_max8r*(c_GLU_ext/Km_8b)*(c_NH3_ext/Km_8c))... 
/(1+(c_GLN_ext/Km_8a)+(c_GLU_ext/Km_8b)+(c_NH3_ext/Km_8c)+(c_GLU_ext/Km_8b)*(c_NH3_ext/Km_8c));
% CHECKED BY:
% DATE:

% 10 %nmol/(10^6 cells)/day
V10 = (v_max10f/TC_10b*(c_ASN_ext/Km_10a)-v_max10r/TC_10b*(c_ASP_ext/Km_10b)*(c_NH3_ext/Km_10c))...
/(1+(c_ASN_ext/Km_10a)+(c_ASP_ext/Km_10b)+(c_NH3_ext/Km_10c)+(c_ASP_ext/Km_10b)*(c_NH3_ext/Km_10c)); 
% CHECKED BY:
% DATE:

% 11 %nmol/(10^6 cells)/day
V11 = v_max11*V10; 
% CHECKED BY:
% DATE:

% 9 %nmol/(10^6 cells)/day
V9 = (v_max9f*V3*(c_NH3_ext/Km_9)-v_max9r*V11)...
/(1+(c_NH3_ext/Km_9)); 
% CHECKED BY:
% DATE:

% 13 %nmol/(10^6 cells)/day
if shift==0                     % MUST DEFINE SHIFT BOOLEAN
    r = R;                % This could probably be initialized as R and this simplified
elseif shift ==1 
    r = 1;
else
    'error(Non-boolean shift)'
end
V13 = r*v_max13*(c_C_C_ext)/Km_13...
/(1+(c_C_C_ext)/Km_13);               % WHAT IS GOING ON HERE
% CHECKED BY:
% DATE:

% 16 %nmol/(10^6 cells)/day
if shift == 0 
    V16 = v_max16/TC_16b*(c_GLN_ext/Km_16a)*(c_ASN_ext/Km_16b)...
    /(1+(c_GLN_ext/Km_16a)+(c_ASN_ext/Km_16b)+(c_GLN_ext/Km_16a)*(c_ASN_ext/Km_16b));
elseif shift==1
    V16 = v_max16/TC_16b*(c_ASN_ext/Km_16b)*(c_GLC_ext/Km_16c)...
    /(1+(c_ASN_ext/Km_16b)+(c_GLC_ext/Km_16c)+(c_ASN_ext/Km_16b)*(c_GLC_ext/Km_16c));
else
    'error(Non-boolean shift)'
end
% CHECKED BY:
% DATE:

% 12 %nmol/(10^6 cells)/day
V12 = (v_max12f*V16*(c_SER_ext/Km_12a)-v_max12r*(c_GLY_ext/Km_12b)^2)...
/(1+(c_SER_ext/Km_12a)+(c_GLY_ext/Km_12b)+(c_GLY_ext/Km_12b)^2); 
% CHECKED BY:
% DATE:


% 17 %nmol/(10^6 cells)/day
V17 = v_max17/(1+(c_LAC_ext/Ki_17)^exp_17); 
% CHECKED BY:
% DATE:

% 33 %nmol/(10^6 cells)/day
V33 = v_max33a*V8+v_max33b*V11;
% CHECKED BY:
% DATE:

rates = [
    V1; V2; V3; 0; 0; 0; 0; V8; V9; V10; V11;
    V12; V13; 0; 0; V16; V17; 0; 0; 0; 0; 0;
    0; 0; 0; 0; 0; 0; 0; 0; 0; 0; V33; 0
    ];
