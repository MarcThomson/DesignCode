function rates = rateLaws(...
    c_AKG, c_ALA_ext, c_ALA_int, c_ANTI_ext, c_ANTI_int, c_ASN_ext, ...
    c_ASN_int, c_ASP_ext, c_ASP_int, c_ATP, c_BIOM_ext, c_BIOM_int, ...
    c_C_C_ext, c_C_C_int, c_CO2_ext, c_CO2_int, c_CYS, c_FADH2, c_G6P, ...
    c_GLC_ext, c_GLN_ext, c_GLN_int, c_GLU_ext, c_GLU_int, c_GLY_ext, ...
    c_GLY_int, c_LAC_ext, c_LAC_int, c_MAL, c_NADH_cyto, c_NADH_mito, ...
    c_NH3_ext, c_NH3_int, c_O2_ext, c_O2_int, c_OXA, c_PYR, c_SER_ext, ...
    c_SER_int, shift, R, cellType,parameters, t, shiftDay)

% Outputs the kinetic rates as a function of the current system
% Rate laws provided by: Lee et al. 2012
% inputs:
% C_i = concentration of i in mmol
% shift = 1 if the shift has occured, 0 otherwise
% R = current redox parameter
% cellType = cell line. 1 = A (standard), 0 = B
% parameters = 47x1 vector of exponentiated kinetic parameters
% t =current system time, days
% shiftDay = day the shift occurs on 
% outputs:
% rates = a 34x1 vector of kinetic reaction rates in nmol/(10^6 cells)/day

%define the parameters
parameterDefinitions

% if the shift has occured, use a weighted average to get the parameters
% this assumes that the shift takes 3 hours. After that, preMultiplier will
% be essentiall 0;
if shift==1
    preMultiplier = exp(-33.9400974951428*(t-shiftDay));
    postMultiplier = 1- exp(-33.9400974951428*(t-shiftDay));
else
    preMultiplier = 1;
    postMultiplier = 0;
end


% take the weighted average of parameters
TC_1 = 1*preMultiplier + TC_1b*postMultiplier;
TC_3 = 1*preMultiplier + TC_3b*postMultiplier;
TC_8 = 1*preMultiplier + TC_8b*postMultiplier;
TC_10 = 1*preMultiplier + TC_10b*postMultiplier;
TC_16 = 1*preMultiplier + TC_16b*postMultiplier;
exp_17 = exp_17a*preMultiplier + exp_17b*postMultiplier;
exp_1 = exp_1a*preMultiplier + exp_1b*postMultiplier;

% value of R that changes the redox state of the cell. Lee uses 1, though
% we experimented with other values
Rcrit = 1;

% RXN 1 %nmol/(10^6 cells)/day
V1 = v_max1f/TC_1*(c_GLC_ext/Km_1)/ (1+(c_LAC_ext/Ki_1)^exp_1)/...
     (1+c_GLC_ext/Km_1); 
 
% RXN 2 % nmol/(10^6 cells)/day
if  R >= Rcrit
    V2 = (R-Rcrit)*v_max2;  
else
    V2 = (R-Rcrit)*v_max2*(c_LAC_ext/Km_2)/(1+c_LAC_ext/Km_2);
end

% RXN 3 %nmol/(10^6 cells)/day
if R > Rcrit
    V3 = (v_max3f*(c_GLC_ext/Km_3c)-v_max3r*c_ALA_ext/Km_3b)...
     /(1+(c_GLC_ext/Km_3c)+(c_ALA_ext/Km_3b));
else 
    V3 = (v_max3f*(c_LAC_ext/Km_3a)-v_max3r*(c_ALA_ext/Km_3b))...
    /(1+(c_LAC_ext/Km_3a)+(c_ALA_ext/Km_3b)); 
end
if V3<0
     V3 = -v_max3r/TC_3*(c_ALA_ext/Km_3b)...
   /(1+(c_ALA_ext/Km_3b));    
end


% RXN 8 %nmol/(10^6 cells)/day
V8 = (v_max8f/TC_8*(c_GLN_ext/Km_8a)-v_max8r*(c_GLU_ext/Km_8b)*(c_NH3_ext/Km_8c))... 
/(1+(c_GLN_ext/Km_8a)+(c_GLU_ext/Km_8b)+(c_NH3_ext/Km_8c)+(c_GLU_ext/Km_8b)*(c_NH3_ext/Km_8c));

% RXN 10 %nmol/(10^6 cells)/day
V10 = (v_max10f/TC_10*(c_ASN_ext/Km_10a)-v_max10r/TC_10*(c_ASP_ext/Km_10b)*(c_NH3_ext/Km_10c))...
/(1+(c_ASN_ext/Km_10a)+(c_ASP_ext/Km_10b)+(c_NH3_ext/Km_10c)+(c_ASP_ext/Km_10b)*(c_NH3_ext/Km_10c)); 

% RXN 11 %nmol/(10^6 cells)/day
V11 = v_max11*V10; 

% RXN 9 %nmol/(10^6 cells)/day
V9 = (v_max9f*V3*(c_NH3_ext/Km_9)-v_max9r*V11)...
/(1+(c_NH3_ext/Km_9)); 


% RXN 13 %nmol/(10^6 cells)/day

r = R*preMultiplier + Rcrit*postMultiplier;

V13 = r*v_max13*(c_C_C_ext)/Km_13...
/(1+(c_C_C_ext)/Km_13);               

% RXN 16 %nmol/(10^6 cells)/day
V16_pre = v_max16/TC_16*(c_GLN_ext/Km_16a)*(c_ASN_ext/Km_16b)...
/(1+(c_GLN_ext/Km_16a)+(c_ASN_ext/Km_16b)+(c_GLN_ext/Km_16a)*(c_ASN_ext/Km_16b));

V16_post = v_max16/TC_16*(c_ASN_ext/Km_16b)*(c_GLC_ext/Km_16c)...
/(1+(c_ASN_ext/Km_16b)+(c_GLC_ext/Km_16c)+(c_ASN_ext/Km_16b)*(c_GLC_ext/Km_16c));

V16 = preMultiplier*V16_pre + postMultiplier*V16_post;


% RXN 12 %nmol/(10^6 cells)/day
V12 = (v_max12f*V16*(c_SER_ext/Km_12a)-v_max12r*(c_GLY_ext/Km_12b)^2)...
/(1+(c_SER_ext/Km_12a)+(c_GLY_ext/Km_12b)+(c_GLY_ext/Km_12b)^2); 

% RXN 17 %nmol/(10^6 cells)/day
V17 = v_max17/(1+(c_LAC_ext/Ki_17)^exp_17); 

% RXN 33 %nmol/(10^6 cells)/day
V33 = v_max33a*V8 + v_max33b*V11;

% combine the rates together
rates = [
    V1; V2; V3; 0; 0; 0; 0; V8; V9; V10; V11;
    V12; V13; 0; 0; V16; V17; 0; 0; 0; 0; 0;
    0; 0; 0; 0; 0; 0; 0; 0; 0; 0; V33; 0
    ];
