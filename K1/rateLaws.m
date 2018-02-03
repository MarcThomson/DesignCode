function rates = rateLaws(...
            c_AKG_int, c_ALA_int, c_ANTI_int, c_ANTI_ext, c_ASN_int,...
            c_ASN_ext, c_ASP_int, c_ATP_int, c_BIOM_int, c_BIOM_ext,...
            c_C_C_int, c_C_C_ext, c_CO2_int, c_CYS_int, c_FADH2_int,...
            c_G6P_int, c_GLC_ext, c_GLN_int, c_GLN_ext, c_GLU_int,...
            c_GLY_int, c_LAC_int, c_MAL_int, c_NADH_cyto, c_NADH_mito,...
            c_NH3_int, c_O2_int, c_O2_ext, c_OXA_int, c_PYR_int, ...
            c_SER_int, c_SER_ext)
% 1
V1 = ; %nmol/(10^6 cellls)/day
% 2
V2 = ; %nmol/(10^6 cellls)/day
% 3
V3 = ; %nmol/(10^6 cellls)/day
% 8
V8 = ; %nmol/(10^6 cellls)/day
% 9
V9 = ; %nmol/(10^6 cellls)/day
% 10
V10 = ; %nmol/(10^6 cellls)/day
% 11
V11 = ; %nmol/(10^6 cellls)/day
% 12
V12 = ; %nmol/(10^6 cellls)/day
% 13
V13 = ; %nmol/(10^6 cellls)/day
% 16
V16 = ; %nmol/(10^6 cellls)/day
% 17
V17 = ; %nmol/(10^6 cellls)/day
% 33
V33 = ; %nmol/(10^6 cellls)/day

rates = [
    V1 
    V2 
    V3
    0
    0
    0
    0
    V8
    V9
    V10
    V11
    V12
    V13
    0
    0
    V16
    V17
    0
    0
    0
    0
    0
    0
    0
    0
    0
    0
    0
    0
    0
    0
    0
    V33
    0
    ];