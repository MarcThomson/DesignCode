%initial conditions for CELL LINE A:

c_AKG_0         =    0;    %  mM  =  1
c_ALA_ext_0     = 0.67;    %  mM  =  2
c_ALA_int_0     =    0;    %  mM  =  3
c_ANTI_ext_0    =   0.14;    %  mM  =  4
c_ANTI_int_0    =   0;    %  mM  =  5
c_ASN_ext_0     =  20.27;    %  mM  =  6
c_ASN_int_0     =    0;    %  mM  =  7
c_ASP_ext_0     =   2.05;    %  mM  =  8
c_ASP_int_0     =    0;    %  mM  =  9
c_ATP_0         =    0;    %  mM  = 10
c_BIOM_ext_0    =   4.158;    %  mM  =  11 RELATED TO VCD % 5.4*2.3 is probably what to use
c_BIOM_int_0    = 0;    %  mM  =  12
c_C_C_ext_0     = 1.50;    %  mM  =  13
c_C_C_int_0     =    0;    %  mM  =  14
c_CO2_ext_0     =    1.7481;    %  mM  =  15
c_CO2_int_0     =    0;    %  mM  =  16
c_CYS_0         =    0;    %  mM  =  17
c_FADH2_0       =    0;    %  mM  =  18
c_G6P_0         =    0;    %  mM  =  19
c_GLC_ext_0     =  73;    %  mM  =  20 % GLUCOSE
c_GLN_ext_0     =  3.9;    %  mM  =  21
c_GLN_int_0     =    0;    %  mM  =  22
c_GLU_ext_0     = 0.62;    %  mM  =  23
c_GLU_int_0     =    0;    %  mM  =  24
c_GLY_ext_0     = 3.3;    %  mM  =  25
c_GLY_int_0     =    0;    %  mM  =  26
c_LAC_ext_0     =   5.55;    %  mM  =  27
c_LAC_int_0     =    0;    %  mM  =  28
c_MAL_0         =    0;    %  mM  =  29
c_NADH_cyto_0   =    0;    %  mM  =  30
c_NADH_mito_0   =    0;    %  mM  =  31
c_NH3_ext_0     =   0.68;    %  mM  =  32
c_NH3_int_0     =    0;    %  mM  =  33
c_O2_ext_0      =   0.4189*1.2;    %  mM  =  34
c_O2_int_0      =    0;    %  mM  =  35
c_OXA_0         =    0;    %  mM  =  36
c_PYR_0         =    0;    %  mM  =  37
c_SER_ext_0     =  11.06;    %  mM  =  38
c_SER_int_0     =   0;    %  mM  =  39
shift_0         =   0;    %    
R_0             =   6.2;    %units

initialConditions_vec = [c_AKG_0, c_ALA_ext_0, c_ALA_int_0, c_ANTI_ext_0,...
    c_ANTI_int_0, c_ASN_ext_0, c_ASN_int_0, c_ASP_ext_0, c_ASP_int_0,...
    c_ATP_0, c_BIOM_ext_0, c_BIOM_int_0, c_C_C_ext_0, c_C_C_int_0,...
    c_CO2_ext_0, c_CO2_int_0, c_CYS_0, c_FADH2_0, c_G6P_0, c_GLC_ext_0,...
    c_GLN_ext_0, c_GLN_int_0, c_GLU_ext_0, c_GLU_int_0, c_GLY_ext_0,...
    c_GLY_int_0, c_LAC_ext_0, c_LAC_int_0, c_MAL_0, c_NADH_cyto_0,...
    c_NADH_mito_0, c_NH3_ext_0, c_NH3_int_0, c_O2_ext_0, c_O2_int_0,...
    c_OXA_0, c_PYR_0, c_SER_ext_0, c_SER_int_0, shift_0, R_0]';
