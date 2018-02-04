%initial conditions for CELL LINE A:

c_AKG_0         =    0;    %  mM
c_ALA_ext_0     = 0.9;    %  mM
c_ALA_int_0     =    0;    %  mM
c_ANTI_ext_0    =   0;    %  mM
c_ANTI_int_0    =   0;    %  mM
c_ASN_ext_0     =  21;    %  mM
c_ASN_int_0     =    0;    %  mM
c_ASP_ext_0     =   2;    %  mM
c_ASP_int_0     =    0;    %  mM
c_ATP_0         =    0;    %  mM
c_BIOM_ext_0    =   5;    %  mM RELATED TO VCD % 5.4*2.3 is probably what to use
c_BIOM_int_0    = 0;    %  mM
c_C_C_ext_0     = 1.5;    %  mM
c_C_C_int_0     =    0;    %  mM
c_CO2_ext_0     =    0;    %  mM
c_CO2_int_0     =    0;    %  mM
c_CYS_0         =    0;    %  mM
c_FADH2_0       =    0;    %  mM
c_G6P_0         =    0;    %  mM
c_GLC_ext_0     =  73;    %  mM % GLUCOSE
c_GLN_ext_0     = 3.9;    %  mM
c_GLN_int_0     =    0;    %  mM
c_GLU_ext_0     = 0.6;    %  mM
c_GLU_int_0     =    0;    %  mM
c_GLY_ext_0     = 3.4;    %  mM
c_GLY_int_0     =    0;    %  mM
c_LAC_ext_0     =   6;    %  mM
c_LAC_int_0     =    0;    %  mM
c_MAL_0         =    0;    %  mM
c_NADH_cyto_0   =    0;    %  mM
c_NADH_mito_0   =    0;    %  mM
c_NH3_ext_0     = 0.7;    %  mM
c_NH3_int_0     =    0;    %  mM
c_O2_ext_0      =    0.1666;    %  mM
c_O2_int_0      =    0;    %  mM
c_OXA_0         =    0;    %  mM
c_PYR_0         =    0;    %  mM
c_SER_ext_0     =  11;    %  mM
c_SER_int_0     =   0;    %  mM
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
