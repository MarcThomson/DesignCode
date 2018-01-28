% Initial Conditions 
c_ACCOA_0   = 1.40E-07; % mmol/(10^6 cell) 
c_ADP_0     = 6.70E-07; % mmol/(10^6 cell) %typo assumed
c_AKG_0     = 1.00E-07; % mmol/(10^6 cell) %typo assumed
c_ALA_0     = 9.60E-01; % mM
c_AMP_0     = 3.10E-07; % mmol/(10^6 cell) 
c_ARG_0     = 2.4;      % mM
c_ASN_0     = 4.80E-01; % mM
c_ASP_0     = 2;        % mM
c_ATP_0     = 7.00E-06; % mmol/(10^6 cell) 
c_CIT_0     = 1.00E-06; % mmol/(10^6 cell) 
c_CYS_0     = 1.00E-01; % mM
c_EGLC_0    = 26;       % mM
c_EGLN_0    = 3.7;      % mM
c_EGLU_0    = 8.00E-01; % mM
c_F6P_0     = 2.20E-07; % mmol/(10^6 cell) 
c_G6P_0     = 1.60E-07; % mmol/(10^6 cell) 
c_GAP_0     = 4.00E-07; % mmol/(10^6 cell) 
c_GLU_0     = 2.20E-04; % mmol/(10^6 cell) 
c_GLY_0     = 7.80E-01; % mM
c_HIS_0     = 8.50E-01; % mM
c_ILE_0     = 2.2;      % mM
c_LAC_0     = 7.50E-01; % mM
c_LEU_0     = 0;        % mM
c_LYS_0     = 2.1;      % mM
c_mAb_0     = 0;        % mg.L-1
c_MAL_0     = 1.50E-06; % mmol/(10^6 cell) 
c_MET_0     = 4.30E-01; % mM
c_NAD_0     = 6.90E-07; % mmol/(10^6 cell) 
c_NADH_0    = 8.50E-07; % mmol/(10^6 cell) 
c_NADP_0    = 3.50E-07; % mmol/(10^6 cell) 
c_NADPH_0   = 2.30E-07; % mmol/(10^6 cell) 
c_NH4_0     = 2.00E-01; % mM
c_OXA_0     = 1.00E-06; % mmol/(10^6 cell) 
c_PEP_0     = 8.80E-07; % mmol/(10^6 cell) 
c_PHE_0     = 1.1;      % mM
c_PRO_0     = 6.50E-01; % mM
c_PYR_0     = 2.00E-06; % mmol/(10^6 cell) 
c_R5P_0     = 5.40E-08; % mmol/(10^6 cell) 
c_SER_0     = 8.30E-01; % mM
c_SUC_0     = 8.00E-08; % mmol/(10^6 cell) 
c_THR_0     = 1.1;      % mM
c_TYR_0     = 1.3;      % mM
c_VAL_0     = 7.00E-01; % mM
c_Cell_0    = 0.3;      % 10^6 cells/mL


y0 = ...
    [c_ACCOA_0, c_ADP_0, c_AKG_0, c_ALA_0, c_AMP_0, c_ARG_0, c_ASN_0, c_ASP_0,...
    c_ATP_0, c_CIT_0, c_CYS_0, c_EGLC_0, c_EGLN_0, c_EGLU_0, c_F6P_0, c_G6P_0,...
    c_GAP_0, c_GLU_0, c_GLY_0, c_HIS_0, c_ILE_0, c_LAC_0, c_LYS_0, c_mAb_0,...
    c_MAL_0, c_MET_0, c_NAD_0, c_NADH_0, c_NADP_0, c_NADPH_0, c_NH4_0, c_OXA_0,...
    c_PEP_0, c_PHE_0, c_PRO_0, c_PYR_0, c_R5P_0, c_SER_0, c_SUC_0, c_THR_0,....
    c_TYR_0, c_VAL_0, c_Cell_0];




