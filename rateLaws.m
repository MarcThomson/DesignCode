function rates = rateLaws(c_ACCOA,c_ADP,c_AKG,c_ALA,c_AMP,c_ARG,c_ARK,c_ASN,...
                          c_ASP,c_ATP,c_Cell,c_CIT,c_EGLC,c_EGLN,c_EGLU,c_F6P,...
                          c_G6P,c_GAP,c_GLN,c_GLU,c_GLY,c_HIS,c_ILE,c_LAC,c_LEU,...
                          c_LYS,c_mAb,c_MAD,c_MAL,c_MET,c_NAD,c_NADH,c_NADP,...
                          c_NADPH,c_NH4,c_O2,c_OXA,c_PEP,c_PHE,c_PRO,c_PYR,...
                          c_R5P,c_SER,c_SUC,c_THR,c_TYR,c_VAL,c_X5P)
                      
parameterDefinitions;

V_KH = V_maxHK * c_EGLC * (1 + beta_AMP_ATP * c_AMP/c_ATP/(alpha_AMP_ATP*K_A_AMP_ATP)) ...
     / (K_mEGLC*(1+c_AMP/c_ATP/K_A_AMP_ATP)+c_EGLC * (1 + c_AMP/c_ATP / (alpha_AMP_ATP * K_A_AMP_ATP)))...
     * c_AMP/c_ATP / (K_mATP_ADP+c_AMP/c_ATP) * c_EGLC/(K_mEGLC + c_EGLC) * K_dG6P/(K_dG6P+c_G6P);
      %mmol/1E6 cells/hr
      
V_PGI = V_fmaxPGI * c_G6P/(K_mG6P+c_G6P)* K_dPEP/(K_dPEP+c_PEP) ...
        -V_rmaxPGI*c_F6P/(K_mG6P+c_F6P) ;     
        %mmol/1E6 cells/hr

V_PFK = V_maxPFK * c_F6P*(1 + beta_AMP_ATP * c_AMP/c_ATP / (alpha_AMP_ATP*K_A_AMP_ATP)) ...
        / (K_mF6P * (1 + c_AMP/c_ATP / (K_A_AMP_ATP)) + c_F6P * (1+c_AMP/c_ATP /(alpha_AMP_ATP*K_A_AMP_ATP))) ...
        * c_AMP/c_ATP / (K_mATP_ADP+c_AMP/c_ATP ) * K_dLAC/(K_dLAC+c_LAC);
        %mmol/1E6 cells/hr
        
V_PGK = V_maxPGK * c_GAP/ (K_mG6P+c_GAP) * c_ADP/c_ATP/(K_mADP_ATP+c_ADP/c_ATP) ... 
        * c_NAD/c_NADH / (K_mNAD_NADH+c_NAD/c_NADH); 
        %mmol/1E6 cells/hr
        
V_PK = V_maxPK * c_PEP * (1 + beta_F6P * c_F6P/ (alpha_F6P * K_A_F6P)) ...
       / (K_mPEP * (1 + c_F6P/K_A_F6P)+c_F6P*(1+c_F6P/K_A_F6P)) ...
       * c_ADP/c_ATP /(K_mADP_ATP+c_ADP/c_ATP)*K_dALA/(K_dALA*c_ALA);
        %mmol/1E6 cells/hr
 
V_LDH = V_fmaxLDH* c_PYR*(1+beta_AMP_ATP*c_AMP/c_ATP/(alpha_AMP_ATP*K_A_AMP_ATP))...
     /(K_mPYR*(1+c_AMP/c_ATP/ K_A_AMP_ATP)+ c_PYR*(1+c_AMP/c_ATP/(alpha_AMP_ATP*K_A_AMP_ATP)))...
     *c_NADH/c_NAD/(K_mNADH_NAD+c_NADH/c_NAD)...
     -V_rmaxLDH * c_LAC/(K_mLAC+c_LAC)*c_NAD/c_NADH/(K_mNAD_NADH+c_NAD/c_NADH)...
     *K_dPYR/(K_dPYR+c_PYR);
     %mmol/1E6 cells/hr

V_G6PDH = V_maxG6PDH * c_G6P/(K_mG6P+c_G6P) * c_NADP/c_NADPH/ (K_mNADP_NADPH+c_NADP/c_NADPH);
     %mmol/1E6 cells/hr
     
V_EP = V_maxEP * c_R5P/(K_mR5P+c_R5P);
    %mmol/1E6 cells/hr
 
V_TK = V_maxTK * c_R5P/(K_mR5P+c_R5P) * c_X5P/(K_mX5P+c_X5P);
    %mmol/1E6 cells/hr
    
V_PDH = V_maxPDH * c_PYR/(K_mPYR+c_PYR) * c_NAD/c_NADH/(K_mNAD_NADH+c_NAD/NADH); % EQ 10
   
V_CS=V_maxCS * c_OXA/(K_mOXA+c_OXA) * c_ACCOA/(K_mACCOA+c_ACCOA);

%V_ not sure about equation 12

V_AKGDH=V_max_AKGDH * c_AKG/(K_mAKG+c_AKG) * c_NAD/c_NADH/(K_mNAD_NADH+c_NAD/c_NADH)

    
    
    
    
    
    
    
    
     %rateLaws( 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
        
