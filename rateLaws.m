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

%12
V_CITS = V_maxCITS * c_CIT/(K_mCIT+c_CIT) * c_NAD/c_NADH/(K_mNAD_NADH+c_NAD/c_NADH);
%13
V_AKGDH = V_maxAKGDH * c_AKG/(K_mAKG+c_AKG) * c_NAD/c_NADH/(K_mNAD_NADH+c_NAD/c_NADH); 
%14
V_SUC = V_maxSUC * c_SUC/(K_mSUC+c_SUC) * (c_NAD/c_NADH)/(K_mNAD_NADH+c_NAD/c_NADH) * (c_ADP/c_ATP)/(K_mADP_ATP+c_ADP_ATP);
%15
V_MLD = V_maxMLD * c_MAL/(K_mMAL+c_MAL) * (c_NAD/c_NADH)/(K_mNAD_NADH+c_NAD/c_NADH);
%16
V_PC = V_maxPC * c_PYR/(K_mPYR+c_PYR) * (c_ADP/c_ATP)/(K_mADP_ATP+c_ADP/c_ATP);
%17
V_ME = V_maxME * c_MAL/(K_mMAL+c_MAL) * c_NAD/c_NADH/(K_mADP_ATP+c_ADP/c_ATP);
%18
V_GlnT = V_fmaxGlnT * c_EGLN/(K_mEGLN+c_EGLN) * c_ATP/c_ADP/(K_mATP_ADP+c_ATP/c_ADP) - V_rmaxGlnT * c_GLU/(K_mGLU+c_GLU) * c_ADP/c_ATP/(K_mADP_ATP+c_ADP/c_ATP) * c_NH4/(K_mNH4+c_NH4);
%19
V_GLDH = V_fmaxGLDH * c_GLU/(K_mGLU+c_GLU) * c_NAD/c_NADH/(K_mNAD_NADH+c_NAD/c_NADH) - V_rmaxGLDH * c_AKG/(K_mAKG+c_AKG) * c_NADH/c_NAD/(K_mNADH_NAD+c_NADH/c_NAD) * c_NH4/(K_mNH4+c_NH4);
%20
V_AlaTA = V_fmaxAlaTA * c_GLU/(K_mGLU+c_GLU) * c_PYR/(K_mPYR+c_PYR) - V_rmaxAlaTA * c_ALA/(K_mALA+c_ALA) * c_AKG/(K_mAKG+c_AKG) * K_dEGLN/(K_dEGLN+c_EGLN);
%21
V_ASN = V_maxASN * c_ASN/(K_mASN+c_ASN);
%22
V_ASTA = V_maxASTA * c_AKG/(K_mAKG+c_AKG) * c_ASP/(K_mASP+c_ASP);
%23
V_AAtoSUC = V_maxAAtoSUC * c_LYS/(K_mLYS+c_LYS) * c_ILE/(K_mILE+c_ILE) * c_AKG/(K_mAKG+c_AKG) * c_VAL/(K_mVAL+c_VAL) * c_NAD/c_NADH/(K_mNAD_NADH+c_NAD/c_NADH) * c_NAPD/c_NAPDH/(K_mNADP_NAPDH+c_NAPD/c_NAPDH) * c_ATP/c_ADP/(K_mATP_ADP+c_ATP/c_ADP);
%24
V_HISARGTA = V_maxHISARGTA * c_HIS/(K_mHIS+c_HIS) * c_ARG/(K_mARG+c_ARG) * c_LEU/(K_mLEU+c_LEU) * c_AKG/(K_mAKG+c_AKG);
%25
V_GluT = Vf_maxGluT * c_GLU/(K_mGLU+c_GLU)-V_rmaxGluT * c_EGLU/(K_mEGLU+c_EGLU);
%26
V_SDHH = V_maxSDHH * c_SER/(K_mSER+c_SER);
%27
V_ATPase = V_maxATPase * c_ATP/(K_mATP+c_ATP);
%28
V_NADPHox = V_maxNADPHox * c_NADPH/(K_mNADPH+c_NADPH);
%29
V_resp = V_maxresp * c_ADP/c_ATP/(K_mADP_ATP+c_ADP/c_ATP) * c_NADH/(K_mNADH+c_NADH);
%30
v_leak = V_maxleak * c_NADH/(K_mNADH+c_NADH);
%31
V_AK = V_fmaxAK * c_AMP/(K_mAMP+c_AMP) * c_ATP/(K_mATP+c_ATP) -V_rmaxAK * c_ADP/(K_mADP+c_ADP);
%32
V_PPRibP = V_maxPPRibP * c_R5P/(K_mR5P+c_R5P) * c_EGLN/(K_mEGLN+c_EGLN) * c_ASP/(K_mASP+c_ASP) * c_GLY/(K_mGLY+c_GLY) * c_ATP/(K_mATP+c_ATP);
%33 There might be typos for dLAC and dNH4
V_growth = V_maxgrowth * c_G6P/(K_mG6P_growth+c_G6P) * c_R5P/(K_mR5P_growth+c_R5P) * c_LYS/(K_mLYS_growth+c_LYS) * c_CIT/(K_mCIT_growth+c_CIT) * c_ILE/(K_mILE_growth+c_ILE) * c_VAL/(K_mVAL_growth+c_VAL) * c_TYR/(K_mTYR_growth+c_TYR) * c_ASN/(K_mASN_growth+c_ASN) * c_ASP/(K_mASP_growth+c_ASP) * c_ALA/(K_mALA_growth+c_ALA) * c_ARG/(K_mARG_growth+c_ARG) * c_EGLN/(K_mGLN_growth+c_EGLN) * c_EGLU/(K_mEGLU_growth+c_EGLU) * c_GLY/(K_mGLY_growth+c_GLY) * c_HIS/K_mHIS_growth+c_HIS) * c_PHE/(K_mPHE_growth+c_PHE) * c_PRO/(K_mPRO_growth+c_PRO) * c_TYR/(K_mTYR_growth+c_TYR) * c_THR/(K_mTHR_growth+c_THR) * K_dLAC_growth/(K_dLAC_growth+c_LAC) * K_dNH4_growth/(K_dNH4_growth+c_NH4) * c_ATP/(K_mATP_growth+c_ATP); 
%34
V_mAb = V_maxmAb * c_LYS/(K_mLYS_mAb+c_LYS) * c_ILE/(K_mILE_mAb+c_ILE) * c_VAL/(K_mVAL_mAb+c_VAL)* c_TYR/(K_mTYR_mAb+c_TYR) * c_ASN/(K_mASN_mAb+c_ASN) * c_ASP/(K_mASP_mAb+c_ASP) * c_ALA/(K_mALA_mAb+c_ALA) * c_ARG/(K_mARG_mAb+c_ARG) * c_EGLN/(K_mEGLN_mAb+c_EGLN) * c_EGLU/(K_mEGLU_mAb+c_EGLU) * c_GLY/(K_mGLY_mAb+c_GLY) * c_HIS/(K_mHIS_mAb+c_HIS) * c_PHE/(K_mPHE_mAb+c_PHE) * c_PRO/(K_mPRO_mAb+c_PRO) * c_TYR/(K_mTYR_mAb+c_TYR) * c_THR/(K_mTHR_mAb+c_THR) * c_ATP/(K_mATP_mAb+c_ATP);            %EGLN typo 
    
    
    
     %rateLaws( 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
        
