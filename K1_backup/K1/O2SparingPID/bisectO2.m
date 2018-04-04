function QO2 = bisectO2(QCO2,Qtotal,C_O2,C_CO2,Dt,Nrad,t,shift,drop)

tol = 0.001;

%left bound
yCO2 = QCO2/Qtotal;
minQO2 = 0.2095*(1-yCO2)*Qtotal;
yO2 = minQO2/Qtotal;


CO2Frac = (yCO2-400E-6)/(1-400E-6);
airFrac = (1-yO2-CO2Frac)/(1-0.2095);

[FO2, CO2] = ...
    spargerV2(C_O2, C_CO2, airFrac, CO2Frac, Qtotal, Dt ,Nrad, t, shift);
minLoss = FO2+drop;    

%right bound
maxQO2 = (1-yCO2)*Qtotal;
yO2 = minQO2/Qtotal;
CO2Frac = (yCO2-400E-6)/(1-400E-6);
airFrac = (1-yO2-CO2Frac)/(1-0.2095);
[FO2, CO2] = ...
    spargerV2(C_O2, C_CO2, airFrac, CO2Frac, Qtotal, Dt ,Nrad, t, shift);
maxLoss = FO2+drop;   


if maxLoss*minLoss>0
    if minLoss > 0
        QO2 = minQO2;
    elseif minLoss<0
        QO2 = maxLoss;
    end
else
    a = minQO2;
    b = maxQO2;
    m = (maxQO2+minQO2)/2;
    while abs(FO2 + drop)>tol
        m = (a+b)/2;
        %A
        yO2 = a/Qtotal;
        CO2Frac = (yCO2-400E-6)/(1-400E-6);
        airFrac = (1-yO2-CO2Frac)/(1-0.2095);
        [FO2, CO2] = ...
        spargerV2(C_O2, C_CO2, airFrac, CO2Frac, Qtotal, Dt ,Nrad, t, shift);
        fa = FO2+drop;  
        
        %b
        yO2 = b/Qtotal;
        CO2Frac = (yCO2-400E-6)/(1-400E-6);
        airFrac = (1-yO2-CO2Frac)/(1-0.2095);
        [FO2, CO2] = ...
        spargerV2(C_O2, C_CO2, airFrac, CO2Frac, Qtotal, Dt ,Nrad, t, shift);
        fb =  FO2+drop; 
        %m
        yO2 = m/Qtotal;
        CO2Frac = (yCO2-400E-6)/(1-400E-6);
        airFrac = (1-yO2-CO2Frac)/(1-0.2095);
        [FO2, CO2] = ...
        spargerV2(C_O2, C_CO2, airFrac, CO2Frac, Qtotal, Dt ,Nrad, t, shift);
        fm =  FO2+drop; 
        
        if fm*fa>0
            a = m;
        else
            b = m;
        end 
        
    end
    QO2 = m;
end
    