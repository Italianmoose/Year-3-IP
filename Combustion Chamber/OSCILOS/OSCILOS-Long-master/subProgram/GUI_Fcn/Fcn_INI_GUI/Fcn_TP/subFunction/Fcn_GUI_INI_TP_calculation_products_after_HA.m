function [TRatio, cMean2, DeltaHr, gamma2, Cp2] =...
         Fcn_GUI_INI_TP_calculation_products_after_HA(handles,indexHA_num,TMean1,pMean1)
% This function is used to calculate the mean properties after the heat
% addition interface
%
% indexHA_num is the order of heat addition interface from the inlet to the
% outlet 
% TMean1 denotes the incident mean temperature, 
% pMean1 denotes the incident mean pressure
% TRatio is the temperature jump ratio
% cMean2 denotes the speed of sound after the interface
% DeltaHr denotes the heat addtion, or heat release rate from flame
% gamma2 denotes the specific heat ratio
% Cp2 denotes the heat capacity at constant pressure after the interface
%
% first created: 2014-12-03
% last modified: 2015-06-03
% author: Jingxuan LI (jingxuan.li@imperial.ac.uk)
global CI
ss          = indexHA_num;
HA_style    = CI.TP.HA_style(ss);                           % heat addition style of selected heating interface
switch HA_style
    case 1
        % temperature determined
        TRatio  = CI.TP.TRatio(ss);
        TMean2  = TMean1*TRatio;
        [cMean1,gamma1,Delta_hsMass1,cpMole1,cpMass,Rg1, meanMW]    = FcnCalHeatPropertiesAir(TMean1);
        [cMean2,gamma2,Delta_hsMass2,cpMole2,Cp2,Rg2, meanMW]       = FcnCalHeatPropertiesAir(TMean2);
        DeltaHr = Delta_hsMass2 - Delta_hsMass1;  
    case 2
        % fuel determined
        load THERMO.mat;
        % ------------------
        % hot gases after the heat addtion interface are considered to be
        % combustion products when FlagAir == 0, otherwise hot gases are considered
        % to be hot air
        FlagAir     = 0;
        TMean2      = FcnTadFlameCal(CI.TP.indexFuel(ss),CI.TP.Phi(ss),TMean1,CI.TP.eff(ss),pMean1,THERMO);
        x           = THERMO.Fuel.MolNumXY(CI.TP.indexFuel(ss),1);
        y           = THERMO.Fuel.MolNumXY(CI.TP.indexFuel(ss),2);
        [chi,nt]    = FcnMole_Fraction_6Species(x,y,CI.TP.Phi(ss),TMean2,pMean1,THERMO);
        [cpMass,gamma,Delta_hsMass,c,Rg,meanMW] = FcnCalHeatProperties(TMean1,TMean2,chi,FlagAir,2,THERMO);
        %
        DeltaHr     = diff(Delta_hsMass);
        cMean2      = c(2);
        Cp2         = cpMass(2);
        WProd       = meanMW(2);
        gamma2      = gamma(2);
        TRatio      = TMean2./TMean1;                           % Temperature jump ratio
    otherwise
        % Code for when there is no match.
end
assignin('base','CI',CI)
% -------------------------------------------------------------------------