function [cpMass,gamma,Delta_hsMass,c,Rg,meanMW] = FcnCalHeatProperties(T1,T2,chi,FlagAir,numOut,THERMO)
% FcnCalHeatProperties
% [cpMass,gamma,Delta_hsMass,c,Rg,meanMW] = FcnCalHeatProperties(T1,T2,chi,FlagAir)
% This function is used to calculate the heat properties before and after
% the combustion.
% The fresh mixture is considered as air
% The combustion product is considered as burned gases or air
% Example: [cpMass,gamma,Delta_hsMass,c,Rg,meanMW] = FcnCalHeatProperties(300,2000)
% Author : Jingxuan Li jingxuan.li@imperial.ac.uk
% last revision: 2015-06-02
%
if nargin < 1
    T1 = 300;
    disp('T1 = 300 K')
end
if nargin < 2
    T2 = T1;
    disp('T2 = T1')
end
if nargin < 3
    chi = zeros(1,6);
    chi(5) = 1/4.76;
    chi(6) = 3.76/4.76;
    disp('Air composition')
    FlagAir = 1;  % combustion product is considered as air
end
if nargin < 4
    FlagAir = 1;
    disp('Air assumption')
end
if nargin < 5
    numOut = 2;
    disp('number of output set')
end
if nargin < 6
   load 'THERMO.mat' 
end
if numOut == 1
    [c(1),gamma(1),Delta_hsMass(1),cpMole(1),cpMass(1),Rg(1), meanMW(1)] = FcnCalHeatPropertiesAir(T1);
else
    switch FlagAir
        case 0
            %%%%%%%
            % fresh mixture
            [c(1),gamma(1),Delta_hsMass(1),cpMole(1),cpMass(1),Rg(1), meanMW(1)] = FcnCalHeatPropertiesAir(T1);
            % combustion products
            for ss = 1:6
                Delta_hsProd(ss) = FcnCalDelta_hs(T2,1,ss,THERMO);
                cpMoleProd(ss)  = FcnCalCp(T2,1,ss,THERMO);
            end
            R                   = THERMO.R;
            WProd               = THERMO.Prod.W(1:6);
            meanMW(2)           = sum(WProd(1:6)'.*chi(1:6));
            Rg(2)               = R*1e3/meanMW(2);
            cpMole(2)           = sum(cpMoleProd(1:6).*chi(1:6));
            cvMole(2)           = cpMole(2) - R;
            cpMass(2)           = cpMole(2)./meanMW(2).*1e3;
            gamma(2)            = cpMole(2)./cvMole(2);
            c(2)                = (gamma(2)*R./meanMW(2).*1e3.*T2).^0.5;
            Delta_hsMole        = sum(Delta_hsProd(1:6).*chi(1:6));
            Delta_hsMass(2)    = Delta_hsMole./meanMW(2).*1e3;
        case 1
            [c(1),gamma(1),Delta_hsMass(1),cpMole(1),cpMass(1),Rg(1), meanMW(1)] = FcnCalHeatPropertiesAir(T1);
            [c(2),gamma(2),Delta_hsMass(2),cpMole(2),cpMass(2),Rg(2), meanMW(2)] = FcnCalHeatPropertiesAir(T2);
    end
end
%
%
% ------------------------------- end -------------------------------------