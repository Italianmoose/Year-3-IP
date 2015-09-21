function [c,gamma,Delta_hsMass,cpMole,cpMass,Rg, meanMW] = FcnCalHeatPropertiesAir(T)
% FcnCalHeatPropertiesAir
% This function is used to calculate the thermodynamic properties of air
% [c,gamma,Delta_hsMass,cpMole,cpMass,Rg, meanMW] = FcnCalHeatPropertiesAir(T)
% Example :
% [c,gamma,Delta_hsMass,cpMole,cpMass,Rg, meanMW] = FcnCalHeatPropertiesAir(300)
%
% Author: Jingxuan LI jingxuan.li@imperial.ac.uk
% last revision: 2015-06-02
%
if nargin < 1
    T = 300;
    disp('T = 300 K')
end
if nargin < 2
   load 'THERMO.mat' 
end
%
R       = THERMO.R;
W       = THERMO.Prod.W(7);
Rg      = R./W.*1000;
[Delta_hsMole, Delta_hsMass] = FcnCalDelta_hs(T,1,7,THERMO);
[cpMole,cpMass] = FcnCalCp(T,1,7,THERMO);
gamma   = cpMole./(cpMole - R);
c       = (gamma.*Rg.*T).^0.5;
meanMW  = W;
%
% ------------------------------end----------------------------------------