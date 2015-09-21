function Kp = FcnKp_CO2_Dissociation(T,THERMO)
% FcnKp_CO2Dissociation --- calculate Kp for the reaction CO2 = CO + 0.5
% O2
% >> Kp = FcnKp_CO2_Dissociation(2000,THERMO)
if nargin < 2
   load 'THERMO.mat' 
end
if nargin < 1
   T = 2000;
   %
   disp('T = 2000 K')
end
R   = THERMO.R;
g   = THERMO.gKp;
Kp  = exp((g(1,1)*T - g(1,2))./R./T);
%
% ----------------------end------------------------------------------------
