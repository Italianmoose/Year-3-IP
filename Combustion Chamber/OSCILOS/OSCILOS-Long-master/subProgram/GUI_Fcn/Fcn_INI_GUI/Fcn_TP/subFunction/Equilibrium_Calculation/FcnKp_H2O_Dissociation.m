function Kp = FcnKp_H2O_Dissociation(T,THERMO)
% FcnKp_H2O_Dissociation --- calculate Kp for the reaction H2O = H2 + 0.5
% O2
% >> Kp = FcnKp_H2O_Dissociation(2000,THERMO)
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
Kp  = exp((g(2,1)*T - g(2,2))./R./T);
%
% ----------------------end------------------------------------------------
