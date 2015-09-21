function Fcn_TD_INI_constant_variables
% This function is used to pre-calculate some constant values
%
% last edited: 2014-11-13 09:22
%
global CI
for ss = 1:CI.TP.numSection
    CI.TP.RhoCU(ss)     = CI.TP.rho_mean(1,ss).*CI.TP.c_mean(1,ss).*CI.TP.u_mean(1,ss);  % rho*c*u, in case to calculate uRatio = (AP - AM)/(rho*c*u)
end
% -----------------------------end-----------------------------------------