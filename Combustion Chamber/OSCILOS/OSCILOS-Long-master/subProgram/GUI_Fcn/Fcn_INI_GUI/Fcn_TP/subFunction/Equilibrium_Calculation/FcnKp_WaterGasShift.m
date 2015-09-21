function Kp = FcnKp_WaterGasShift(T,THERMO)
%
%  FcnKp_WaterGasShift -- calculate Kp of water-gas shift
%  reaction for a given temperature T
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
Kp  = exp((g(3,1)*T - g(3,2))./R./T);
%
% ----------------------end------------------------------------------------