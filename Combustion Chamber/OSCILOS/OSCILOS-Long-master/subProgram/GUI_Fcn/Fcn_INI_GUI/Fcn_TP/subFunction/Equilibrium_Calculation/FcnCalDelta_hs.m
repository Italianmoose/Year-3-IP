function [Delta_hsMole, Delta_hsMass, hsMole, hsMass] = FcnCalDelta_hs(T,flagGasType,index,THERMO)

% FcnCalDelta_hs  --- Calculation of sensible enthalpy
% T:            temperature
% flagGasType:  type of gas. 
%               1: combustion products. 2: fuel
% index:        index of selected type of fuel 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 8 fuels can be precribed:
% CH4,
% C2H4,
% C2H6,
% C3H8,
% 1-butene (C4H8),
% n-butane (C4H10),
% isobutane (C4H10),
% Jet-A(C12H23, gas)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 7 products can be prescribed:
% CO2
% CO
% H2O
% H2
% O2
% N2
% Air
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THERMO: data for calculation
% 
%
% >>  [Delta_hsMole, Delta_hsMass, hsMole, hsMass] = FcnCalDelta_hs(2000,1,1,THERMO)
%
% author: Jingxuan LI (jingxuan.li@imperial.ac.uk) and Aimee S.Morgans (a.morgans@imperial.ac.uk)
%
% last revision: 2015-06-02
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 4
   load 'THERMO.mat' 
end
if nargin < 3
   index = 1;
   %
   disp('index = 1')
end
if nargin < 2
   flagGasType = 1;
   disp('CO2')
end
if nargin < 1
   T = 2000;
   disp('T = 2000 K')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R   = THERMO.R;           % R = 8.3145
T0  = THERMO.T0;          % T0 = 298.15 
if T<=1000
    indexN=2*index-1;
elseif 1000<T
    indexN=2*index;
end
switch flagGasType  
    case 1
        aT0 = THERMO.Prod.a(2*index-1,:);
        aT  = THERMO.Prod.a(indexN,:);
        bT0 = THERMO.Prod.b(2*index-1,:);
        bT  = THERMO.Prod.b(indexN,:);
        W   = THERMO.Prod.W(index);
    case 2
        aT0 = THERMO.Fuel.a(2*index-1,:);
        aT  = THERMO.Fuel.a(indexN,:);
        bT0 = THERMO.Fuel.a(2*index-1,:);
        bT  = THERMO.Fuel.a(indexN,:);
        W   = THERMO.Fuel.W(index);
    otherwise
    % Code for when there is no match.
end
%
% Delta_hs/RT = -a1 T^-2 + a2 lnT T^-1 + a3 + a4 T/2 + a5 T^2/3 + a6 T^3/4 + a7 T^4/5 + b1/T
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% hs at input temperature
temp = 0;
for ss=1:5
    temp = temp + aT(ss+2).*T.^ss./ss;
end
hs = (temp - aT(1)./T + aT(2).*log(T) + bT(1)).*R;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% hs at standard temperature
temp=0;
for ss=1:5
    temp=temp+aT0(ss+2).*T0.^ss./ss;
end
hs0             = (temp - aT0(1)./T0 + aT0(2).*log(T0) + bT0(1)).*R;
Delta_hsMole    = hs - hs0;
hsMole          = hs0;
Delta_hsMass    = Delta_hsMole/W*1e3;
hsMass          = hsMole/W*1e3;
%
% -------------------------- end ------------------------------------------