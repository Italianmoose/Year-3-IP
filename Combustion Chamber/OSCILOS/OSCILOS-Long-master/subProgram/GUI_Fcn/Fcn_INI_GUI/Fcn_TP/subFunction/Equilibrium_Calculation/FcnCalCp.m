function [cpMole,cpMass] = FcnCalCp(T,flagGasType,index,THERMO)
% FcnCalCp  --- Calculation of heat capacity at constant pressure
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
% THERMO:           data for calculation
% 
%
% >>  [cpMole,cpMass] = FcnCalCp(2000,1,1,THERMO)
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
R = THERMO.R;
if      T <= 1000
    indexN=2*index-1;
elseif  T > 1000
    indexN=2*index;
end
switch flagGasType  
    case 1
        aT  = THERMO.Prod.a(indexN,:);
        W   = THERMO.Prod.W(index);
    case 2
        aT  = THERMO.Fuel.a(indexN,:);
        W   = THERMO.Fuel.W(index);
    otherwise
    % Code for when there is no match.
end
%         
% Cp/Ru = a1 T^-2 + a2 T^-1 + a3 + a4 T + a5 T^2 + a6 T^3 + a7 T^4
%
cpMole = R./T.^2.*polyval(aT(7:-1:1),T);
cpMass = cpMole./W.*1e3;
% 
% -------------------------- end ------------------------------------------
%