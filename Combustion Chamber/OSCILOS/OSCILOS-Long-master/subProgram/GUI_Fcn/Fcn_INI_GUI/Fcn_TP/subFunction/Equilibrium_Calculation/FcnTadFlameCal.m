function Tf = FcnTadFlameCal(indexFuel,phi,T1,eff,p1,THERMO)
% 
% FcnTadFlameCal
%
% This function is used to calculate the flame temperature of
% hydrocarbon-air flame
%
% 
% This code accounts for 6 species in the combustion products
% consider a reaction
% CxHy + a(O2 + 3.76N2) --> b CO2 + c CO + d H20 + e H2 + f O2 + 3.76a N2
%
% -------------------- Fuels and products ---------------------------------
% 8 fuels can be precribed:
% CH4,
% C2H4,
% C2H6,
% C3H8,
% 1-butene (C4H8),
% n-butane (C4H10),
% isobutane (C4H10),
% Jet-A(C12H23, gas)
%
% 
% --------------------------- Inputs --------------------------------------
% THERMO contains all the information to be used
% indexFuel is the index of the fuel
% phi is the equivalence ratio
% T1 is the incident mean temperature
% pi0 is the standard presure
% eff is the combustion efficiency
%
% >> Tf = FcnTadFlameCal(1,1,300,1,1e5,THERMO)
%
% author: Jingxuan LI (jingxuan.li@imperial.ac.uk) and Aimee S.Morgans (a.morgans@imperial.ac.uk)
% last edited: 2015-6-2

if nargin < 6
   load 'THERMO.mat' 
end
if nargin < 5
    p1 = 101325;
    disp('ambient pressure = 101325 Pa')
end
if nargin < 4
    eff = 1;
    disp('combustion efficiency eta = 1')
end
if nargin < 3
    T1 = 298.15;
    disp('initial temperature T1 = 298.15 K')
end
if nargin < 2
    phi = 1;
    disp('equivalence ratio phi = 1')
end
if nargin < 1
    indexFuel = 1;
    disp('CH4 is used as the fuel')
end
%
Tf = Fcn_precise_calculation(indexFuel,phi,T1,eff,p1,THERMO);
%

function Tf = Fcn_precise_calculation(indexFuel,phi,T1,eff,p1,THERMO)

% Rough prediction of flame temperature
T2Guess = 2e3;
[Tf0,chi0,nt0] = FcnTadFlameCal_MajorSpecies(indexFuel,phi,T1,eff,T2Guess,THERMO);
% -------------------------------
Mix.indexFuel   = indexFuel;
Mix.phi         = phi;
Mix.T1          = T1;
Mix.p1          = p1;
Mix.eff         = eff;
% -------------------------------
x = THERMO.Fuel.MolNumXY(indexFuel,1);               
y = THERMO.Fuel.MolNumXY(indexFuel,2);
Mix.x = x;
Mix.y = y;
Mix.chi0 = chi0;
Mix.nt0 = nt0;
% -------------------------------------------------------------------------
dTini = 1; % Guess a temperature
F = @(dT)MyFun(dT,Tf0,Mix,THERMO);

options = optimset('Display','off');        % the calculation results by fsolve will not be shown in the workspace
[DeltaTf,   fval,  exitflag] = fsolve(F,dTini,options);
Tf = Tf0 + DeltaTf;
%


function F = MyFun(dT,Tf0,Mix,THERMO)

nu0     = Mix.chi0.*Mix.nt0;
nut0    = Mix.nt0;
alpha1  = 2.78e5;
alpha2  = 2.51e5;
CpAppr  = [60, 36, 51, 34, 37.7, 35.8];
alpha3  = sum(nu0.*CpAppr);
R       = THERMO.R;
Kp10    = exp((THERMO.gKp(1,1)*Tf0 - THERMO.gKp(1,2))./R./Tf0);
Kp30    = exp((THERMO.gKp(3,1)*Tf0 - THERMO.gKp(3,2))./R./Tf0);
k1      = (2*nut0.*nu0(1).^2.*THERMO.p0/Mix.p1).^(1/3);
d1      = k1*Kp10.^(2/3);
d2      = nu0(3)/nu0(1)*Kp30;
% d2      = Mix.y/2/Mix.x*Kp30;
% d2      = Mix.y/2/Mix.x*0.2;
beta    = d1*exp(2/3*THERMO.gKp(1,2)*dT./(R*Tf0.^2));
if Mix.phi <=1
    xi1 = beta.*(3*beta + 2*nu0(5))./(3*beta + 4*nu0(5));
    xi2 = d2*xi1;
else 
    xi1 = 2*beta.^2./(2*beta + nu0(2)); 
    xi2 = d2*xi1;
end
F       = alpha1*xi1 + alpha2*xi2 + alpha3*dT;
%
% --------------------------- end -----------------------------------------




