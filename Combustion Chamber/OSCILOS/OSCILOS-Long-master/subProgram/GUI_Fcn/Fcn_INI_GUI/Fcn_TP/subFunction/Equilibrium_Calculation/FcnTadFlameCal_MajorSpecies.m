function [T2,chi,nt] = FcnTadFlameCal_MajorSpecies(indexFuel,phi,T1,eff,T2Guess,THERMO)
% 
% FcnEquilibrium_MajorSpecies -- calculate the combustion equilibrium
% state based on water-gas shift reaction.
% 
% This code accounts for 6 species in the combustion products
% consider a reaction
% CxHy + a(O2 + 3.76N2) --> b CO2 + c CO + d H20 + e H2 + f O2 + 3.76aN2
% when phi <= 1: It is considered that all the fuels are consumed and there
% are no productions of CO and H2.
% We thus have:
% a = ast/phi
% ast = x + y/4;
% b = x;
% c = 0;
% d = y/2;
% e = 0;
% f = (1-phi)*a;
% 
% when phi > 1, the oxygen is not enough to ensure all the CO and H2 to be
% consumed to CO2 and H2O. The ratio of c/e needs to be determined. We
% employe a single equlibrium reaction, CO + H2O = CO2 + H2, the so-called
% waer-gas shift reaction, to acount for the simultaneous presence of the
% incomplete products of combustion, CO and H2. 
% We thus have the relations:
% b = (2*a.*(Kp-1)+x+y./2)./2./(Kp-1)-...
%     1/(2*(Kp-1)).*((2*a.*(Kp-1)+x+y./2).^2-4*Kp.*(Kp-1).*(2*a.*x-x.^2)).^0.5;
% c = x-b;
% d = 2*a-b-x;
% e = -2*a+b+x+y./2;
% f =0;
% Kp = (b*e)/(c*d);
% where, Kp is the equilibrium constant for the water-gas shift reaction. 
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
% One can use the indexFuel to choose the fuel
% Fuel :
%    WFuel=[    16.04246 28.05316 30.06904 ...
%               44.09562 56.10632 58.12220 ...
%               58.12220 167.31102];
% Products: 
% NameProdLatex = {'CO$_2$'    'CO'    'H$_2$O'    'H$_2$'    'O$_2$'    'N$_2$'   'Air'};
% WProd=[ 44.00950 28.01010 18.01528 2.01588 31.99880 28.01340 28.96512];
% 
% --------------------------- Inputs --------------------------------------
% THERMO contains all the information to be used
% indexFuel is the index of the fuel
% phi is the equivalence ratio
% T1 is the incident mean temperature
% pi0 is the standard presure
% eff is the combustion efficiency
% T2Guess is a first guess of the flame temperature
%
% >> Tad = FcnTadFlameCal_MajorSpecies(1,1,300,1,1500,THERMO)
%
% author: Jingxuan LI (jingxuan.li@imperial.ac.uk) and Aimee S.Morgans (a.morgans@imperial.ac.uk)
% last edited: 2015-6-2
if nargin < 6
   load 'THERMO.mat' 
end
if nargin < 5
    T2Guess = 1500;
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
    disp('CH4')
end
% -------------------------------
Mix.indexFuel   = indexFuel;
Mix.phi         = phi;
Mix.T1          = T1;
Mix.p1          = THERMO.p0;
Mix.eff         = eff;
% -------------------------------
x = THERMO.Fuel.MolNumXY(indexFuel,1);               
y = THERMO.Fuel.MolNumXY(indexFuel,2);
Mix.x = x;
Mix.y = y;
% -------------------------------------------------------------------------
% calculate the flame temperature
T2ini = T2Guess; % Guess a temperature
F = @(T2)MyFcn_Tf(T1,T2,Mix,THERMO);
options = optimset('Display','off');        % the calculation results by fsolve will not be shown in the workspace
[T2,fval,exitflag] = fsolve(F,T2ini,options);
%
% -------------------------------------------------------------------------
%
[chi,nt] = FcnMole_Fraction_MajorSpecies(x,y,phi,T2,THERMO);
%
% -----------------------end of main code---------------------------------


function DeltaH = MyFcn_Tf(T1,T2,Mix,THERMO)
a = (Mix.x + Mix.y/4)./Mix.phi;
[chi,nt] = FcnMole_Fraction_MajorSpecies(Mix.x,Mix.y,Mix.phi,T2,THERMO);
% combustion products
% CO2,CO,H2O,H2,O2,N2
% -------------------- Delta hs ---------------------------
Dhs_Prod    = zeros(2,6);
for ss = 1:6
    Dhs_Prod(2,ss)  = FcnCalDelta_hs(T2,1,ss,THERMO);   
end
for ss = 5:6
    Dhs_Prod(1,ss)  = FcnCalDelta_hs(T1,1,ss,THERMO);   
end
Dhs_Fuel    = FcnCalDelta_hs(T1,2,Mix.indexFuel,THERMO); % at temperature T1
%
% --------------------- Delta hf ---------------------------
Dhf_Prod(1:6)   = THERMO.Prod.Delta_hf(1:6)';
Dhf_Fuel        = THERMO.Fuel.Delta_hf(Mix.indexFuel);
%
% --------------------- final DH ----------------------------
DeltaH  = Mix.eff.*(sum(chi(1:6).*Dhf_Prod(1:6))-1/nt.*Dhf_Fuel)...
        + (sum(chi(1:6).*Dhs_Prod(2,1:6))...
        -1/nt.*(Dhs_Fuel+a*Dhs_Prod(1,5)+3.76*a*Dhs_Prod(1,6)));
%
% ---------------------end of MyFcn_Tf ------------------------------------