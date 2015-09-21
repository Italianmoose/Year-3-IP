function [chi,nt] = FcnMole_Fraction_MajorSpecies(x,y,phi,T,THERMO)
%
% FcnMole_Fraction_MajorSpecies -- calculate mole fractions of the 
% combustion products of water-gas shift reaction for a given temperature T
% {'CO2','CO','H2O','H2','O2','N2'};
% >> FcnMole_Fraction_MajorSpecies(1,4,1.5,2000,TF)
if nargin < 5
   load 'THERMO.mat'
end
if nargin < 4
    T = 2000;
    disp('T = 2000 K')
end
if nargin < 3
    phi = 1.5;
    disp('phi = 1.5')
end
if nargin < 2
   y = 4;
   %
   disp('y = 4')
end
if nargin < 1
   x = 1;
   %
   disp('x = 1')
   disp('CH4 is used as the fuel')
end
%------------------------------------------------
ast = x+y/4;   % stochiometric mole number         
a   = ast/phi;
if phi<=1.0
% lean flame
    nt = x+y./2+(1-phi+3.76).*a;
    %-
    chi(1) = x./nt;
    chi(2) = 0;
    chi(3) = (y./2)./nt;
    chi(4) = 0;
    chi(5) = (1-phi).*a./nt;
    chi(6) = 3.76*a./nt;
else
% rich flame
    Kp = FcnKp_WaterGasShift(T,THERMO);
    nt = x+y./2+3.76*a;
    b = (2*a.*(Kp-1)+x+y./2)./2./(Kp-1)-...
        1/(2*(Kp-1)).*((2*a.*(Kp-1)+x+y./2).^2-4*Kp.*(Kp-1).*(2*a.*x-x.^2)).^0.5;
    chi(1) = b./nt;
    chi(2) = (x-b)./nt;
    chi(3) = (2*a-b-x)./nt;
    chi(4) = (-2*a+b+x+y./2)./nt;
    chi(5) = 0;
    chi(6) = 3.76*a./nt;
end
% ----------------------------------end------------------------------------