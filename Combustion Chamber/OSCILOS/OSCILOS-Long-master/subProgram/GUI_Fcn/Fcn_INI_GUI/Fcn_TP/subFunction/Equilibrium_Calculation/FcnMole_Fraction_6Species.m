function [chi,nt] = FcnMole_Fraction_6Species(x,y,phi,T,p1,THERMO)
%
% FcnMole_Fraction_6Species -- calculate mole fractions of the 
% combustion products.
% A temperature is required to calculate Kp.
% The ambient pressure is also required. 
%  accounted gases:
% {'CO2','CO','H2O','O2','N2'};
% >> [chi,nt] = FcnMole_Fraction_6Species(x,y,phi,T,p1,THERMO)
if nargin < 6
   load 'THERMO.mat' 
end
if nargin < 5
    p1 = 101325; 
    disp('p = 101325 Pa')
end
if nargin < 4
    T = 2000; 
    disp('T = 2000 K')
end
if nargin < 3
    phi = 1.0;
    disp('phi = 1.0')
end
if nargin < 2
    y = 4;
    disp('y = 4')
end
if nargin < 1
    x = 1;
    disp('x = 1')
    disp('CH4 is used as the fuel')
end
%------------------------------------------------
ast     = x+y/4;   % stochiometric mole number         
a       = ast/phi;
pRatio  = p1/THERMO.p0;
if phi <= 1
% lean flame
% mole fractions from the major species assumption
    b = x;              % CO2
    c = 0;              % CO
    d = y/2;            % H20
    e = 0;              % H2
    f = (1 - phi)*a;    % O2
    g = 3.76*a;         % N2
    nt= x + y/2 + (4.76-phi).*a;
% add correction
    Kp1 = FcnKp_CO2_Dissociation(T,THERMO);
    h = (Kp1*b)^2.*nt/pRatio;
    Delta1 = (2*h)^2/4 - 8/27*f.^3.*(2*h);
    eta1 = -2*f/3 + (-8/27*f.^3+h + Delta1.^0.5).^(1/3) + (-8/27*f.^3+h - Delta1.^0.5).^(1/3);
    Kp3 = FcnKp_WaterGasShift(T,THERMO);
    eta2 = eta1*d*Kp3/b;
    b = b - eta1;
    c = c + eta1;
    d = d - eta2;
    e = e + eta2;
    f = f + 0.5*eta1+ 0.5*eta2;
    nt = b + c + d + e + f + g;
    % ------------------------------
    chi(1) = b/nt;
    chi(2) = c/nt;
    chi(3) = d/nt;
    chi(4) = e/nt;
    chi(5) = f/nt;
    chi(6) = g/nt;
else
% rich flame
% mole fractions from the major species assumption
    Kp3 = FcnKp_WaterGasShift(T,THERMO);
    nt  = x+y./2+3.76*a;
    b   = (2*a.*(Kp3-1)+x+y./2)./2./(Kp3-1)-...
        1/(2*(Kp3-1)).*((2*a.*(Kp3-1)+x+y./2).^2-4*Kp3.*(Kp3-1).*(2*a.*x-x.^2)).^0.5;
    c   = x - b;
    d   = (2*a-b-x);
    e   = (-2*a+b+x+y./2); 
    f   = 0;
    g   = 3.76*a;
% add correction 
    Kp1 = FcnKp_CO2_Dissociation(T,THERMO);
    h = (Kp1*b)^2.*nt/pRatio;
    % -----------------------------
    Delta1 = (2*h)^2/4 + 1/27*c.^3.*(2*h);
    eta1 = -2*c/3 + (1/27*c.^3+h + Delta1.^0.5).^(1/3) + (1/27*c.^3+h - Delta1.^0.5).^(1/3);
    eta2 = (Kp3*d)/(b)*eta1;
    % ------------------------------
    b = b - eta1;
    c = c + eta1;
    d = d - eta2;
    e = e + eta2;
    f = f + 0.5*eta1 + 0.5*eta2;
    nt = b + c + d + e + f + g;
%     % -------------------------------
    chi(1) = b/nt;
    chi(2) = c/nt;
    chi(3) = d/nt;
    chi(4) = e/nt;
    chi(5) = f/nt;
    chi(6) = g/nt;
end
% ----------------------------------end------------------------------------