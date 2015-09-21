function [ci,co,DeltaHr,Cp_o]=Fcn_calculation_c_q_air(Ti,To)
if nargin == 1
   To=Ti;
end
% global TF                                   % global value TF, it is recommended to clear this cell once calculation is finished
% addpath('GUI_Data')                         % add directory to search path, it includes some supportment documents
%---------------------------
%
% Fcn_initialization
load('GUI_Data.mat')
global TF                                   % global value TF, it is recommended to clear this cell once calculation is finished
%
%---------------------------
R=TF.Ru;
Delta_hs_i=Fcn_calculation_Deltahs(Ti);
Delta_hs_o=Fcn_calculation_Deltahs(To);
cp_i=Fcn_calculation_cp(Ti);
cp_o=Fcn_calculation_cp(To);
Cp_o=cp_o./TF.WProd(7).*1000;
ci=(cp_i./(cp_i-R).*R./TF.WProd(7).*1000.*Ti).^0.5;
co=(cp_o./(cp_o-R).*R./TF.WProd(7).*1000.*To).^0.5;
DeltaHr=-(Delta_hs_i-Delta_hs_o)./(TF.WProd(7)./1000);       % represent the heat value used to heat the air from Ti to To per kg
assignin('base','TF',TF); 
%
%--------------------------------------------------------------------------
%-- Delta_hs
function Delta_hs=Fcn_calculation_Deltahs(T)
global TF
R=TF.Ru;
T0=TF.T0;
Tnum=7;
if      T<=1000
    indexN=2*Tnum-1;
elseif  1000<T
    indexN=2*Tnum;
end
%
%
aT0 = TF.aCoeff_Prod(2*Tnum-1,:);
aT  = TF.aCoeff_Prod(indexN,:);
bT0 = TF.bCoeff_Prod(2*Tnum-1,:);
bT  = TF.bCoeff_Prod(indexN,:);
%
% H/RT = -a1 T^-2 + a2 lnT T^-1 + a3 + a4 T/2 + a5 T^2/3 + a6 T^3/4 + a7 T^4/5 + b1/T
temp=0;
for j=1:5
    temp=temp+aT(j+2).*T.^j./j;
end
hs=(temp-aT(1)./T+aT(2).*log(T)+bT(1)).*R;
temp=0;
for j=1:5
    temp=temp+aT0(j+2).*T0.^j./j;
end
hs0=(temp-aT0(1)./T0+aT0(2).*log(T0)+bT0(1)).*R;
Delta_hs=hs-hs0;
%
function cp=Fcn_calculation_cp(T)
global TF
R=TF.Ru;
T0=TF.T0;
Tnum=7;
if      T<=1000
    indexN=2*7-1;
elseif  1000<T
    indexN=2*7;
end
%
%
aT0 = TF.aCoeff_Prod(2*Tnum-1,:);
aT  = TF.aCoeff_Prod(indexN,:);
bT0 = TF.bCoeff_Prod(2*Tnum-1,:);
bT  = TF.bCoeff_Prod(indexN,:);
%         
% Cp/R = a1 T^-2 + a2 T^-1 + a3 + a4 T + a5 T^2 + a6 T^3 + a7 T^4
cp=R./T.^2.*polyval(aT(7:-1:1),T);
%
%-----------------
function Fcn_initialization
global TF
TF.p0 = 101325;
TF.T0 = 298.15;
TF.R  = 8.3145;
% The NASA polynomials have the form:
% 
%       Cp/R = a1 T^-2 + a2 T^-1 + a3 + a4 T + a5 T^2 + a6 T^3 + a7 T^4
%     	H/RT = -a1 T^-2 + a2 lnT T^-1 + a3 + a4 T/2 + a5 T^2/3 + a6 T^3/4 + a7 T^4/5 + b1/T
%     	S/R  = -a1 T^-2 - a2 T^-1 + a3 lnT + a4 T + a5 T^2/2 + a6 T^3/3 + a7 T^4/4 + b2
% 
% where a1, a2, a3, a4, a5, a6, and a7 are the numerical coefficients 
% supplied in NASA thermodynamic files. The first 7 numbers starting on the
% second line of each species entry (five of the second line and the first 
% two of the third line) are the seven coefficients (a1 through a7, 
% respectively) for the high-temperature range (above 1000 K, the upper 
% boundary is specified on the first line of the species entry). The
% following seven numbers are the coefficients (a1 through a7, respectively)
% for the low-temperature range (below 1000 K, the lower boundary is specified 
% on the first line of the species entry).
% 
% H in the above equation is defined as
% 
%     H(T) = Delta Hf(298) + [ H(T) - H(298) ]
% 
% so that, in general, H(T) is not equal to Delta Hf(T) and one needs to 
% have the data for the reference elements to calculate Delta Hf(T). 
% load polyfit coefficients-----------------------------------------------
filename='NASA_polyfit_species_2002.txt';      
fid=fopen(filename);
frewind(fid);                           % index the first row
cellStart = textscan(fid, '%s', 1, 'headerLines', 23);
StrStart=cellStart{1};
sprintf(StrStart{1});
%
% Products
TF.NameProdPre={'CO$_2$'    'CO'    'H$_2$O'    'H$_2$'    'O$_2$'    'N$_2$'   'Air'};
TF.WProd=[ 44.00950 28.01010 18.01528 2.01588 31.99880 28.01340 28.96512];
TF.Delta_hf_Prod=  [      -393510.000;...
                                -110535.196;...
                                -241826.000;...
                                0;...
                                0;...
                                0;...
                                -125.530];   
   
for i=1:7
    textscan(fid, '%s', 1, 'headerLines', 1);
    cellName = textscan(fid, '%s', 1, 'headerLines', 1);
    StrName=cellName{1};
    NameProd{i}=StrName{1};
    sprintf(StrName{1});
    cellGas = textscan(fid, '%f%f%f%f%f', 2, 'headerLines', 3);       % read title
    dataGas =cell2mat(cellGas);
    TF.aCoeff_Prod(2*i-1,1:5)=dataGas(1,1:5);
    TF.aCoeff_Prod(2*i-1,6:7)=dataGas(2,1:2);
    TF.bCoeff_Prod(2*i-1,1:2)=dataGas(2,4:5);
    cellGas = textscan(fid, '%f%f%f%f%f%f', 2, 'headerLines', 2);       % read title
    dataGas =cell2mat(cellGas);
    TF.aCoeff_Prod(2*i,1:5)=dataGas(1,1:5);
    TF.aCoeff_Prod(2*i,6:7)=dataGas(2,1:2);
    TF.bCoeff_Prod(2*i,1:2)=dataGas(2,4:5);
end
fclose(fid);
assignin('base','TF',TF); 

