function [Tf,chi,DeltaHr,c_mean,Cp_mean,WProd]=Fcn_calculation_combustion_products(indexFuel,Phi,Ti,pi0,eff,dil)
global TF                                   % global value TF, it is recommended to clear this cell once calculation is finished
% addpath('GUI_Data')                         % add directory to search path, it includes some supportment documents
%---------------------------
%s
% Fcn_initialization
load('GUI_Data.mat')
%
%---------------------------
TF.indexMethod=2;                           % the second method accounts for the dissociation
TF.Phi=Phi./(1+dil);                        % equivalence ratio, corrected by the dilution rate
TF.indexFuel=indexFuel;                     % index of fuel
TF.Ti=Ti;                                   % incident flow temperature
TF.pi0=pi0;                                 % incident pressure
TF.eff=eff;                                 % combustion efficiency
TF.dil=dil;                                 % dilution rate by the downstream air flow or pilot air flow
TF.x=TF.CHxy(TF.indexFuel,1);               
TF.y=TF.CHxy(TF.indexFuel,2);
TF.a=(TF.x+TF.y./4)./TF.Phi;                % The air used for the combustion
TF.a_dil=TF.dil.*TF.a;                      % The air used to dilute the combustion products
Ru=TF.Ru;
assignin('base','TF',TF); 
%
%-- calculate the flame temperature without dilution
Tini=2000;
options= optimset('Display','off');
TF.Tf_wodil = fsolve(@MyFcn_Tf,Tini, options);               % flame temperature without dilution
assignin('base','TF',TF); 

Tf=TF.Tf_wodil;
chi=TF.chi;
%--
Delta_hr=Fcn_Hr_calculation;
n_air=TF.a*4.76;
W_air=TF.WProd(7);
% DeltaHr(1)=Delta_hr./((n_air.*W_air+TF.WFuel(indexFuel))./1000);
% DeltaHr(2)=Delta_hr./(n_air.*W_air./1000);
% DeltaHr(3)=Delta_hr./(TF.WFuel(indexFuel)./1000);
DeltaHr=Delta_hr./(n_air.*W_air./1000);
%DeltaHr represents the heat release per kg of (mixture, air, fuel)

% TF.WProd=[ 44.00950 28.01010 18.01528 2.01588 31.99880 28.01340 28.96512];
% TF.WFuel=[    16.04246 28.05316 30.06904 ...
%                     44.09562 56.10632 58.12220 ...
%                     58.12220 167.31102];
%
cp(1)=Fcn_calculation_cp(Tf,1,1);
cp(2)=Fcn_calculation_cp(Tf,1,2);
cp(3)=Fcn_calculation_cp(Tf,1,3);
cp(4)=Fcn_calculation_cp(Tf,1,4);
cp(5)=Fcn_calculation_cp(Tf,1,5);
cp(6)=Fcn_calculation_cp(Tf,1,6);
%
WProd=sum(chi(1:6).*TF.WProd(1:6));
cp_mean=sum(chi(1:6).*cp(1:6));
c_mean=(cp_mean./(cp_mean-Ru).*Ru.*Tf./WProd.*1000).^0.5; % speed of sound in the combustion products
Cp_mean=cp_mean./WProd.*1000;
%
assignin('base','TF',TF); 
%
%
function DeltaH=MyFcn_Tf(Tf)
global TF
indexMethod=TF.indexMethod;
if indexMethod==1
    Fcn_calculation_species_mole_fraction_method1(Tf);
elseif indexMethod==2
    Fcn_calculation_species_mole_fraction_method2(Tf);
end
assignin('base','TF',TF); 
chi=TF.chi;
nt=TF.nt;
a=TF.a;
Delta_hs_Prod(1)=Fcn_calculation_Deltahs(Tf,1,1);   % CO2
Delta_hs_Prod(2)=Fcn_calculation_Deltahs(Tf,1,2);   % CO
Delta_hs_Prod(3)=Fcn_calculation_Deltahs(Tf,1,3);   % H2O
Delta_hs_Prod(4)=Fcn_calculation_Deltahs(Tf,1,4);   % H2
Delta_hs_Prod(5)=Fcn_calculation_Deltahs(Tf,1,5);   % O2
Delta_hs_Prod(6)=Fcn_calculation_Deltahs(Tf,1,6);   % N2
%
Delta_hs_ProdTi(5)=Fcn_calculation_Deltahs(TF.Ti,1,5);   % O2
Delta_hs_ProdTi(6)=Fcn_calculation_Deltahs(TF.Ti,1,6);   % N2
%
Delta_hf_Prod=TF.Delta_hf_Prod;
%
Delta_hs_Fuel   =Fcn_calculation_Deltahs(Tf,2,TF.indexFuel);
Delta_hs_FuelTi =Fcn_calculation_Deltahs(TF.Ti,2,TF.indexFuel);
%
Delta_hf_Fuel=TF.Delta_hf_Fuel(TF.indexFuel);
%
DeltaH  = TF.eff.*(sum(chi(1:6).*Delta_hf_Prod(1:6)')-1/nt.*Delta_hf_Fuel)...
        + (sum(chi(1:6).*Delta_hs_Prod(1:6))...
        -1/nt.*(Delta_hs_FuelTi+a*Delta_hs_ProdTi(5)+3.76*a*Delta_hs_ProdTi(6)));
    
function DeltaHr=Fcn_Hr_calculation
global TF
assignin('base','TF',TF); 
chi=TF.chi;
nt=TF.nt;
a=TF.a;
%
Delta_hf_Prod=TF.Delta_hf_Prod;
%
Delta_hf_Fuel=TF.Delta_hf_Fuel(TF.indexFuel);
%
DeltaHr = -nt.*TF.eff.*(sum(chi(1:6).*Delta_hf_Prod(1:6)')-1/nt.*Delta_hf_Fuel);
        
%   
function Fcn_calculation_species_mole_fraction_method1(Tf)
global TF
a=TF.a;
x=TF.x;
y=TF.y;
Phi=TF.Phi;
pi0=TF.pi0;
%------------------------------------------------
if TF.Phi<=1.0
% lean flame
    TF.nt=x+y./2+(1-Phi+3.76).*a;
    nt=TF.nt;
    %-
    TF.chi(1)=x./nt;
    TF.chi(2)=0;
    TF.chi(3)=(y./2)./nt;
    TF.chi(4)=0;
    TF.chi(5)=(1-Phi).*a./nt;
    TF.chi(6)=3.76*a./nt;
else
% rich flame
    Kp=Fcn_calculation_Kp_water_gas_shift_reaction(Tf);
    TF.nt=x+y./2+3.76*a;
    nt=TF.nt;
    b=(2*a.*(Kp-1)+x+y./2)./2./(Kp-1)-...
        1/(2*(Kp-1)).*((2*a.*(Kp-1)+x+y./2).^2-4*Kp.*(Kp-1).*(2*a.*x-x.^2)).^0.5;
    TF.chi(1)=b./nt;
    TF.chi(2)=(x-b)./nt;
    TF.chi(3)=(2*a-b-x)./nt;
    TF.chi(4)=(-2*a+b+x+y./2)./nt;
    TF.chi(5)=0;
    TF.chi(6)=3.76*a./nt;
end
assignin('base','TF',TF); 
%--
function Fcn_calculation_species_mole_fraction_method2(Tf)
global TF
a=TF.a;
x=TF.x;
y=TF.y;
Phi=TF.Phi;
pi0=TF.pi0;
[CI.varep1 CI.varep2]=Fcn_Kp_Lean_Flame(Tf);
varep1=CI.varep1;
varep2=CI.varep2;
%------------------------------------------------
% for all flame
    TF.nt=x+y./2+(1-Phi+3.76).*a+0.5*x*varep1+0.25*y*varep2;
    nt=TF.nt;
    %-
    TF.chi(1)=x.*(1-varep1)./nt;
    TF.chi(2)=x.*varep1./nt;
    TF.chi(3)=(y./2).*(1-varep2)./nt;
    TF.chi(4)=(y./2).*varep2./nt;
    TF.chi(5)=((1-Phi).*a+0.5*x*varep1+0.25*y*varep2)./nt;
    TF.chi(6)=3.76*a./nt;
% else
assignin('base','TF',TF); 
%
function [vp1_solu, vp2_solu]=Fcn_Kp_Lean_Flame(Tf)
global TF
a=TF.a;
x=TF.x;
y=TF.y;
Phi=TF.Phi;
pi0=TF.pi0;
p0=TF.p0;
p_ratio=pi0./p0;
%
Kp1=Fcn_calculation_Kp_CO2_dissociation(Tf);
Kp2=Fcn_calculation_Kp_H2O_dissociation(Tf);
syms vp1 vp2 
[vp1, vp2] =...
 solve(vp1/(1-vp1)*(((1-Phi)*a+0.5*x*vp1+0.25*y*vp2)...
                    /(x+0.5*y+(4.76-Phi)*a+0.5*x*vp1+0.25*y*vp2))^0.5.*p_ratio.^0.5==Kp1,...
       vp2/(1-vp2)*(((1-Phi)*a+0.5*x*vp1+0.25*y*vp2)...
                    /(x+0.5*y+(4.76-Phi)*a+0.5*x*vp1+0.25*y*vp2))^0.5.*p_ratio.^0.5==Kp2, vp1, vp2);                
vp1_soluP=double(vp1);
vp2_soluP=double(vp2);
vp1_solu=vp1_soluP(1);
vp2_solu=vp2_soluP(1); 
%
%------------------------------------------------
function [vp1, vp2, b]=Fcn_Kp_rich_Flame(Tf)
global TF
a=TF.a;
x=TF.x;
y=TF.y;
Phi=TF.Phi;
%
Kp1=Fcn_calculation_Kp_CO2_dissociation(Tf);
Kp2=Fcn_calculation_Kp_H2O_dissociation(Tf);
Kp=Fcn_calculation_Kp_water_gas_shift_reaction(Tf);
syms vp1 vp2 b
[vp1, vp2, b] =...
 solve((x-b*(1-vp1))/(b*(1-vp1))*...
        ((0.5*b*vp1+0.5*(2*a-b-x)*vp2)/(x+0.5*y+3.76*a+0.5*b*vp1+0.5*(2*a-b-x)*vp2))^0.5==Kp1,...
       ((0.5*y-(2*a-b-x)*(1-vp2))/((2*a-b-x)*(1-vp2)))*...
        ((0.5*b*vp1+0.5*(2*a-b-x)*vp2)/(x+0.5*y+3.76*a+0.5*b*vp1+0.5*(2*a-b-x)*vp2))^0.5==Kp2,...
        ((0.5*y-(2*a-b-x)*(1-vp2))/((2*a-b-x)*(1-vp2)))/((x-b*(1-vp1))/(b*(1-vp1)))==Kp,...
        vp1, vp2, b);
%--------------------------------------------------------------------------
%-- Delta_hs
function Delta_hs=Fcn_calculation_Deltahs(T,TSpec,Tnum)
global TF
Ru=TF.Ru;
T0=TF.T0;
if      T<=1000
    indexN=2*Tnum-1;
elseif  1000<T
    indexN=2*Tnum;
end
aa=indexN;
switch TSpec  
    case 1
        aT0 = TF.aCoeff_Prod(2*Tnum-1,:);
        aT  = TF.aCoeff_Prod(indexN,:);
        bT0 = TF.bCoeff_Prod(2*Tnum-1,:);
        bT  = TF.bCoeff_Prod(indexN,:);
    case 2
        aT0 = TF.aCoeff_Fuel(2*Tnum-1,:);
        aT  = TF.aCoeff_Fuel(indexN,:);
        bT0 = TF.bCoeff_Fuel(2*Tnum-1,:);
        bT  = TF.bCoeff_Fuel(indexN,:);
    otherwise
    % Code for when there is no match.
end
%
% H/RuT = -a1 T^-2 + a2 lnT T^-1 + a3 + a4 T/2 + a5 T^2/3 + a6 T^3/4 + a7 T^4/5 + b1/T
temp=0;
for j=1:5
    temp=temp+aT(j+2).*T.^j./j;
end
hs=(temp-aT(1)./T+aT(2).*log(T)+bT(1)).*Ru;
temp=0;
for j=1:5
    temp=temp+aT0(j+2).*T0.^j./j;
end
hs0=(temp-aT0(1)./T0+aT0(2).*log(T0)+bT0(1)).*Ru;
Delta_hs=hs-hs0;
%
function cp=Fcn_calculation_cp(T,TSpec,Tnum)
global TF
Ru=TF.Ru;
T0=TF.T0;
if      T<=1000
    indexN=2*Tnum-1;
elseif  1000<T
    indexN=2*Tnum;
end
aa=indexN;
switch TSpec  
    case 1
        aT0 = TF.aCoeff_Prod(2*Tnum-1,:);
        aT  = TF.aCoeff_Prod(indexN,:);
        bT0 = TF.bCoeff_Prod(2*Tnum-1,:);
        bT  = TF.bCoeff_Prod(indexN,:);
    case 2
        aT0 = TF.aCoeff_Fuel(2*Tnum-1,:);
        aT  = TF.aCoeff_Fuel(indexN,:);
        bT0 = TF.bCoeff_Fuel(2*Tnum-1,:);
        bT  = TF.bCoeff_Fuel(indexN,:);
    otherwise
    % Code for when there is no match.
end
%         
% Cp/Ru = a1 T^-2 + a2 T^-1 + a3 + a4 T + a5 T^2 + a6 T^3 + a7 T^4
cp=Ru./T.^2.*polyval(aT(7:-1:1),T);
%
%-----------------
function Fcn_initialization
global TF
TF.p0 = 101325;
TF.T0 = 298.15;
TF.Ru  = 8.3145;
% The NASA polynomials have the form:
% 
%       Cp/Ru = a1 T^-2 + a2 T^-1 + a3 + a4 T + a5 T^2 + a6 T^3 + a7 T^4
%     	H/RuT = -a1 T^-2 + a2 lnT T^-1 + a3 + a4 T/2 + a5 T^2/3 + a6 T^3/4 + a7 T^4/5 + b1/T
%     	S/Ru  = -a1 T^-2 - a2 T^-1 + a3 lnT + a4 T + a5 T^2/2 + a6 T^3/3 + a7 T^4/4 + b2
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
%
% Fuels
TF.NameFuelPre={      'CH$_4$','C$_2$H$_4$','C$_2$H$_6$','C$_3$H$_8$',...
                            'C$_4$H$_8$,1-butene','C$_4$H$_{10}$,n-butane',...
                            'C$_4$H$_{10}$,isobutane','C$_{12}$H$_{23}$,Jet-A(gas)'};
TF.WFuel=[    16.04246 28.05316 30.06904 ...
                    44.09562 56.10632 58.12220 ...
                    58.12220 167.31102];
TF.Delta_hf_Fuel=  [      -74800.000;...
                                52500.000;...
                                -83851.544;...
                                -104680.000;...
                                -540.000;...
                                -125790.000;...
                                -92200.000;...
                                -249657.000];
TF.CHxy=         [      1,4;...
                        2,4;...
                        2,6;...
                        3,8;...
                        4,8;...
                        4,10;...
                        4,10;...
                        12,23];
assignin('base','TF',TF); 
cellStart = textscan(fid, '%s', 1, 'headerLines', 2);
StrStart=cellStart{1};
sprintf(StrStart{1});
for i=1:8
    textscan(fid, '%s', 1, 'headerLines', 1);
    cellName = textscan(fid, '%s', 1, 'headerLines', 1);
    StrName=cellName{1};
    NameFuel{i}=StrName{1};
    sprintf(StrName{1});
    cellGas = textscan(fid, '%f%f%f%f%f', 2, 'headerLines', 3);       % read title
    dataGas =cell2mat(cellGas);
    TF.aCoeff_Fuel(2*i-1,1:5)=dataGas(1,1:5);
    TF.aCoeff_Fuel(2*i-1,6:7)=dataGas(2,1:2);
    TF.bCoeff_Fuel(2*i-1,1:2)=dataGas(2,4:5);
    cellGas = textscan(fid, '%f%f%f%f%f%f', 2, 'headerLines', 2);       % read title
    dataGas =cell2mat(cellGas);
    TF.aCoeff_Fuel(2*i,1:5)=dataGas(1,1:5);
    TF.aCoeff_Fuel(2*i,6:7)=dataGas(2,1:2);
    TF.bCoeff_Fuel(2*i,1:2)=dataGas(2,4:5);
end
fclose(fid);
Fcn_Load_JANAF_data
assignin('base','TF',TF); 
%
%--------------------------------------------------------------------------
% -- Kp for water-gas shift reaction
function Kp=Fcn_calculation_Kp_water_gas_shift_reaction(T)
global TF
% FileNameProdJANAF=  {'CO2','CO','H2O','H2','O2','N2'};
Ru= 8.3145;
n=40;
Delta_g=TF.Delta_g;
TJ=TF.TJ_Prod(1,1:n);
Delta_G_Water=Delta_g(1,1:n)-Delta_g(2,1:n)-Delta_g(3,1:n);
x=1000*Delta_G_Water(10:end);
pCoeff=polyfit(TJ(10:end),x,2);
y=polyval(pCoeff,T);
Kp=exp(-y./T./Ru);
%
%--------------------------------------------------------------------------
% -- Kp for CO2 dissociation
function Kp=Fcn_calculation_Kp_CO2_dissociation(T)
global TF
FileNameProdJANAF=  {'CO2','CO','H2O','H2','O2','N2'};% load combustion products from JANAF table
Ru= 8.3145;
n=40;
Delta_g=TF.Delta_g;
TJ=TF.TJ_Prod(1,1:n);
Delta_G=Delta_g(2,1:n)-Delta_g(1,1:n);
x=1000*Delta_G(10:end);
pCoeff=polyfit(TJ(10:end),x,2);
y=polyval(pCoeff,T);
Kp=exp(-y./T./Ru);
%
%--------------------------------------------------------------------------
% -- Kp for H2O dissociation
function Kp=Fcn_calculation_Kp_H2O_dissociation(T)
global TF
FileNameProdJANAF=  {'CO2','CO','H2O','H2','O2','N2'};% load combustion products from JANAF table
Ru= 8.3145;
n=40;
Delta_g=TF.Delta_g;
TJ=TF.TJ_Prod(1,1:n);
Delta_G=-Delta_g(3,1:n);
x=1000*Delta_G(10:end);
pCoeff=polyfit(TJ(10:end),x,2);
y=polyval(pCoeff,T);
Kp=exp(-y./T./Ru);
%
% Load the JANAF data of the main combustion products
function Fcn_Load_JANAF_data
global TF
FileNameProdJANAF=  {'CO2','CO','H2O','H2','O2','N2'};% load combustion products from JANAF table
Ru=TF.Ru;
cellGas=cell(6,1);
for i=1:6
    filename    = [FileNameProdJANAF{i} '.txt'];
    dataCas     = Fcn_Load_JANAF_data_From_txt(filename);
    dataSize    = size(dataCas);
    cellGas{i}  = dataCas;
    TF.TJ_Prod(i,1:dataSize(1))        = dataCas(1:end,1);
    TF.cpJ_Prod(i,1:dataSize(1))       = dataCas(1:end,2);
    TF.Delta_hsJ_Prod(i,1:dataSize(1)) = dataCas(1:end,5);
    TF.Delta_g(i,1:dataSize(1))        = dataCas(1:end,7);
end
assignin('base','TF',TF); 
%--------------------------------------------------------------------------
%-- Color set
function color_type=Fcn_color_set(color_number)
cmap = colormap(hot);
n_sample=round(linspace(1,64,color_number));
for ss=1:length(n_sample)
    color_type(ss,1:3)=cmap(n_sample(ss),1:3);
end
%
%--------------------------------------------------------------------------
%-- Load data from JANAF table
function dataGas=Fcn_Load_JANAF_data_From_txt(filename)     
fid=fopen(filename);
frewind(fid);                           % index the first row
cellGas = textscan(fid, '%f%f%f%f%f%f%f%f', inf, 'headerLines', 3);       % read title
dataGas=cell2mat(cellGas);
fclose(fid);