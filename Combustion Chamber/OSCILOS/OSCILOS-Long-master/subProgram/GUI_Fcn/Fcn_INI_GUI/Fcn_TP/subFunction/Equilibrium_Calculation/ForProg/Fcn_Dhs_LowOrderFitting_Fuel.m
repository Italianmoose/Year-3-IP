function FIT = Fcn_Dhs_LowOrderFitting_Fuel
close all
%  Fcn_Dhs_LowOrderFitting -- low order fitting of Dhs at a
% temperature range
% T1: lower temperature
% T2: upper temperature
% N: fitting order 
%
% if nargin < 3
%    N = 2;
%    %
%    disp('fitting order = 2')
% end
% if nargin < 2
%    T2 = 1000;
%    disp('T2 = 3000 K')
% end
% if nargin < 1
%    T1 = 300;
%    disp('T1 = 1000 K')
% end
N = 2;
load 'GUI_Data.mat'
TRange = [300 1000 3000];
indexRange = 1;
TSp1 = linspace(TRange(indexRange), TRange(indexRange+1), 100);
TSp2 = linspace(TRange(indexRange+1), TRange(indexRange+2), 200);
Dhs1 = zeros(6, length(TSp1));
for kk = 1 : 8
    for ss = 1: length(TSp1)
        Dhs1(kk,ss) =  Fcn_Cal_Deltahs(TSp1(ss),2,kk,TF);
    end
    p1(kk,:) = polyfit(TSp1, Dhs1(kk,:),N);
    DhsFit1(kk,:) = polyval(p1(kk,:), TSp1);
    error1(kk) = mean( abs(Dhs1(kk,:) - DhsFit1(kk,:))./Dhs1(kk,:));
end
for kk = 1 : 8
    for ss = 1: length(TSp2)
        Dhs2(kk,ss) =  Fcn_Cal_Deltahs(TSp2(ss),2,kk,TF);
    end
    p2(kk,:) = polyfit(TSp2, Dhs2(kk,:),N);
    DhsFit2(kk,:) = polyval(p2(kk,:), TSp2);
    error2(kk) = mean( abs(Dhs2(kk,:) - DhsFit2(kk,:))./Dhs2(kk,:));
end
FIT.gas = {'CH4', 'C2H4', 'C2H6', 'C3H8', '1-butene', 'n-butane', 'isobutane', 'Jet-A'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('THERMO.mat')
THERMO.Fuel.LowOrder.pDhs       = [p1,p2];
THERMO.Fuel.LowOrder.errorDhs   = [error1',error2'];
THERMO.Fuel.LowOrder.OrderDhs   = N;
save('THERMO.mat','THERMO');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nStep = 5;
if nargout < 1
    h = figure;
    axes
    hold on
    title('Dhs with temperature')
    TSpSum = [TSp1,TSp2];
    DhsFit = [DhsFit1, DhsFit2];
    for kk = 1 : 8
        plot([TSp1(1:nStep:end),TSp2(1:nStep:end)], [Dhs1(kk,1:nStep:end),Dhs2(kk,1:nStep:end)], '+r')
        plot(TSp1, DhsFit1(kk,:), '-k')
        plot(TSp2, DhsFit2(kk,:), '-k')
        index = ceil(kk*length(TSpSum)./9);
        text(TSpSum(index), DhsFit(kk,index), FIT.gas{kk},'backgroundcolor','w');
    end
    xlabel('Temperature [K]')
    ylabel('Dhs [J/(mol)]')
end
%
% -------------------------------end---------------------------------------