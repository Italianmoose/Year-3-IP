function Fcn_test_time_delay_in_area_change_tube
% close all

dx = 1e-3;
c0 = 340;
r0 = 20e-3;
r1 = 40e-3;
l = 100e-3;
xSp = 0:dx:l;
rSp = r0 +(r1-r0).*xSp./l;
ASp = pi*rSp.^2;
M0 = 1e-1;
u0 = M0.*c0;
uSp = u0*ASp(1)./ASp;
N =length(xSp);
tauSp(1,2:N) = dx./(c0 + uSp(1:N-1));
tauSp(2,2:N) = dx./(c0 - uSp(1:N-1));
tauSp(1:2,1) = 0;
DeltatSp(1,:) = cumsum(tauSp(1,:));
DeltatSp(2,:) = cumsum(tauSp(2,:));

Deltat0 = l./c0;
alpha = (r1/r0 - 1)./l;
beta = (u0/c0)^0.5;
DeltatSpTHEO(1,:) = xSp./c0 -  u0./(alpha*c0.^2*beta)*(atan((1+alpha*xSp)./beta) - atan(1/beta)); 

DeltatSpTHEO(2,:) = xSp./c0 -  u0./(2*alpha*c0.^2*beta)*( log(abs((1+alpha*xSp+beta)./(1+alpha*xSp-beta)))...
                                                        -log(abs((1+beta)./(1-beta)))); 




figure
hold on
plot(xSp,DeltatSp(1,:),'-r')
plot(xSp,DeltatSpTHEO(1,:),'s','color','r')
plot(xSp,DeltatSp(2,:),'-k')
plot(xSp,DeltatSpTHEO(2,:),'s','color','k')
xlabel('x [m]')
ylabel('tau [s]')