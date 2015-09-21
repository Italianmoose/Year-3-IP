function FcnKpCoeff

g(1,1) = 84;
g(1,2) = 2.79e5;
g(2,1) = 57.8;
g(2,2) = 2.51e5;
g(3,1) = g(2,1) - g(1,1);
g(3,2) = g(2,2) - g(1,2);
R = 8.3145;
load('THERMO.mat')
THERMO.gKp = g;
save('THERMO.mat','THERMO');




function Kp3 = Fcn_Kp3(T)
g11 = 84;
g12 = 2.79e5;
g21 = 57.8;
g22 = 2.51e5;
g31 = g21 - g11;
g32 = g22 - g12;
R = 8.3145;
Kp3 = exp((g31*T - g32)./R./T);

function Kp1 = Fcn_Kp1(T)
g11 = 84;
g12 = 2.79e5;
g21 = 57.8;
g22 = 2.51e5;
g31 = g21 - g11;
g32 = g22 - g12;
R = 8.3145;
Kp1 = exp((g11*T - g12)./R./T);