function Fcn_test_nonlinear_hot_wire_model

uAsp = 0:1e-3:1;
for ss = 1:length(uAsp)
yAsp(ss) = Fcn_yA(uAsp(ss));
end

figure
plot(uAsp,yAsp,'-s')


function yA = Fcn_yA(uA)
fs = 1e4;
dt = 1/fs;
tSp = 0:dt:1;
% uA = 0.33;
f = 10;
tau = 1e-3;
omega = 2*pi*f;
u = uA.*sin(omega*(tSp-tau));
y = (abs(1/3 + u)).^0.5 - (1/3)^0.5;
yA = max(y) - min(y);

% figure
% plot(tSp,y)
