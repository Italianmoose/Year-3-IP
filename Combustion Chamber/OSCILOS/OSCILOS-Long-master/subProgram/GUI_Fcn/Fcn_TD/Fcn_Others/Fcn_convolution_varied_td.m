function z = Fcn_convolution_varied_td(y,x,Var,tauf,dt)
% This function is used to calculate the convolution in case the time delay
% of the Green's function is changing with time
% x denotes the Green's function
% y denotes the array to be convoluted.
% tauf denotes the time delay of the transfer function
% z = conv(y,x(t-tauf)),
% which can also be changed to 
% z = conv(y(t-tauf),x);
% Discrete form:
% z(Var(1):Var(2)) = conv(y((Var(1):Var(2)) - nf), x)
% where, nf = tauf./dt;
% dt is the time step
% author: Jingxuan LI (jingxuan.li@imperial.ac.uk)
% first created:    2014-11-18
% last edited:      2014-11-19
% 
Nx          = length(x);
taufMin     = min(tauf);
taufMax     = max(tauf);
nfMax       = floor(taufMax./dt);
nfMin       = floor(taufMin./dt);
%
yindex      = Var(1) - (nfMax + 1) - (Nx - 1) : Var(2) - nfMin;
yA          = y(yindex);
index0      = Var(1) - (nfMax + 1) : Var(2) - nfMin;
indexIntp   = (Var(1):Var(2)) - tauf./dt;
%
zConv       = conv(yA,x,'valid').*dt;
z           = interp1(index0,zConv,indexIntp,'linear','extrap');
%
% ------------------------ end --------------------------------------------
