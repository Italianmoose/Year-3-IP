function z = Fcn_interp1_varied_td(y,Var,tau,dt)
% This function is used to calculate an array of values y(t-tau) based on
% 1-D interpolation
% td/dt are not integers
% interpolation is necessary
% dt is the time step
% tau can be an array
% if tau is an array, the length of tau should be equal to Var
% author: Jingxuan LI (jingxuan.li@imperial.ac.uk)
% first created:    2014-11-11
% last edited:      2014-11-19
%
tauMin      = min(tau);
tauMax      = max(tau);
nMin        = floor(tauMin/dt);
nMax        = floor(tauMax/dt);
%
index0      = Var(1) - (nMax + 1) : Var(2) - nMin;
indexIntp   = (Var(1):Var(2)) - tau./dt;
z0          = y(index0);
z           = interp1(index0,z0,indexIntp,'linear','extrap');
%
% ------------------------ end --------------------------------------------
