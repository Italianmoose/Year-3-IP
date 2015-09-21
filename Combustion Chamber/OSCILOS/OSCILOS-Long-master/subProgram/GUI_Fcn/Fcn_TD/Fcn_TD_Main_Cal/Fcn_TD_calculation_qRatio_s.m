function qRatio = Fcn_TD_calculation_qRatio_s(uRatio,FTFSys,tFTFSys,Var,tdRem,dt,Lf, indexHP)
% The flame transfer function is the function of s. 
% This function is used to calculate qRatio of flame
% uRatio is the velocity ratio before the flame
% The function Fcn_lsim_varied_td(y,Sys,tSys,Var,td,dt) is used to
% calculate qRatio
% Please refer to this function for detailed information
% index of time steps Var(1):Var(2)
% dt is the time step
% Lf is the nonlinear gain ratio
% time delay of flame transfer function is inclued in uRatio
%
% author: Jingxuan LI (jingxuan.li@imperial.ac.uk)
% first created:    2014-11-11
% last edited:      2014-12-15
%
global CI

qRatio = Fcn_lsim_varied_td(uRatio,FTFSys,tFTFSys,Var,tdRem,dt);

switch CI.FM.HP{1,indexHP}.NL.style
    case 2
        qRatio = Fcn_saturation(qRatio,CI.FM.HP{1,indexHP}.NL.Model2.alpha);
    otherwise
        qRatio = qRatio.*Lf((Var(1):Var(2)));
end
%
% --------------------------------end--------------------------------------
