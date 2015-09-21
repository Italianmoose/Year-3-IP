function [uRatio,taufRem] = Fcn_TD_calculation_uRatio(Var,indexSection, indexHP)
% This function is used to calculate the velocity ratio before the flame,
% uRatio = (AP(t-tauP-tauf) - AM(t-tauf))/(rho*c*u);
% The time delay of the flame describing function is included.
% index of time steps Var(1):Var(2)
% indexSection denotes the index of section before the flame
% author: Jingxuan LI (jingxuan.li@imperial.ac.uk)
% first created:    2014-11-11
% last edited:      2014-11-19
global CI
ss  = indexSection;
% herein, when length(Var) > 1, CI.TD.taufMin should change with time, but
% to increase the calculation speed, we assume it is a constant value
AP  = Fcn_interp1_varied_td(    CI.TD.AP(ss,:),...
                                Var,...
                                CI.TD.taufMin(indexHP) + CI.TP.tau_plus(ss),...
                                CI.TD.dt);
AM  = Fcn_interp1_varied_td(    CI.TD.AM(ss,:),...
                                Var,...
                                CI.TD.taufMin(indexHP),...
                                CI.TD.dt);
uRatio = (AP - AM)./CI.TP.RhoCU(ss);
if CI.FM.indexFM(indexHP) < 3
% if CI.FM.NL.style ~= 4
    % CI.TP.RhoCU denotes rho*c*u
    % velocity ratio limit : CI.TD.uRatioMax, which is first defined in the
    % program : Fcn_TD_initialization
    uRatio = Fcn_saturation(uRatio,CI.TD.uRatioMax);
end
%
% ----------------------------end------------------------------------------