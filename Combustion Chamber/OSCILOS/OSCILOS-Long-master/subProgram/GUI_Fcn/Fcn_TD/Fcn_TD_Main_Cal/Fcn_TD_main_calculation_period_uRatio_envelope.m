function Fcn_TD_main_calculation_period_uRatio_envelope(mm)
% calculate the envelope of a period
% author: Jingxuan LI (jingxuan.li@imperial.ac.uk)
% first created:    2014-11-18
% last edited:      2014-11-19
%
global CI
Var(1)     = CI.TD.nPadding + (CI.TD.indexPeriod(mm,1)-1)*CI.TD.nGap + 1;
Var(2)     = CI.TD.nPadding +  CI.TD.indexPeriod(mm,2)*CI.TD.nGap;
N = Var(2) - Var(1) + 1;
CI.TD.uRatioEnv(Var(1):Var(2))...
    = Fcn_calculation_envelope(CI.TD.uRatio(Var(1):Var(2)),N);
%
% ------------------------------end----------------------------------------
