function Fcn_TD_calculation_one_period(indexPeriod)
% This function is used to calculate the wave values in one period
% one period means nGap*numGap number of time steps, which are calculated 
% in one father loop of Gap loop (Fcn_TD_calculation_one_gap), which is named 
% as Period loop (Fcn_TD_calculation_one_period)
%
% last edited: 2014-11-12 16:16
%
global CI
for nn = CI.TD.indexPeriod(indexPeriod,1):CI.TD.indexPeriod(indexPeriod,2)
    Var(1:2)     = CI.TD.nPadding + [(nn-1)*CI.TD.nGap + 1,...
                                         nn*CI.TD.nGap];

    % --------------------------
    Fcn_TD_calculation_one_gap(Var)
    % --------------------------
end
%
% ------------------------------end----------------------------------------

