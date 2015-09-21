function [qRatioLinear,Lfq] = Fcn_GUI_INI_FM_Nonlinear_model_Dowling(alpha,DqRatio,qRatioMax)
% Dowling's nonlinear model
qRatioLinear    = DqRatio:DqRatio:qRatioMax;
beta            = qRatioLinear./alpha;
Lfq             = zeros(1,length(beta));
for ss = 1:length(Lfq)
    if beta(ss) <= 1
        Lfq(ss) = 1;
    else
        psi = acos(1./beta(ss));
        Lfq(ss) = 1 - 2*psi./pi...
            + 2*(1 - 1./beta(ss).^2).^0.5./(pi*beta(ss));
    end
end
%
% ----------------------------end------------------------------------------

