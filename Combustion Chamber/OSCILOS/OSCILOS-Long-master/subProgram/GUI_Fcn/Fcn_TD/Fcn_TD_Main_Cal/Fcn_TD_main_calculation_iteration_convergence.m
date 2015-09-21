function Fcn_TD_main_calculation_iteration_convergence(mm,Iteration)
% This function is used to calculate the true envelope of uRatio by
% iterative method.
% The criteria is: 
% abs(std(uRatioEnv_pre - uRatioEnv_now)) < Tolerance. 
% The nonlinear flame transfer function is updated in every loop, based on 
% the new envelope of the velocity ratio.
% Since Hilbert transform is used to calculate the envelope, head and end
% padding signals are necessary. Moreover, zero head and end padding signals
% of the padded signal are necessary. 
% 
% In every period calculation, the nonlinear flame describing functions of
% current period and end padding period are updated, however, those of 
% previous period are not updated.
%
% Iteration is a field which contains:
% tol : tolerance
% num : maximum iteration times
%
% author: Jingxuan LI (jingxuan.li@imperial.ac.uk)
% first created:    2014-11-18
% last edited:      2014-11-19
%
global CI
%
Var(1)      = CI.TD.nPadding + (CI.TD.indexPeriod(mm,1)-1)*CI.TD.nGap + 1;   % Start of the period
Var(2)      = CI.TD.nPadding +  CI.TD.indexPeriod(mm,2)*CI.TD.nGap;          % End of the period                                           % Length of the period
N1          = CI.TD.numGap.*CI.TD.nGap;
% We only focus on the envelope from samples Var(1) to Var(1)-1+N1, 
% the rest from Var(1) + N1 to Var(2) is the padding signal which is used
% to avoid errors in Hilbert Transform calculation
% 
%
VarPad(1)   = CI.TD.nPadding + 1;   % Start of the signal to be processed by Hilbert transform
VarPad(2)   = Var(2);               % End of the signal to ...
N           = VarPad(2) - VarPad(1) + 1;
% ------------------------------------------
hFig        = figure;
hAxes       = axes;
titleText = ['Finished:' num2str(mm) '/' num2str(CI.TD.nPeriod) ];
title(hAxes,titleText)
set(hAxes,'yscale','log')

hFig1 = figure;
hAxes1 = axes;

% ------------------------------------
uRatioEnv2  = zeros(1,N);
for ss = 1 : Iteration.num
    % --------------------------
    Fcn_TD_calculation_one_period(mm)
    % --------------------------
    uRatioEnv1 = Fcn_calculation_envelope(CI.TD.uRatio(VarPad(1):VarPad(2)),N);       % calculate uRatio
    uRatioEnv1 = 0.5*(uRatioEnv1 + uRatioEnv2);
    
    [   CI.TD.Lf(Var(1):VarPad(2)),...
        CI.TD.tauf(Var(1):VarPad(2))]...
        = Fcn_flame_model_NL_JLi_AMorgans(uRatioEnv1(Var(1)-CI.TD.nPadding:end));       % only update the nonlinear value from Var(1) to VarPad(2)

        CI.TD.taufRem(Var(1):VarPad(2))  = CI.TD.tauf(Var(1):VarPad(2)) - CI.TD.taufMin;
    index = [Var(1), Var(1)-1+N1] - CI.TD.nPadding;
    DeltaDiff = abs(std(uRatioEnv1(index(1) : index(2)) - uRatioEnv2(index(1) : index(2))))   
    hold on
    plot(hAxes,ss,DeltaDiff,'s','markersize',10)
    drawnow
    if DeltaDiff < Iteration.tol
        break;
    end
    uRatioEnv2 = uRatioEnv1;
    %----------------
    hPlot = plot(hAxes1,CI.TD.tSpTotal,CI.TD.uRatio,'-k');
    drawnow
    delete(hPlot)
    % ----------------------------
end
close(hFig)
close(hFig1)
CI.TD.uRatioEnv(Var(1):Var(2)) = uRatioEnv2((Var(1):Var(2))-CI.TD.nPadding);




