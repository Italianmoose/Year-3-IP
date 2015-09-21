function Fcn_TD_INI_Green_function
% -------------------------------------------------------------------------
% -------------------------Green's function--------------------------------
% -------------------------------------------------------------------------
% This function is used to recalculate the Green's function based on the
% new sampling rate
% index meaning:
% 1. inlet boundary condition
% 2. outlet boundary condition
% 3. flame transfer function for weak perturbations
% The cut-off time delays are from the GUI function "GUI_TD_GreenFcn"
% last edited: 2014-11-14 17:26
%
global CI
for ss = 1:length(CI.TD.Green.tauConv2)
    if CI.TD.Green.tauConv2(ss) == 0
        CI.TD.Green.y{ss} = CI.TD.Green.y2{ss}.*CI.TD.fs;
        CI.TD.Green.t{ss} = 0;
    else
        CI.TD.Green.tauConv2(ss) = ceil(CI.TD.Green.tauConv2(ss).*CI.TD.fs)./CI.TD.fs;
        CI.TD.Green.t{ss} = 0:CI.TD.dt:CI.TD.Green.tauConv2(ss);
        CI.TD.Green.y{ss}...
        = impulse(CI.TD.Green.sys{ss},CI.TD.Green.t{ss});
    end
end
%
assignin('base','CI',CI);
%
% ------------------------ end --------------------------------------------