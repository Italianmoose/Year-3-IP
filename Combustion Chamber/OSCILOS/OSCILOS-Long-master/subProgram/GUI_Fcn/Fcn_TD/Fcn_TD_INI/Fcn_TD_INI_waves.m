function Fcn_TD_INI_waves
% This function is used to initialize the waves
% Its parent program is "Fcn_TD_initialization".
% This function should be run after "Fcn_TD_initialization_samples_information"
%
% last edited: 2014-11-12 16:41 
global CI

CI.TD.AP                = zeros(CI.TP.numSection,CI.TD.nTotal);             % Downward acoustic waves
CI.TD.AM                = zeros(CI.TP.numSection,CI.TD.nTotal);             % Upward acoustic waves
CI.TD.E                 = zeros(CI.TP.numSection,CI.TD.nTotal);             % Entropy waves
%
if ~isempty(CI.CD.indexHP)    % if there are heat perturbations                   
    numHP                   = length(CI.CD.indexHP);                                % number of heat perturbations
    CI.TD.qRatio            = zeros(numHP,CI.TD.nTotal);                            % heat release ratio
    CI.TD.uRatio            = zeros(numHP,CI.TD.nTotal);                            % velocity ratio before the flame
    CI.TD.uRatioEnv         = zeros(numHP,CI.TD.nTotal);                            % envelope of velocity ratio, 
end
                                                                            % in the first step, it is set to zero
%
assignin('base','CI',CI);
%
% --------------------------end--------------------------------------------