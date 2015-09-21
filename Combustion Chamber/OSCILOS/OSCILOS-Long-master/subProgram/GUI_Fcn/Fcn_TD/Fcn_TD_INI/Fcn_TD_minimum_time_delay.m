function [tauMin,taufMin] = Fcn_TD_minimum_time_delay(uRatioMax)
% This function is used to calculate the minimum time delay 
% In this version, the convection of entropy is not accounted for
% the function ``Fcn_PreProcessing '' must be run before this program
% last edit: 2014-12-15
%
global CI
% --------------------------------
% Convection time of the acoustic waves in different sections
tauMin = min(   min(CI.TP.tau_plus),...
                min(CI.TP.tau_minus));
%
% --------------------------------
% Boundary condition
% Inlet boundary:
tauMin = min(tauMin, CI.TP.tau_minus(1) + CI.BC.tau_d1);
% Outlet boundary:
tauMin = min(tauMin, CI.TP.tau_plus(1)  + CI.BC.tau_d2);
%
% --------------------------------
% flame transfer function
if ~isempty(CI.CD.indexHP)    % if there are heat perturbations                   
    numHP                   = length(CI.CD.indexHP);                                % number of heat perturbations
    for ss = 1 : numHP
        HP = CI.FM.HP{ss};
        switch CI.FM.indexFM(ss)
            case 1
                taufMin(ss) = HP.FTF.tauf;
            case 2
                switch HP.NL.style
                    case {1,2}
                        taufMin(ss) = HP.FTF.tauf;
                    case 3
                        % Flame describing function model proposed by J.LI
                        % and A.Morgans JSV
                        % We must define a maximum velocity ratio value
                        uRatioLimit     = [0 uRatioMax];
                        LfLimit         = interp1(  HP.NL.Model3.uRatio,...
                                                    HP.NL.Model3.Lf,...
                                                    uRatioLimit,'linear','extrap');
                        taufNLimit      = HP.NL.Model3.taufN.*(1 - LfLimit);
                        taufLimit       = HP.FTF.tauf + taufNLimit;
                        taufMin(ss)     = min(taufLimit);
                end
            case 3
                error('The current version does not support this situation!')
                
            case 4 % G-equation, the time delay here is fixed through GUI  
                taufMin(ss) = HP.GEQU.tau_f;
                
        end
        taufMin0 = min(taufMin);
    end
    tauMin = min(tauMin, taufMin0);
else  % in case without heat perturbation
    taufMin = [];
end
%
% ------------------------ end --------------------------------------------