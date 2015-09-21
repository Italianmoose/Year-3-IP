function tauPadding = Fcn_TD_Padding_time_delay
% This function is used to calculate the minimum padding time delay
% last edited: 2014-12-15
global CI
% --------------------------------
% Convection time of the acoustic waves in different sections
tauPadding      = max(  max(CI.TP.tau_plus),...
                        max(CI.TP.tau_minus));
tauPadding      = max(tauPadding,max(CI.TP.tau_c));
%
% --------------------------------
% Boundary condition
% Inlet boundary:
tauPadding = max(tauPadding, CI.TP.tau_minus(1) + CI.BC.tau_d1);
% Outlet boundary:
tauPadding = max(tauPadding, CI.TP.tau_plus(1)  + CI.BC.tau_d2);
%
% flames
% --------------------------------
if ~isempty(CI.CD.indexHP)    % if there are heat perturbations                   
    numHP                   = length(CI.CD.indexHP);                                % number of heat perturbations
    for ss = 1 : numHP
        HP = CI.FM.HP{ss};
        switch CI.FM.indexFM(ss)
            case 1
                taufNLimit      = 0;
                taufLimit       = HP.FTF.tauf + taufNLimit;
            case 2
                switch HP.NL.style
                    case {1,2}
                        taufNLimit      = 0;
                        taufLimit       = HP.FTF.tauf + taufNLimit;
                    case 3
                        % Flame describing function
                        % Herein, we only use the flame describing function model proposed by J.LI
                        % and A.Morgans JSV
                        % We must define a maximum velocity ratio value
                        uRatioLimit     = [0 1];
                        LfLimit         = interp1(  HP.NL.Model3.uRatio,...
                                                    HP.NL.Model3.Lf,...
                                                    uRatioLimit,'linear','extrap');
                        taufNLimit      = HP.NL.Model3.taufN.*(1 - LfLimit);
                        taufLimit       = HP.FTF.tauf + taufNLimit;
                    otherwise
                end
            case 3
                error('The current version does not support this situation!')
            case 4
                 % ************* need further change **************
                taufNLimit      = 0;  % G-equation
                taufLimit = HP.GEQU.tau_f + taufNLimit;
                 % ************* need further change **************s
        end
        tauPadding = max(tauPadding, max(taufLimit));
    end
end
% Green's function
tauPadding = max(tauPadding, max(CI.TD.Green.tauConv2));

%
% ------------------------ end --------------------------------------------
