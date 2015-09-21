function Fcn_TD_INI_flame_model
% -------------------------------------------------------------------------
% ---------------------flame transfer function-----------------------------
% -------------------------------------------------------------------------
% For the first loop, the time delay and nonlinear gain ratio are set to
% the time delay of flame transfer function for weak perturbations and unit,
% respectively
%
% first created: 2014-11-11
% Last edit: 2014-12-15
%
global CI
uRatio0 = 0;                                                                % for weak velocity perturbations

if ~isempty(CI.CD.indexHP)    % if there are heat perturbations                   
    numHP                   = length(CI.CD.indexHP);                        % number of heat perturbations
    CI.TD.tauf      = zeros(numHP,CI.TD.nTotal);                             % Nonlinear time delay, which varies with velocity perturbations
    CI.TD.Lf        = zeros(numHP,CI.TD.nTotal);                               % Nonlinear model, which describes the saturation of heat release rate with velocity perturbations
    CI.TD.taufRem   = zeros(numHP,CI.TD.nTotal); 
    for ss = 1 : numHP
        HP = CI.FM.HP{ss};
        switch CI.FM.indexFM(ss)
            case 1 
                Lf(ss)      = 1;
                tauf(ss)    = HP.FTF.tauf;
            case 2
                switch HP.NL.style
                    case {1,2}
                        Lf(ss)      = 1;
                        tauf(ss)    = HP.FTF.tauf;
                    case 3
                        [Lf(ss),tauf(ss)] = Fcn_flame_model_NL_JLi_AMorgans(uRatio0, ss);
                    otherwise
                end
            case 3
                error('The current version does not support this situation!');
                % ************* need further change **************
            case 4 % G-equation case, these are filled for coherency with other models. 
                Lf(ss) = 1;
                tauf(ss) = HP.GEQU.tau_f;
                % ************* need further change **************
            otherwise
        end
        CI.TD.tauf(ss,:)    = CI.TD.tauf(ss,:) + tauf(ss);
        CI.TD.taufRem(ss,:) = CI.TD.tauf(ss,:) - CI.TD.taufMin(ss);
        CI.TD.Lf(ss,:)      = CI.TD.Lf(ss,:) + Lf(ss);
    end
end    
%
CI.TD.indexFlame = find(CI.CD.SectionIndex == 11);
%
assignin('base','CI',CI);
%
% --------------------------end--------------------------------------------