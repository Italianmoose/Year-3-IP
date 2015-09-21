function Fcn_TD_main_calculation(indexGUI)
% This function is used to lauch the time domain simulation
% The combustor dimensions, mean flow and thermal properties, flame model 
% and boundary conditions are initialized from OSCILOS. Then, user can lauch
% this file from OSCILOS_long.m file or manually lauch this file.
%
% check input 
if nargin == 0
    % default value
    indexGUI = 0;
end
switch indexGUI
    case 0
        Fcn_TD_main_calculation_without_GUI
    case 1
        Fcn_TD_main_calculation_with_GUI
end 
Fcn_TD_results_plot                                     % Plot the result
%
% -------------------------------------------------------------------------
% load from OSCILOS_long.m file
function Fcn_TD_main_calculation_with_GUI
global CI
switch CI.FM.indexFM(1) % at the moment all flame models must be the same for all flames                                                          
    case 3
        GUI_TD_Cal_JLI_AMorgans('OSCILOS_long', handles.figure);   
    otherwise
        Fcn_TD_main_calculation_without_identification_uRatio;
end
%
% -------------------------------------------------------------------------
% manually lauch, in this case, it must be lauched after the Green's function
% of the tranfer functions have been examined. 
function Fcn_TD_main_calculation_without_GUI
global CI
%
addpath(genpath('./'))                                                      % add directories to search path
Fcn_PreProcessing
dt                  = 1e-4;
tEndRaw             = 1;
uRatioMax           = 1;
nGapRatio           = 1;
NoiseInfo.level     = -40;                                                  % pNoise level 0.01 Pa
NoiseInfo.t(1)      = -100;
NoiseInfo.t(2)      = 100;   
%
Fcn_TD_INI(dt,tEndRaw,uRatioMax,nGapRatio,NoiseInfo)
%
nPeriod = 10;               % Time domain is separated to nPeriod sections.
RatioGapPadding = 1;        % Padding area ratio
Fcn_TD_INI_Period_Seperation(nPeriod,RatioGapPadding)
%
hWaitBar = waitbar(0,'Time domain calculations, please wait...');          
%
for mm = 1:CI.TD.nPeriod
    waitbar(mm/CI.TD.nPeriod);
    drawnow
    switch CI.FM.indexFM(1) % at the moment all flame models must be the same for all flames    
        case 3
            Iteration.tol = 1e-5;
            Iteration.num = 50;
            Fcn_TD_main_calculation_iteration_convergence(mm,Iteration)
            
        case 4 % G-equation
            CI.TD.IT = mm;
            % --------------------------
            Fcn_TD_calculation_one_period(mm)
            % --------------------------
            Fcn_TD_main_calculation_period_uRatio_envelope(mm)   % calculate the envelope
        otherwise
            % --------------------------
            Fcn_TD_calculation_one_period(mm)
            % --------------------------
            Fcn_TD_main_calculation_period_uRatio_envelope(mm)   % calculate the envelope
    end
end
close(hWaitBar)
assignin('base','CI',CI);
%
% --------------------------end--------------------------------------------




