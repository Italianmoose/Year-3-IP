function Fcn_TD_INI(dt,tEndRaw,uRatioMax,nGapRatio,NoiseInfo,ExtForceInfo)
global CI
% This function is used to initilize the time domain calculation code
% The sampling frequency and the raw time ends,
% and ...
% first created: 2014-11-12
% last edited: 2014-11-27
%
% -------------------------------------------------------------------------
% check input 
if nargin == 0
    % default value
    dt              = 1e-4;
    tEndRaw         = 2;
    uRatioMax       = 1;
    nGapRatio       = 1;
    NoiseInfo.level = -40;
    NoiseInfo.t     = [-100 100];
end
%
addpath(genpath('./'))                                                      % add directories to search path
% -------------------------------------------------------------------------
%
CI.TD.uRatioMax = uRatioMax;                                                % Set an uplimit of uRatio, which ensures that
%                                                                           % the time domain simulation results are not unreasonable. 
%
Fcn_TD_INI_matrix_preprocessing                                             % preprocessing some matrix connecting the components
%                                                                           % at two sides of an interface
%
Fcn_TD_INI_samples_information(dt,tEndRaw,uRatioMax,nGapRatio)              % initialize some samples
%       
Fcn_TD_INI_waves                                                            % initialize waves
%
Fcn_TD_INI_flame_model                                                      % some flame describing function model information
%
Fcn_TD_INI_background_noise(NoiseInfo)                                      % background noise information
%
Fcn_TD_INI_Green_function                                                   % Green's function information
% 
Fcn_TD_INI_constant_variables                                               % pre-calculate some constant values to improve the calculation speed
%   
% Fcn_TD_INI_Period_Seperation(nPeriod,RatioGapPadding)                     %
%
Fcn_TD_INI_external_forcing(ExtForceInfo)                                   % external forcing from a loudspeaker mounted to the combustor wall
%
assignin('base','CI',CI);
%
% ------------------------ end --------------------------------------------
