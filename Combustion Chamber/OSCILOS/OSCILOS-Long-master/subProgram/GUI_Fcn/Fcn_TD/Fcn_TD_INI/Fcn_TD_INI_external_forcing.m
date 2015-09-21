function Fcn_TD_INI_external_forcing(ExtForceInfo)
% -------------------------------------------------------------------------
% This function is used to set the external forcing information
% The external forcing is from a loudspeaker whose position is specified in
% the parameter configuration panel of OSCILOS. 
% For convient, we only account for the positions at the inlet and outlet,
% and the positions before and after the interfaces connecting neighbour
% modules. 
% Current forcing signal is a sinusoidal singal, which is defined as:
% p' = a*cos(2*pi*f*t) + b*sin(2*pi*f*t)
% where a and b are the magnitudes of the sine and cos waves,respectively.
% f is the frquency of the signal. 
% all the informations are saved in the cell `ExtForceInfo'
% It contains: 
% The position of the loudspeaker, use an index
% ExtForceInfo.indexPos = n
% n = 1: inlet
% n = end: outlet
% n = 2*n1: upstream side of the interface n1
% n = 2*n1+1: downstream side of the interface n1
% a, b, f are also saved
% the start time and end time are saved to t(1:2)
% ExtForceInfo.indexForcing indicates that if the external forcing is used.
% when ExtForceInfo.indexForcing == 1, without forcing
% when ExtForceInfo.indexForcing == 2, sinusoidal forcing
% ... for future development
%
% first created: 2015-01-22 
% last edited: 2015-01-22
%
% check input 
global CI
if nargin == 0
    % default value
    ExtForceInfo.indexForcing   = 0;    
    ExtForceInfo.indexPos       = 1;  
    ExtForceInfo.a = []; 
    ExtForceInfo.b = [];
    ExtForceInfo.f = []; 
    ExtForceInfo.t = zeros(1,2);
end
CI.TD.ExtForceInfo = ExtForceInfo;
switch ExtForceInfo.indexForcing
    case 1
        CI.TD.pExtForce = zeros(1,CI.TD.nTotal);
    case 2
        omega   = 2*pi*ExtForceInfo.f;
        N1      = round((ExtForceInfo.t(1)-CI.TD.tSpTotal(1)).*CI.TD.fs) + 1;
        N2      = round((ExtForceInfo.t(end)-CI.TD.tSpTotal(1)).*CI.TD.fs) + 1;
        CI.TD.pExtForce         = ExtForceInfo.a.*cos(omega.*CI.TD.tSpTotal) + ExtForceInfo.b.*sin(omega.*CI.TD.tSpTotal);
        if N1 >1
            CI.TD.pExtForce(1:N1)   = 0;
        end
        if N2 < CI.TD.nTotal
            CI.TD.pExtForce(N2:end) = 0;   
        end
end
%
assignin('base','CI',CI);
%
% ------------------------end----------------------------------------------