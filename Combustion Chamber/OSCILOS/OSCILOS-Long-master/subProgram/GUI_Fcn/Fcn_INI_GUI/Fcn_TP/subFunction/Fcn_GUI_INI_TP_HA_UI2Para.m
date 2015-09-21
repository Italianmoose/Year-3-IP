function Fcn_GUI_INI_TP_HA_UI2Para(hObject,HA_num)
% This function is used to update the parameters to calculate the heat
% addition for a selected heat addition interface
% Use the value of ui to update the value of different parameters
% first created: 2014-12-04
% last modified: 2015-06-03
% author: Jingxuan LI (jingxuan.li@imperial.ac.uk)
%
global CI
handles = guidata(hObject);
% ----------------get the index number of heat addition interface ---------
%
% HA_num   = get(handles.pop_HA_num,'value');
% HA_num is the index number of HA envisaged
%
% -----------------update the parameters for the calculation of HA --------
CI.TP.HA_style(HA_num)  = get(handles.pop_HA_style,'value');                    % heat addition style 
CI.TP.TRatio(HA_num)    = str2num(get(handles.edit_TD_Tratio,   'string'));     % temperature jump ratio
CI.TP.indexFuel(HA_num) = get(handles.pop_FD_fuel,'value');                     % index of fuel
CI.TP.eff(HA_num)       = str2num(get(handles.edit_FD_effi,     'string'));     % combustion efficiency
CI.TP.Phi(HA_num)       = str2num(get(handles.edit_FD_phi,      'string'));     % equivalence ratio
% CI.TP.dil(HA_num)       = str2num(get(handles.edit_FD_dilute,
% 'string'));     % diluted ratio, which is not used any more
%
assignin('base','CI',CI);
%
% -----------------------------end-----------------------------------------