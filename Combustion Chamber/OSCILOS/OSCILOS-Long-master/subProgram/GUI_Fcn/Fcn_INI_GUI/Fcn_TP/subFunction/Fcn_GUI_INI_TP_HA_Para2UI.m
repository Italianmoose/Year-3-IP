function Fcn_GUI_INI_TP_HA_Para2UI(varargin)
% This funcition is used to update the value of UI by the values of
% parameters
% first created: 2014-12-04
% last modified: 2015-06-03
% author: Jingxuan LI (jingxuan.li@imperial.ac.uk)
%
global CI
hObject     = varargin{1};
handles     = guidata(hObject);
% ----------------get the index number of heat addition interface ---------
%
HA_num      = get(handles.pop_HA_num,'value');  % get the order of current HA 
CI.TP.HA_num= HA_num ;                          % the flag is changed to current HA 
% HA_num is the index number of HA envisaged
%
% -----------------update the parameters for the calculation of HA --------
%
set(handles.pop_HA_style,       'value',    CI.TP.HA_style(HA_num));        % heat addition style 
set(handles.edit_TD_Tratio,     'string',   num2str(CI.TP.TRatio(HA_num))); % temperature jump ratio
set(handles.pop_FD_fuel,        'value',    CI.TP.indexFuel(HA_num));       % index of fuel
set(handles.edit_FD_effi,       'string',   num2str(CI.TP.eff(HA_num)));    % combustion efficiency
set(handles.edit_FD_phi,        'string',   num2str(CI.TP.Phi(HA_num)));    % equivalence ratio
% now, we do not use dilution ratio
% set(handles.edit_FD_dilute,     'string',   num2str(CI.TP.dil(HA_num)));    % diluted ratio
set(handles.edit_HR,            'string',   num2str(CI.TP.Q(HA_num)./1e3)); % heat release rate
%
% ----------------enable and disable some UI -----------------------------
%
pop_HA_style = CI.TP.HA_style(HA_num);
switch pop_HA_style
% {  'heater (Tb/Tu is given)';...
%    'combustion'}
    case 1
        set(handles.edit_TD_Tratio,             'Enable', 'on')
        set(handles.ObjEditEnable_FD,           'Enable', 'off')
    case 2
        set(handles.edit_TD_Tratio,             'Enable', 'off')
        set(handles.ObjEditEnable_FD,           'Enable', 'on')
end
%
guidata(hObject, handles);
assignin('base','CI',CI);
%
% -----------------------------end-----------------------------------------