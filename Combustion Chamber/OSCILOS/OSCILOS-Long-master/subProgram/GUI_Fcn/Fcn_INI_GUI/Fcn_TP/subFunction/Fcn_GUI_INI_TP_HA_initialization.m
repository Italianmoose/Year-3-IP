function Fcn_GUI_INI_TP_HA_initialization(varargin)
% This function is used to check the heat addition in the initialization
% stage and initialize some heat addition properties
% first created: 2014-12-03
% last modified: 2015-06-03
% author: Jingxuan LI (jingxuan.li@imperial.ac.uk)
%
global CI
hObject = varargin{1};
handles = guidata(hObject);
%
%----------------check if there is no heat addition -----------------------
%
Fcn_GUI_INI_TP_HA_visible(hObject)
%
%--------------- set the UI and initialize some heat addition properties -- 
%
if CI.CD.isHA == 0   
    set(handles.pop_HA_num,...
                        'string','No heat addition',...
                        'enable','off',...
                        'value',1);
    set(handles.edit_TD_Tratio,             'Enable', 'off',...
                                            'string', 1);
    set(handles.pop_HA_style,               'value',  1,...
                                            'enable', 'off');
    set(handles.ObjEditEnable_FD,           'Enable', 'off');
    %
    k = 1;                                       % still initialize these properties even there is no HA
    StringHA{k}         = ['Heat source number:' num2str(k)]; 
    CI.TP.indexHA(k)    = 0;                     % location of the heat addition interface
    CI.TP.HA_style(k)   = 1;                     % heat addition style   
    CI.TP.TRatio(k)     = 1;                     % temperature jump ratio
    CI.TP.indexFuel(k)  = 1;                     % index of fuel
    CI.TP.eff(k)        = 1;                     % combustion efficiency
    CI.TP.Phi(k)        = 1;                     % equivalence ratio
%     CI.TP.dil(k)        = 0;                     % diluted ratio
    CI.TP.Q(k)          = 0;                     % heat release rate   
    CI.TP.DeltaHr(k)    = 0;                     % heat release rate per mass flow rate 
    CI.TP.HA_num        = 1;                     % the flag is defaultly set to the 1 to indicate the first HA interface
else
    NumHA = length(CI.CD.indexHA);
    for k = 1:NumHA
        StringHA{k}         = ['Heat source number:' num2str(k)]; 
        CI.TP.indexHA(k)    = CI.CD.indexHA(k);      % location of the heat addition interface
        CI.TP.HA_style(k)   = 1;                     % heat addition style   
        CI.TP.TRatio(k)     = 1;                     % temperature jump ratio
        CI.TP.indexFuel(k)  = 1;                     % index of fuel
        CI.TP.eff(k)        = 1;                     % combustion efficiency
        CI.TP.Phi(k)        = 1;                     % equivalence ratio
%         CI.TP.dil(k)        = 0;                     % diluted ratio
        CI.TP.Q(k)          = 0;                     % heat release rate   
        CI.TP.DeltaHr(k)    = 0;                     % heat release rate per mass flow rate 
        CI.TP.HA_num        = 1;                     % the flag is defaultly set to the 1 to indicate the first HA interface
    end
    set(handles.pop_HA_num,                 'string',StringHA,...
                                            'visible','on',...
                                            'enable', 'on',...
                                            'value',   1);
    set(handles.pop_HA_style,               'value',   1,...
                                            'enable', 'on');
    set(handles.ObjEditEnable_FD,           'enable', 'off');
    set(handles.edit_TD_Tratio,             'enable', 'on',...
                                            'string', 1);
end
guidata(hObject, handles);
assignin('base','CI',CI);
%
% -----------------------------end-----------------------------------------