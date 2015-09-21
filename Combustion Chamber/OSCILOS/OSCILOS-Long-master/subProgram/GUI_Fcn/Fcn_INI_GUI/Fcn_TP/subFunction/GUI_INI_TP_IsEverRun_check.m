function GUI_INI_TP_IsEverRun_check(varargin)
% This function is used to check if the function GUI_INI_TP has ever been
% run
% first created: 2014-12-03
% last modified: 2014-12-04
% author: Jingxuan LI (jingxuan.li@imperial.ac.uk)
hObject = varargin{1};
handles = guidata(hObject);
global CI
switch CI.IsRun.GUI_INI_TP    % This value is defined in the OSCILOS_Long.m
    case 1
        % ---------------------enable pushbutton---------------------------
        %
        set(handles.pb_Cal,             'enable','on');
        set(handles.pb_Plot,            'enable','on');
        set(handles.pb_saveFig,         'enable','on');
        set(handles.pb_OK,              'enable','on');
        set(handles.pop_plot,           'enable','on')
        %
        % --------------------inlet mean properties -----------------------
        %
        set(handles.edit_TP_p1,         'string',   num2str(CI.TP.p_mean(1,1)));
        set(handles.edit_TP_T1,         'string',   num2str(CI.TP.T_mean(1,1)));
        set(handles.pop_TP_M1_u1,       'value',    2);                             % set to mean velocity
        set(handles.edit_TP_M1,         'string',   num2str(CI.TP.u_mean(1,1)));    % ...
        try
            set(handles.pop_gamma,          'Value',    CI.TP.index_gamma);
        catch
            set(handles.pop_gamma,          'Value',    1);
        end
        %
        % --------------------update the rest -----------------------------
        %--------------- set the UI and initialize some heat addition properties -- 
        if CI.CD.isHA == 1   
            NumHA = length(CI.CD.indexHA);
            for k = 1:NumHA
                StringHA{k}         = ['Heat source number:' num2str(k)]; 
            end
        end
        try
        set(handles.pop_HA_num,...
                        'string',StringHA,...
                        'enable','on',...
                        'value',1);   
        catch
        end
        % -----------------------------------------------------------------
        guidata(hObject, handles);
        Fcn_GUI_INI_TP_HA_Para2UI(hObject)
        handles = guidata(hObject);
        guidata(hObject, handles);
        Fcn_GUI_INI_TP_HA_visible(hObject)
        handles = guidata(hObject);
        guidata(hObject, handles);
        Fcn_GUI_INI_TP_plot(handles.axes1,handles);
        handles = guidata(hObject);
        guidata(hObject, handles);

    case 0
        % ------------------initialize HA -----------------------------------------
        Fcn_GUI_INI_TP_HA_initialization(hObject);
        handles = guidata(hObject);
        guidata(hObject, handles);
end
%
% ------------------------------end----------------------------------------