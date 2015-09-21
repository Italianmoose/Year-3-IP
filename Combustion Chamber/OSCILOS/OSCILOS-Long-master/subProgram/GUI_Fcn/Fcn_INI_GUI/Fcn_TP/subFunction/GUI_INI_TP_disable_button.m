function GUI_INI_TP_disable_button(varargin)
hObject = varargin{1};
handles = guidata(hObject);
set(handles.pb_Plot,            'enable','off');
set(handles.pb_saveFig,         'enable','off');
set(handles.pb_OK,              'enable','off');
guidata(hObject, handles);
% ----------------------------end------------------------------------------
