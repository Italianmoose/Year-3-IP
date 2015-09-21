function handles = Fcn_GUI_default_configuration(handles)
% This function is used to set the default GUI parameters
handles.bgcolor{1} = [0.95, 0.95, 0.95];
handles.bgcolor{2} = [0, 0, 0];
handles.bgcolor{3} = [.75, .75, .75];
handles.bgcolor{4} = [0.90,0.90,1];
%
handles.sW  = 800;
handles.sH  = 600;
%
if ispc
    handles.FontSize(1)=11;                 % set the default fontsize
    handles.FontSize(2)=9;
else
    handles.FontSize(1)=12;                 % set the default fontsize
    handles.FontSize(2)=10;   
end