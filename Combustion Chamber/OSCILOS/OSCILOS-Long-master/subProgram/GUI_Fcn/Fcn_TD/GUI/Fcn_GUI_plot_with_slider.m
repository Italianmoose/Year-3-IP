function Fcn_GUI_plot_with_slider(hObject,data,slider_vec,x_vec,titles)

% This function creates a plot of a matrix with a slider. 
% the data matrix should have the x_values as columns, and the slider
% values as lines
handles = guidata(hObject);

fig = handles.figure;
panel = handles.uipanel_axes;
axes = handles.axes1;
set(panel,...
    'units', 'points')
pannelsize=get(panel,'position');
pW=pannelsize(3);
pH=pannelsize(4);
% set(handles.axes1,      'units', 'points',...
%     'Fontunits','points',...
%     'position',[pW*2/10 pH*2.0/10 pW*7.5/10 pH*6.5/10],...
%     'fontsize',handles.FontSize(1),...
%     'color',handles.bgcolor{1},...
%     'box','on');

% Create plot
f = fig;
h = plot(axes,x_vec,data(1,:));

xlabel(titles(1))
ylabel(titles(2))

% add the slider
handles.slider = uicontrol('Parent',panel,'Style','slider','Position',[pW*2.15/10  pH*0.7/10 pW*9/10 pH*0.5/10],...
     'value',1, 'min',1, 'max',length(slider_vec));
bgcolor = get(f,'Color'); % colour of the figure
handles.slider_bl1 = uicontrol('Parent',panel,'Style','text','Position',[pW*0.4/10  pH*0.25/10 pW*1.5/10 pH*1/10],...
    'String',num2str(slider_vec(1)),'BackgroundColor',bgcolor);
handles.slider_bl2 = uicontrol('Parent',panel,'Style','text','Position',[pW*11.4/10  pH*0.25/10 pW*1.5/10 pH*1/10],...
    'String',num2str(slider_vec(end)),'BackgroundColor',bgcolor);
handles.slider_bl3 = uicontrol('Parent',panel,'Style','text','Position',[pW*2/10  pH*0.05/10 pW*10/10 pH*0.5/10],...
    'String',strcat(titles(3), ' = ' ,num2str(slider_vec(1))),'BackgroundColor',bgcolor);

% add a listener for the slider handle b
try    % R2013b and older, if it is still valid and you get a warning to use ContinuousValueChange, disregard warning
    addlistener(handles.slider,'ActionEvent',@(src,eventdata)slider_callback(src,eventdata,data,slider_vec,h,panel,bgcolor,titles(3),pW,pH));
catch  % R2014a and newer
    addlistener(handles.slider,'ContinuousValueChange',@(src,eventdata)slider_callback(src,eventdata,data,slider_vec,h,panel,bgcolor,titles(3),pW,pH));
end

guidata(hObject, handles);

end

function slider_callback(src,evt,data,slider_vec,plot_handle,panel_handle,bgcolor,slider_title,pW,pH)
slider_value=round(get(src,'value'));
set(plot_handle,'YData',data(slider_value,:))
uicontrol('Parent',panel_handle,'Style','text','Position',[pW*2/10  pH*0.05/10 pW*10/10 pH*0.5/10],...
    'String',strcat(slider_title, ' = ' ,num2str(slider_vec(slider_value))),'BackgroundColor',bgcolor);

end