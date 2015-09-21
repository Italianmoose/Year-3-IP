function varargout = GUI_FREQ_EigCal(varargin)
% GUI_FREQ_EIGCAL MATLAB code for GUI_FREQ_EigCal.fig
%      GUI_FREQ_EIGCAL, by itself, creates a new GUI_FREQ_EIGCAL or raises the existing
%      singleton*.
%
%      H = GUI_FREQ_EIGCAL returns the handle to a new GUI_FREQ_EIGCAL or the handle to
%      the existing singleton*.
%
%      GUI_FREQ_EIGCAL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_FREQ_EIGCAL.M with the given input arguments.
%
%      GUI_FREQ_EIGCAL('Property','Value',...) creates a new GUI_FREQ_EIGCAL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_FREQ_EigCal_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_FREQ_EigCal_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_FREQ_EigCal

% Last Modified by GUIDE v2.5 02-Apr-2015 17:42:48

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_FREQ_EigCal_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_FREQ_EigCal_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before GUI_FREQ_EigCal is made visible.
function GUI_FREQ_EigCal_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_FREQ_EigCal (see VARARGIN)
% --------------------------------------------------------------------------
global CI
handles.indexEdit = 1;
switch handles.indexEdit 
    case 0
        %--------------------------------------------------------------------------
        dontOpen = false;
        mainGuiInput = find(strcmp(varargin, 'OSCILOS_long'));
        if (isempty(mainGuiInput)) ...
            || (length(varargin) <= mainGuiInput) ...
            || (~ishandle(varargin{mainGuiInput+1}))
            dontOpen = true;
        else
            % load from the main GUI
            % handles of main GUI
            handles.MainGUI = varargin{mainGuiInput+1};
            try
                handles.ExampleGUI = varargin{mainGuiInput+2};
            catch
            end
            % Obtain handles using GUIDATA with the caller's handle 
            mainHandles = guidata(handles.MainGUI);
            % background colors
            handles.bgcolor=mainHandles.bgcolor;
            % fontsize
            handles.FontSize=mainHandles.FontSize;
            %
            handles.sW = mainHandles.sW;
            handles.sH = mainHandles.sH;
            % Update handles structure
            guidata(hObject, handles);
            % Initialization
            GUI_FREQ_EigCal_Initialization(hObject, eventdata, handles)
            handles = guidata(hObject);
        end
        guidata(hObject, handles);
        handles.output = hObject;
        guidata(hObject, handles);
        if dontOpen
           disp('-----------------------------------------------------');
           disp('This program cannot be currently lauched due to some reasons.') 
           disp('-----------------------------------------------------');
        else
%            uiwait(hObject);
        end
    case 1
        handles = Fcn_GUI_default_configuration(handles);
        guidata(hObject, handles);  
        GUI_FREQ_EigCal_Initialization(hObject, eventdata, handles)
        handles = guidata(hObject);
        guidata(hObject, handles);
        handles.output = hObject;
        guidata(hObject, handles);
end
%
% -------------------------------------------------------------------------
%
function GUI_FREQ_EigCal_Initialization(varargin)
hObject = varargin{1};
handles = guidata(hObject);
global CI
set(0, 'units', 'points');
screenSize  = get(0, 'ScreenSize');                 % get the screen size
sW          = handles.sW;                           % screen width
sH          = handles.sH;                          % screen height
FigW=sW.*1;                                        % window width
FigH=sH.*3/4;                                        % window height
set(handles.figure,     'units', 'points',...
                        'position',[(screenSize(3)-FigW)./2 (screenSize(4)-FigH)./2 FigW FigH],...
                        'name','Eigenvalue calculation',...
                        'color',handles.bgcolor{3});
%----------------------------------------
% pannel axes
set(handles.uipanel_Axes,...
                        'units', 'points',...
                        'fontUnits','points',...
                        'position',[FigW*9.5/20 FigH*2.25/20 FigW*10/20 FigH*17.25/20],...
                        'Title','Plots',...
                        'visible','on',...
                        'highlightcolor',handles.bgcolor{3},...
                        'borderwidth',1,...
                        'fontsize',handles.FontSize(2),...
                        'backgroundcolor',handles.bgcolor{3});  
pannelsize=get(handles.uipanel_Axes,'position');
pW=pannelsize(3);
pH=pannelsize(4);                
set(handles.axes1,...
                        'units', 'points',...
                        'position',[pW*2/10 pH*4.75/10 pW*6.0/10 pH*3.25/10],...
                        'fontsize',handles.FontSize(1),...
                        'color',handles.bgcolor{1},...
                        'box','on');  
set(handles.axes2,...
                        'units', 'points',...
                        'position',[pW*2/10 pH*1.5/10 pW*6.0/10 pH*3.25/10],...
                        'fontsize',handles.FontSize(1),...
                        'color',handles.bgcolor{1},...
                        'box','on');  
set(handles.pop_numMode,...
                        'units', 'points',...
                        'fontUnits','points',...
                        'position',[pW*5.5/10 pH*8.5/10 pW*4/10 pH*1.0/10],...
                        'fontsize',handles.FontSize(2),...
                        'string',{'Mode number:1';'Mode number:2'},...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','left',...
                        'enable','on',...
                        'visible','off');  
set(handles.pop_PlotType,...
                        'units', 'points',...
                        'fontUnits','points',...
                        'position',[pW*0.5/10 pH*8.5/10 pW*4/10 pH*1.0/10],...
                        'fontsize',handles.FontSize(2),...
                        'string',{'Map of eigenvalues';'Modeshape'; 'Evolution of eigenvalue with velocity ratio'},...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','left',...
                        'enable','off'); 
guidata(hObject, handles);
%----------------------------------------
% pannel initialize velocity range information 
set(handles.uipanel_INI,...
                        'units', 'points',...
                        'fontUnits','points',...
                        'position',[FigW*0.5/20 FigH*12.5/20 FigW*8.5/20 FigH*7/20],...
                        'Title',{'Set perturbation ratio (v''/v before the flame or p''1/p1) range'},...
                        'visible','on',...
                        'highlightcolor',handles.bgcolor{3},...
                        'borderwidth',1,...
                        'fontsize',handles.FontSize(2),...
                        'backgroundcolor',handles.bgcolor{3});  
pannelsize = get(handles.uipanel_INI,'position');
pW=pannelsize(3);
pH=pannelsize(4);  
set(handles.text_uRatio_min,...
                        'units', 'points',...
                        'fontUnits','points',...
                        'position',[pW*0.5/10 pH*7.0/10 pW*6/10 pH*1.0/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Minimum value (u''/u_bar): [-]',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'horizontalalignment','left');                         
set(handles.edit_uRatio_min,...
                        'units', 'points',...
                        'fontUnits','points',...
                        'position',[pW*7/10 pH*7/10 pW*2.5/10 pH*1.25/10],...
                        'fontsize',handles.FontSize(2),...
                        'string',0,...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','right',...
                        'Enable','on');
set(handles.text_uRatio_max,...
                        'units', 'points',...
                        'fontUnits','points',...
                        'position',[pW*0.5/10 pH*5.0/10 pW*6/10 pH*1.0/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Maximum value (u''/u_bar): [-]',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'horizontalalignment','left');                         
set(handles.edit_uRatio_max,...
                        'units', 'points',...
                        'fontUnits','points',...
                        'position',[pW*7/10 pH*5/10 pW*2.5/10 pH*1.25/10],...
                        'fontsize',handles.FontSize(2),...
                        'string',1,...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','right',...
                        'Enable','on');                    
set(handles.text_uRatio_SampNum,...
                        'units', 'points',...
                        'fontUnits','points',...
                        'position',[pW*0.5/10 pH*3.0/10 pW*6/10 pH*1.0/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Number of velocity ratio samples: [-]',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'horizontalalignment','left');                         
set(handles.edit_uRatio_SampNum,...
                        'units', 'points',...
                        'fontUnits','points',...
                        'position',[pW*7/10 pH*3/10 pW*2.5/10 pH*1.25/10],...
                        'fontsize',handles.FontSize(2),...
                        'string',11,...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','right',...
                        'Enable','on');
                    
% Similar settings for pRatio
set(handles.text_pRatio_min,...
                        'units', 'points',...
                        'fontUnits','points',...
                        'position',[pW*0.5/10 pH*7.0/10 pW*6/10 pH*1.0/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Minimum value (p''_1/p_1): [x10^(-5)]',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'horizontalalignment','left');                         
set(handles.edit_pRatio_min,...
                        'units', 'points',...
                        'fontUnits','points',...
                        'position',[pW*7/10 pH*7/10 pW*2.5/10 pH*1.25/10],...
                        'fontsize',handles.FontSize(2),...
                        'string',0,...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','right',...
                        'Enable','on');
set(handles.text_pRatio_max,...
                        'units', 'points',...
                        'fontUnits','points',...
                        'position',[pW*0.5/10 pH*5.0/10 pW*6/10 pH*1.0/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Maximum value (p''_1/p_1): [x10^(-5)]',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'horizontalalignment','left');                         
set(handles.edit_pRatio_max,...
                        'units', 'points',...
                        'fontUnits','points',...
                        'position',[pW*7/10 pH*5/10 pW*2.5/10 pH*1.25/10],...
                        'fontsize',handles.FontSize(2),...
                        'string',1000,...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','right',...
                        'Enable','on');                    
set(handles.text_pRatio_SampNum,...
                        'units', 'points',...
                        'fontUnits','points',...
                        'position',[pW*0.5/10 pH*3.0/10 pW*6/10 pH*1.0/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Number of pressure ratio samples: [-]',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'horizontalalignment','left');                         
set(handles.edit_pRatio_SampNum,...
                        'units', 'points',...
                        'fontUnits','points',...
                        'position',[pW*7/10 pH*3/10 pW*2.5/10 pH*1.25/10],...
                        'fontsize',handles.FontSize(2),...
                        'string',11,...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','right',...
                        'Enable','on');                    
% Settings for the calculating botton and the scan domain setting botton                    
                                       
set(handles.pb_CalEig,...
                        'units', 'points',...
                        'fontUnits','points',...
                        'position',[pW*5.5/10 pH*0.5/10 pW*4/10 pH*1.5/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Calculate eigenvalues',...
                        'backgroundcolor',handles.bgcolor{3});
set(handles.pb_SetScanDomain,...
                        'units', 'points',...
                        'fontUnits','points',...
                        'position',[pW*0.5/10 pH*0.5/10 pW*4/10 pH*1.5/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Set scan range',...
                        'backgroundcolor',handles.bgcolor{3});
%----------------------------------------
% pannel 2
set(handles.uipanel_CAL,...
                        'units', 'points',...
                        'position',[FigW*0.5/20 FigH*2.25/20 FigW*8.5/20 FigH*9.75/20],...
                        'Title','',...
                        'visible','on',...
                        'highlightcolor',handles.bgcolor{3},...
                        'borderwidth',1,...
                        'fontsize',handles.FontSize(2),...
                        'backgroundcolor',handles.bgcolor{3});  
pannelsize = get(handles.uipanel_CAL,'position');
pW=pannelsize(3);
pH=pannelsize(4);  
cnames = {'Eigen-frequency [Hz]','Growth rate [rad/s]'};
%
set(handles.text_Table,...
                        'units', 'points',...
                        'fontUnits','points',...
                        'position',[pW*0.5/10 pH*9.25/10 pW*9.0/10 pH*0.5/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Eigenvalues for selected velocity ratio',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'horizontalalignment','left'); 
                    
set(handles.uitable,    'units', 'points',...
                        'fontUnits','points',...
                        'position',[pW*0.5/10 pH*4/10 pW*9/10 pH*4.5/10],...
                        'fontsize',handles.FontSize(1),...
                        'columnName',cnames,...
                        'data',[],...
                        'backgroundcolor',handles.bgcolor{1},...
                        'RearrangeableColumns','on'); 
set(handles.text_uRatio,...
                        'units', 'points',...
                        'fontUnits','points',...
                        'position',[pW*0.5/10 pH*2.25/10 pW*6/10 pH*0.9/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Selected velocity ratio: [-]',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'horizontalalignment','left');                         
set(handles.edit_uRatio,...
                        'units', 'points',...
                        'fontUnits','points',...
                        'position',[pW*7/10 pH*2.5/10 pW*2.5/10 pH*0.9/10],...
                        'fontsize',handles.FontSize(2),...
                        'string',0,...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','right',...
                        'Enable','on');                    
set(handles.slider_uRatio,...
                        'units', 'points',...
                        'fontUnits','points',...
                        'position',[pW*0.5/10 pH*0.5/10 pW*9/10 pH*1.5/10],...
                        'fontsize',handles.FontSize(2),...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','right',...
                        'Enable','on',...
                        'min',1,...
                        'max',str2double(get(handles.edit_uRatio_SampNum,'string')),...
                        'value',1);
% Similiar settings for pRatio
set(handles.text_pRatio,...
                        'units', 'points',...
                        'fontUnits','points',...
                        'position',[pW*0.5/10 pH*2.25/10 pW*6/10 pH*0.9/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Selected pressure ratio: [-]',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'horizontalalignment','left');                         
set(handles.edit_pRatio,...
                        'units', 'points',...
                        'fontUnits','points',...
                        'position',[pW*7/10 pH*2.5/10 pW*2.5/10 pH*0.9/10],...
                        'fontsize',handles.FontSize(2),...
                        'string',0,...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','right',...
                        'Enable','on');                    
set(handles.slider_pRatio,...
                        'units', 'points',...
                        'fontUnits','points',...
                        'position',[pW*0.5/10 pH*0.5/10 pW*9/10 pH*1.5/10],...
                        'fontsize',handles.FontSize(2),...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','right',...
                        'Enable','on',...
                        'min',1,...
                        'max',str2double(get(handles.edit_pRatio_SampNum,'string')),...
                        'value',1);




%---------------------------
% pannel AOC                   
set(handles.uipanel_AOC,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[FigW*0.5/20 FigH*0/20 FigW*19/20 FigH*1.75/20],...
                        'Title','',...
                        'visible','on',...
                        'highlightcolor',handles.bgcolor{3},...
                        'borderwidth',1,...
                        'fontsize',handles.FontSize(2),...
                        'backgroundcolor',handles.bgcolor{3}); 
pannelsize=get(handles.uipanel_AOC,'position');                    
pW=pannelsize(3);
pH=pannelsize(4);                
set(handles.pb_AOC_Plot,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*1/10 pH*2/10 pW*2.0/10 pH*6/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Plot figure',...
                        'backgroundcolor',handles.bgcolor{3});
set(handles.pb_AOC_SaveFig,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*4/10 pH*2/10 pW*2.0/10 pH*6/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Save figure',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'enable','off',...
                        'visible','off');
set(handles.pb_AOC_OK,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*4/10 pH*2/10 pW*2.0/10 pH*6/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','OK',...
                        'backgroundcolor',handles.bgcolor{3});
set(handles.pb_Cancel,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*7/10 pH*2/10 pW*2.0/10 pH*6/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Cancel',...
                        'backgroundcolor',handles.bgcolor{3});
% --------------------------------
%
handles.ObjEditEnable_AOC       = findobj('-regexp','Tag','_AOC_');
handles.ObjEditVisible_uRatio   = findobj('-regexp','Tag','_uRatio');
handles.ObjEditVisible_pRatio   = findobj('-regexp','Tag','_pRatio');


guidata(hObject, handles);
% lauch this function
Fcn_PreProcessing
%
%
handles = guidata(hObject);
% default enable settings
set(handles.ObjEditEnable_AOC,          'Enable','on');
set(handles.ObjEditVisible_uRatio,      'visible','on');
set(handles.ObjEditVisible_pRatio,      'visible','on');
%
guidata(hObject, handles);
% GUI_FREQ_EigCal_Initialization_check_ever_run(hObject)
handles = guidata(hObject);

GUI_FREQ_EigCal_pannel_appearance(hObject)
handles = guidata(hObject);
%
assignin('base','CI',CI)

%----------------------------------------
guidata(hObject, handles);
%

% -------------------------------------------------------------------------
%
% --- Outputs from this function are returned to the command line.
function varargout = GUI_FREQ_EigCal_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
try
varargout{1} = handles.output;
end

% --- Executes when user attempts to close figure.
function figure_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(hObject);

% ---------------------------Pannel initialization-------------------------
%
function edit_uRatio_min_Callback(hObject, eventdata, handles)
% hObject    handle to edit_uRatio_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uRatio_min      = str2double(get(handles.edit_uRatio_min,'string'));
uRatio_max      = str2double(get(handles.edit_uRatio_max,'string'));
if isnan(uRatio_min)
    set(handles.edit_uRatio_min, 'String', 0);
    errordlg('Input must be a number','Error');
end
if uRatio_min<=0
    set(handles.edit_uRatio_min, 'String', 0);
end
if uRatio_max <= uRatio_min
    set(handles.edit_uRatio_min, 'String', 0);
    set(handles.edit_uRatio_max, 'String', 1);
    errordlg('Minimum value input must be smaller than Maximum value input','Error');
end

% --- Executes during object creation, after setting all properties.
function edit_uRatio_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_uRatio_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_uRatio_max_Callback(hObject, eventdata, handles)
% hObject    handle to edit_uRatio_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uRatio_min      = str2double(get(handles.edit_uRatio_min,'string'));
uRatio_max      = str2double(get(handles.edit_uRatio_max,'string'));
if isnan(uRatio_max)
    set(handles.edit_uRatio_max, 'String', 1);
    errordlg('Input must be a number','Error');
end
if uRatio_max<=0
    set(handles.edit_uRatio_max, 'String', 1);
end
if uRatio_max <= uRatio_min
    set(handles.edit_uRatio_min, 'String', 0);
    set(handles.edit_uRatio_max, 'String', 1);
    errordlg('Minimum value input must be smaller than Maximum value input','Error');
end

% --- Executes during object creation, after setting all properties.
function edit_uRatio_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_uRatio_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_uRatio_SampNum_Callback(hObject, eventdata, handles)
% hObject    handle to edit_uRatio_SampNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uRatio_SampNum  = str2double(get(handles.edit_uRatio_SampNum,'string'));
if isnan(uRatio_SampNum )
    set(handles.edit_uRatio_SampNum, 'String', 11);
    errordlg('Input must be a number','Error');
end
if uRatio_SampNum <=0
    set(handles.edit_uRatio_SampNum, 'String', 11);
    errordlg('Input must be a positive integer','Error');
end
if rem(uRatio_SampNum ,1)~=0
    set(handles.edit_uRatio_SampNum, 'String', num2str(ceil(datEdit)));
end

% --- Executes during object creation, after setting all properties.
function edit_uRatio_SampNum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_uRatio_SampNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pb_SetScanDomain.
function pb_SetScanDomain_Callback(varargin)
hObject     = varargin{1};
handles     = guidata(hObject);
GUI_FREQ_EigCal_AD('GUI_FREQ_EigCal', handles.figure);

% --- Executes on button press in pb_CalEig.
function pb_CalEig_Callback(varargin)
hObject = varargin{1};
handles = guidata(hObject);
GUI_FREQ_EigCal_Eigenvalues_calculation(hObject);
GUI_FREQ_EigCal_PLOT(hObject)
%
% -------------------------------------------------------------------------
%

% -------------------------------------------------------------------------

%
% --- Executes on slider movement.
function slider_uRatio_Callback(hObject, eventdata, handles)
% hObject    handle to slider_uRatio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global CI
ValueSlider     = get(hObject,'Value');
indexShow       = round(ValueSlider);
set(handles.edit_uRatio,'string',num2str(CI.EIG.FDF.uRatioSp(indexShow)));
%
data_num(:,1)   = abs(imag(CI.EIG.Scan.EigValCol{indexShow})./2./pi);
data_num(:,2)   = real(CI.EIG.Scan.EigValCol{indexShow});
data_cell       = num2cell(data_num);
set(handles.uitable,'data',data_cell);         % Update the table
%
% Update the plotType choices
pop_PlotType = get(handles.pop_PlotType, 'Value');
switch pop_PlotType
    case 1
        set(handles.pop_numMode,        'visible','off',...
                                        'enable','on');
    case 2
        set(handles.pop_numMode,        'visible','on',...
                                        'enable','on',...
                                        'Value',1); 
        % Renew the number of modes with this u-ratio
        eigenvalue          = CI.EIG.Scan.EigValCol{indexShow};
        for k = 1:1:length(eigenvalue)
        StringMode{k} = ['Mode number: ' num2str(k)]; 
        end
        set(handles.pop_numMode,    'string',StringMode);
        
    case 3
        set(handles.pop_numMode,        'visible','on',...
                                        'enable','on',...
                                        'Value',1);
        % The number of modes is always the same as that of the first group
        eigenvalue          = CI.EIG.Scan.EigValCol{1};
        for k = 1:1:length(eigenvalue)
        StringMode{k} = ['Mode number: ' num2str(k)]; 
        end
        set(handles.pop_numMode,     'string',StringMode);
end

guidata(hObject, handles);
GUI_FREQ_EigCal_PLOT(hObject)

% --- Executes during object creation, after setting all properties.
function slider_uRatio_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_uRatio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function edit_uRatio_Callback(hObject, eventdata, handles)
% hObject    handle to edit_uRatio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_uRatio as text
%        str2double(get(hObject,'String')) returns contents of edit_uRatio as a double
datEdit = str2double(get(hObject, 'String'));
if isnan(datEdit)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end
set(handles.pb_PlotC_M,     'enable','off');
set(handles.pb_AOC_Plot,  'enable','off');
set(handles.edit_indexMode1,'enable','off');
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_uRatio_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_uRatio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%-----------Similar settings for pressure oscillation ratios---------------
%--------------------------------------------------------------------------

function edit_pRatio_min_Callback(hObject, eventdata, handles)
% hObject    handle to edit_pRatio_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pRatio_min      = str2double(get(handles.edit_pRatio_min,'string'));
pRatio_max      = str2double(get(handles.edit_pRatio_max,'string'));
if isnan(pRatio_min)
    set(handles.edit_pRatio_min, 'String', 0);
    errordlg('Input must be a number','Error');
end
if pRatio_min<=0
    set(handles.edit_pRatio_min, 'String', 0);
end
if pRatio_max <= pRatio_min
    set(handles.edit_pRatio_min, 'String', 0);
    set(handles.edit_pRatio_max, 'String', 1000);
    errordlg('Minimum value input must be smaller than Maximum value input','Error');
end
% Hints: get(hObject,'String') returns contents of edit_pRatio_min as text
%        str2double(get(hObject,'String')) returns contents of edit_pRatio_min as a double



% --- Executes during object creation, after setting all properties.
function edit_pRatio_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_pRatio_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_pRatio_max_Callback(hObject, eventdata, handles)
% hObject    handle to edit_pRatio_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pRatio_min      = str2double(get(handles.edit_pRatio_min,'string'));
pRatio_max      = str2double(get(handles.edit_pRatio_max,'string'));
if isnan(pRatio_max)
    set(handles.edit_pRatio_max, 'String', 1);
    errordlg('Input must be a number','Error');
end
if pRatio_max<=0
    set(handles.edit_pRatio_max, 'String', 1);
end
if pRatio_max <= pRatio_min
    set(handles.edit_pRatio_min, 'String', 0);
    set(handles.edit_pRatio_max, 'String', 1000);
    errordlg('Minimum value input must be smaller than Maximum value input','Error');
end

% Hints: get(hObject,'String') returns contents of edit_pRatio_max as text
%        str2double(get(hObject,'String')) returns contents of edit_pRatio_max as a double


% --- Executes during object creation, after setting all properties.
function edit_pRatio_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_pRatio_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_pRatio_SampNum_Callback(hObject, eventdata, handles)
% hObject    handle to edit_pRatio_SampNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pRatio_SampNum  = str2double(get(handles.edit_pRatio_SampNum,'string'));
if isnan(pRatio_SampNum )
    set(handles.edit_pRatio_SampNum, 'String', 11);
    errordlg('Input must be a number','Error');
end
if pRatio_SampNum <=0
    set(handles.edit_pRatio_SampNum, 'String', 11);
    errordlg('Input must be a positive integer','Error');
end
if rem(pRatio_SampNum ,1)~=0
    set(handles.edit_pRatio_SampNum, 'String', num2str(ceil(datEdit)));
end

% Hints: get(hObject,'String') returns contents of edit_pRatio_SampNum as text
%        str2double(get(hObject,'String')) returns contents of edit_pRatio_SampNum as a double


% --- Executes during object creation, after setting all properties.
function edit_pRatio_SampNum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_pRatio_SampNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider_pRatio_Callback(hObject, eventdata, handles)
% hObject    handle to slider_pRatio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global CI
ValueSlider     = get(hObject,'Value');
indexShow       = round(ValueSlider);
set(handles.edit_pRatio,'string',num2str(CI.EIG.FDF.pRatioSp(indexShow)));
%
data_num(:,1)   = abs(imag(CI.EIG.Scan.EigValCol{indexShow})./2./pi);
data_num(:,2)   = real(CI.EIG.Scan.EigValCol{indexShow});
data_cell       = num2cell(data_num);
set(handles.uitable,'data',data_cell);         % Update the table
%
% Update the plotType choices
pop_PlotType = get(handles.pop_PlotType, 'Value');
switch pop_PlotType
    case 1
        set(handles.pop_numMode,        'visible','off',...
                                        'enable','on');
    case 2
        set(handles.pop_numMode,        'visible','on',...
                                        'enable','on',...
                                        'Value',1); 
        % Renew the number of modes with this u-ratio
        eigenvalue          = CI.EIG.Scan.EigValCol{indexShow};
        for k = 1:1:length(eigenvalue)
        StringMode{k} = ['Mode number: ' num2str(k)]; 
        end
        set(handles.pop_numMode,    'string',StringMode);
        
    case 3
        set(handles.pop_numMode,        'visible','on',...
                                        'enable','on',...
                                        'Value',1);
        % The number of modes is always the same as that of the first group
        eigenvalue          = CI.EIG.Scan.EigValCol{1};
        for k = 1:1:length(eigenvalue)
        StringMode{k} = ['Mode number: ' num2str(k)]; 
        end
        set(handles.pop_numMode,     'string',StringMode);
end
guidata(hObject, handles);
GUI_FREQ_EigCal_PLOT(hObject)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider_pRatio_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_pRatio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit_pRatio_Callback(hObject, eventdata, handles)
% hObject    handle to edit_pRatio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit_pRatio as text
%        str2double(get(hObject,'String')) returns contents of edit_pRatio as a double
datEdit = str2double(get(hObject, 'String'));
if isnan(datEdit)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end
set(handles.pb_PlotC_M,     'enable','off');
set(handles.pb_AOC_Plot,  'enable','off');
set(handles.edit_indexMode1,'enable','off');
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit_pRatio_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_pRatio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% ---------------------------Pannel axes-----------------------------------
%
% --- Executes on selection change in pop_numMode.
function pop_numMode_Callback(hObject, eventdata, handles)
% hObject    handle to pop_numMode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pop_numMode contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_numMode

% --- Executes during object creation, after setting all properties.
function pop_numMode_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_numMode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in pop_PlotType.
function pop_PlotType_Callback(hObject, eventdata, handles)
% hObject    handle to pop_PlotType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pop_PlotType contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_PlotType
pop_PlotType = get(hObject,'Value');
set(handles.pb_AOC_SaveFig,     'enable','OFF');
% set(handles.pb_AOC_SaveFig,'enable','off');
switch pop_PlotType
    case 1
        set(handles.pop_numMode,        'visible','off',...
                                        'enable','on');
    otherwise
        set(handles.pop_numMode,        'visible','on',...
                                        'enable','on');
end
        

% --- Executes during object creation, after setting all properties.
function pop_PlotType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_PlotType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% ---------------------------Pannel AOC-----------------------------------
%
% --- Executes on button press in pb_AOC_Plot.
function pb_AOC_Plot_Callback(hObject, eventdata, handles)
% hObject    handle to pb_AOC_Plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% indexMode=str2num(get(handles.edit_indexMode1,'string'));
% Fcn_plot_eigenmode(indexMode)
GUI_FREQ_EigCal_PLOT(hObject)
set(handles.pb_AOC_SaveFig,     'enable','on');


% --- Executes on button press in pb_AOC_SaveFig.
function pb_AOC_SaveFig_Callback(varargin)
%
hObject = varargin{1};
handles = guidata(hObject);
pop_PLot            = get(handles.pop_PlotType,'Value');
Fig = figure;
set(Fig,        'units','points')
posFig = get(handles.figure,'position');
copyobj(handles.axes1, Fig);
copyobj(handles.axes2, Fig);

hAxes = get(Fig,'children');
set(hAxes(2),       'units','points',...
                    'position',[80 60 200 150],...
                    'ActivePositionProperty','position')
set(hAxes(1),       'units','points',...
                    'position',[80 210 200 150],...
                    'ActivePositionProperty','position')
pos1 = get(hAxes(1),'position');
pos2 = get(hAxes(2),'position');
switch pop_PLot
    case 1
%         try
            posAxesOuter = [0 0 450 400];
            set(hAxes(2),       'position',[80 80 300 300])
            set(hAxes(1),       'position',[80 80 300 300])
            pos1CT=get(hAxes(1),'position');
            pos2CT=get(hAxes(2),'position');
            colormap(hot);
            hcb=colorbar;
            ylimhcb  = handles.hColorbar.ylimit;       
            hcb_ylim = get(hcb,'ylim');
            hcb_yticklabel = round(linspace(ylimhcb(1),ylimhcb(2),6).*10)./10;
            hcb_ytick = linspace(hcb_ylim(1),hcb_ylim(2),6);
            set(hcb,'ytick',hcb_ytick,'yticklabel',hcb_yticklabel);
            set(hcb,'Fontsize',handles.FontSize(2),'box','on','Unit','points')
            set(hcb,'position',[380,80,10,300])
            set(hAxes(1),       'position',pos1CT)
            set(hAxes(2),       'position',pos2CT) 
%         catch
%         end
    otherwise
        posAxesOuter = [0 0 300 400];
end
set(Fig,        'units','points',...
                'position', [posFig(1)+0.5*posFig(3)-0.5*posAxesOuter(3),...
                            posFig(2)+0.5*posFig(4)-0.5*posAxesOuter(4),...
                            posAxesOuter(3:4)])  

% --- Executes on button press in pb_AOC_OK.
function pb_AOC_OK_Callback(hObject, eventdata, handles)
% hObject    handle to pb_AOC_OK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global CI
CI.IsRun.GUI_FREQ_EigCal = 1;
assignin('base','CI',CI);
delete(handles.figure);

% --- Executes on button press in pb_Cancel.
function pb_Cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pb_Cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(handles.figure);



function handles = Fcn_GUI_default_configuration(handles)
% This function is used to set the default GUI parameters
handles.bgcolor{1} = [1, 1, 1];
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


% --------------------------------------------------------------------
function uipushtool1_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pop_PLot            = get(handles.pop_PlotType,'Value');
Fig = figure;
set(Fig,        'units','points')
posFig = get(handles.figure,'position');
copyobj(handles.axes1, Fig);
copyobj(handles.axes2, Fig);

hAxes = get(Fig,'children');
set(hAxes(2),       'units','points',...
                    'position',[80 60 200 150],...
                    'ActivePositionProperty','position')
set(hAxes(1),       'units','points',...
                    'position',[80 210 200 150],...
                    'ActivePositionProperty','position')
pos1 = get(hAxes(1),'position');
pos2 = get(hAxes(2),'position');
switch pop_PLot
    case 1
%         try
            posAxesOuter = [0 0 450 400];
            set(hAxes(2),       'position',[80 80 300 300])
            set(hAxes(1),       'position',[80 80 300 300])
            pos1CT=get(hAxes(1),'position');
            pos2CT=get(hAxes(2),'position');
            colormap(hot);
            hcb=colorbar;
            ylimhcb  = handles.hColorbar.ylimit;       
            hcb_ylim = get(hcb,'ylim');
            hcb_yticklabel = round(linspace(ylimhcb(1),ylimhcb(2),6).*10)./10;
            hcb_ytick = linspace(hcb_ylim(1),hcb_ylim(2),6);
            set(hcb,'ytick',hcb_ytick,'yticklabel',hcb_yticklabel);
            set(hcb,'Fontsize',handles.FontSize(2),'box','on','Unit','points')
            set(hcb,'position',[380,80,10,300])
            set(hAxes(1),       'position',pos1CT)
            set(hAxes(2),       'position',pos2CT) 
%         catch
%         end
    otherwise
        posAxesOuter = [0 0 300 400];
end
set(Fig,        'units','points',...
                'position', [posFig(1)+0.5*posFig(3)-0.5*posAxesOuter(3),...
                            posFig(2)+0.5*posFig(4)-0.5*posAxesOuter(4),...
                            posAxesOuter(3:4)])  
