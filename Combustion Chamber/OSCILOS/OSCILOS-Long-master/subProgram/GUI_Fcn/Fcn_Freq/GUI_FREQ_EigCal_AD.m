function varargout = GUI_FREQ_EigCal_AD(varargin)
% GUI_FREQ_EIGCAL_AD MATLAB code for GUI_FREQ_EigCal_AD.fig
%      GUI_FREQ_EIGCAL_AD, by itself, creates a new GUI_FREQ_EIGCAL_AD or raises the existing
%      singleton*.
%
%      H = GUI_FREQ_EIGCAL_AD returns the handle to a new GUI_FREQ_EIGCAL_AD or the handle to
%      the existing singleton*.
%
%      GUI_FREQ_EIGCAL_AD('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_FREQ_EIGCAL_AD.M with the given input arguments.
%
%      GUI_FREQ_EIGCAL_AD('Property','Value',...) creates a new GUI_FREQ_EIGCAL_AD or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_FREQ_EigCal_AD_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_FREQ_EigCal_AD_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_FREQ_EigCal_AD

% Last Modified by GUIDE v2.5 10-Oct-2014 10:57:35

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_FREQ_EigCal_AD_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_FREQ_EigCal_AD_OutputFcn, ...
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
%
%--------------------------------------------------------------------------
%
% --- Executes just before GUI_FREQ_EigCal_AD is made visible.
function GUI_FREQ_EigCal_AD_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_FREQ_EigCal_AD (see VARARGIN)
indexEdit = 0;
switch indexEdit 
    case 0
        %--------------------------------------------------------------------------
        dontOpen = false;
        mainGuiInput = find(strcmp(varargin, 'GUI_FREQ_EigCal'));
        if (isempty(mainGuiInput)) ...
            || (length(varargin) <= mainGuiInput) ...
            || (~ishandle(varargin{mainGuiInput+1}))
            dontOpen = true;
        else % load from the main GUI
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
            GUI_FREQ_EigCal_AD_Initialization(hObject, eventdata, handles)
        end
        guidata(hObject, handles);
        handles.output = hObject;
        guidata(hObject, handles);
        if dontOpen
           disp('-----------------------------------------------------');
           disp('This is a subprogram. It cannot be run independently.') 
           disp('Please load the program "GUI_FREQ_EigCal'' from the ')
           disp('parent directory!')
           disp('-----------------------------------------------------');
        else
%            uiwait(hObject);
        end
    case 1
        global CI
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
        CI.Ru = 8.3145;
        CI.IsRun.GUI_FREQ_EigCal_AD = 0;
        assignin('base','CI',CI);                   % save the current information to the works
        guidata(hObject, handles); 
        GUI_FREQ_EigCal_AD_Initialization(hObject, eventdata, handles)
        guidata(hObject, handles);
        handles.output = hObject;
        guidata(hObject, handles);
end
%
%--------------------------------------------------------------------------
%
function Fcn_check_ever_run(varargin)
hObject = varargin{1};
handles = guidata(hObject);
global CI
switch CI.IsRun.GUI_FREQ_EigCal_AD 
    case 0
        CI.EIG.Scan.FreqMin = 0;
        CI.EIG.Scan.FreqMax = 1000;
        CI.EIG.Scan.FreqNum = 10;
        CI.EIG.Scan.GRMin   = -200;
        CI.EIG.Scan.GRMax   = 200;
        CI.EIG.Scan.GRNum   = 10;
    case 1
end
%
guidata(hObject, handles);
assignin('base','CI',CI);


%
%--------------------------------------------------------------------------
%
function GUI_FREQ_EigCal_AD_Initialization(varargin)
hObject = varargin{1};
handles = guidata(hObject);
guidata(hObject, handles);
Fcn_check_ever_run(hObject)
handles = guidata(hObject);
global CI
% ------------------------------------------
% positions reconfiguration
set(0, 'units', 'points');
screenSize  = get(0, 'ScreenSize');                     % get the screen size
sW          = handles.sW;                               % screen width
sH          = handles.sH;                               % screen height
FigW=sW*2/5;                                                % window width
FigH=sH.*1/2;                                           % window height
set(handles.figure,     'units', 'points',...
                        'position',[(screenSize(3)-FigW)./2 (screenSize(4)-FigH)./2 FigW FigH],...
                        'name','Scan range configuration',...
                        'color',handles.bgcolor{3});
%----------------------------------------
%
%
set(handles.uipanel_FREQ,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[FigW*0.5/20 FigH*11.25/20 FigW*19/20 FigH*8.5/20],...
                        'Title','Frequency scan domain setting',...
                        'visible','on',...
                        'highlightcolor',handles.bgcolor{3},...
                        'borderwidth',1,...
                        'fontsize',handles.FontSize(2),...
                        'backgroundcolor',handles.bgcolor{3}); 
pannelsize = get(handles.uipanel_FREQ,'position');
pW = pannelsize(3);
pH = pannelsize(4);   
%
set(handles.text_FREQ_a1,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*0.5/10 pH*6.0/10 pW*5/10 pH*1.25/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Frequency min. value [Hz]:',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'horizontalalignment','left',...
                        'visible','on');                   
set(handles.edit_FREQ_a1,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*6/10 pH*6/10 pW*3/10 pH*1.5/10],...
                        'fontsize',handles.FontSize(2),...
                        'string',CI.EIG.Scan.FreqMin,...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','right',...
                        'visible','on'); 
set(handles.text_FREQ_a2,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*0.5/10 pH*3.5/10 pW*5/10 pH*1.25/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Frequency max. value [Hz]:',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'horizontalalignment','left',...
                        'visible','on');                   
set(handles.edit_FREQ_a2,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*6/10 pH*3.5/10 pW*3/10 pH*1.5/10],...
                        'fontsize',handles.FontSize(2),...
                        'string',CI.EIG.Scan.FreqMax,...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','right',...
                        'visible','on');  
set(handles.text_FREQ_a3,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*0.5/10 pH*1.0/10 pW*5/10 pH*1.25/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Frequency samples [-]:',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'horizontalalignment','left',...
                        'visible','on');                   
set(handles.edit_FREQ_a3,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*6/10 pH*1/10 pW*3/10 pH*1.5/10],...
                        'fontsize',handles.FontSize(2),...
                        'string',CI.EIG.Scan.FreqNum,...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','right',...
                        'visible','on'); 
% ---------------
set(handles.uipanel_GR,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[FigW*0.5/20 FigH*2.5/20 FigW*19/20 FigH*8.5/20],...
                        'Title','Growth rate scan domain setting:',...
                        'visible','on',...
                        'highlightcolor',handles.bgcolor{3},...
                        'borderwidth',1,...
                        'fontsize',handles.FontSize(2),...
                        'backgroundcolor',handles.bgcolor{3}); 
pannelsize = get(handles.uipanel_GR,'position');
pW = pannelsize(3);
pH = pannelsize(4);   
%
set(handles.text_GR_a1,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*0.5/10 pH*6.0/10 pW*5/10 pH*1.25/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Growth rate min. value [1/s]:',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'horizontalalignment','left',...
                        'visible','on');                   
set(handles.edit_GR_a1,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*6/10 pH*6/10 pW*3/10 pH*1.5/10],...
                        'fontsize',handles.FontSize(2),...
                        'string',CI.EIG.Scan.GRMin,...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','right',...
                        'visible','on');  
set(handles.text_GR_a2,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*0.5/10 pH*3.5/10 pW*5/10 pH*1.25/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Growth rate max. value [1/s]:',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'horizontalalignment','left',...
                        'visible','on');                   
set(handles.edit_GR_a2,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*6/10 pH*3.5/10 pW*3/10 pH*1.5/10],...
                        'fontsize',handles.FontSize(2),...
                        'string',CI.EIG.Scan.GRMax,...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','right',...
                        'visible','on');  
set(handles.text_GR_a3,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*0.5/10 pH*1.0/10 pW*5/10 pH*1.25/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Growth samples [-]:',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'horizontalalignment','left',...
                        'visible','on');                   
set(handles.edit_GR_a3,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*6/10 pH*1/10 pW*3/10 pH*1.5/10],...
                        'fontsize',handles.FontSize(2),...
                        'string',CI.EIG.Scan.GRNum,...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','right',...
                        'visible','on');   
% ----------------------------------------
% pannel AOC                   
set(handles.uipanel_AOC,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[FigW*0.5/20 FigH*0/20 FigW*19/20 FigH*2.0/20],...
                        'Title','',...
                        'visible','on',...
                        'highlightcolor',handles.bgcolor{3},...
                        'borderwidth',1,...
                        'fontsize',handles.FontSize(2),...
                        'backgroundcolor',handles.bgcolor{3}); 
pannelsize=get(handles.uipanel_AOC,'position');                    
pW=pannelsize(3);
pH=pannelsize(4);                
set(handles.pb_Reset,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*1/10 pH*2/10 pW*2.0/10 pH*6/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Reset',...
                        'backgroundcolor',handles.bgcolor{3});
set(handles.pb_OK,...
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
                    
handles.ObjVisible_FREQ         = findobj('-regexp','Tag','_FREQ_');
handles.ObjVisible_GR           = findobj('-regexp','Tag','_GR_');
handles.ObjEditEnable_FREQ      = findobj('-regexp','Tag','edit_FREQ');
handles.ObjEditEnable_GR        = findobj('-regexp','Tag','edit_GR');
set(handles.ObjVisible_FREQ,                'visible','on')
set(handles.ObjVisible_GR,                  'visible','on')
set(handles.ObjEditEnable_FREQ,             'Enable','on')
set(handles.ObjEditEnable_GR,               'Enable','on')
guidata(hObject, handles);
%
%--------------------------------------------------------------------------
%
% --- Executes on button press in pb_OK.
function pb_OK_Callback(hObject, eventdata, handles)
% hObject    handle to pb_OK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global CI
CI.EIG.Scan.FreqMin = str2double(get(handles.edit_FREQ_a1,'String'));
CI.EIG.Scan.FreqMax = str2double(get(handles.edit_FREQ_a2,'String'));
CI.EIG.Scan.FreqNum = str2double(get(handles.edit_FREQ_a3,'String'));
CI.EIG.Scan.GRMin   = str2double(get(handles.edit_GR_a1,'String'));
CI.EIG.Scan.GRMax   = str2double(get(handles.edit_GR_a2,'String'));
CI.EIG.Scan.GRNum   = str2double(get(handles.edit_GR_a3,'String'));
CI.IsRun.GUI_FREQ_EigCal_AD = 1;
assignin('base','CI',CI); 
delete(handles.figure);
%
%--------------------------------------------------------------------------
%
% --- Executes on button press in pb_Reset.
function pb_Reset_Callback(hObject, eventdata, handles)
% hObject    handle to pb_Reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global CI
handles = guidata(hObject);
set(handles.edit_FREQ_a1,   'string',CI.EIG.Scan.FreqMin);
set(handles.edit_FREQ_a2,   'string',CI.EIG.Scan.FreqMax);
set(handles.edit_FREQ_a3,   'string',CI.EIG.Scan.FreqNum);
set(handles.edit_GR_a1,     'string',CI.EIG.Scan.GRMin);
set(handles.edit_GR_a2,     'string',CI.EIG.Scan.GRMax);
set(handles.edit_GR_a3,     'string',CI.EIG.Scan.GRNum);
guidata(hObject, handles);
%
%--------------------------------------------------------------------------
%
% --- Executes on button press in pb_Cancel.
function pb_Cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pb_Cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(handles.figure);
%
%--------------------------------------------------------------------------
%
% --- Outputs from this function are returned to the command line.
function varargout = GUI_FREQ_EigCal_AD_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
try
varargout{1} = handles.output;
end
%
%--------------------------------------------------------------------------
%
function edit_GR_a1_Callback(hObject, eventdata, handles)
% hObject    handle to edit_GR_a1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
a1 = str2double(get(handles.edit_GR_a1, 'string'));
a2 = str2double(get(handles.edit_GR_a2, 'string'));
a3 = str2double(get(handles.edit_GR_a3, 'string'));
if isnan(a1)
    set(handles.edit_GR_a1, 'String', -200);
    errordlg('Input must be a number','Error');
end
if a2 <= a1
    set(handles.edit_GR_a1, 'String', -200);
    set(handles.edit_GR_a2, 'String', 100);
    errordlg('Minimum value input must be smaller than Maximum value input','Error');
end

%
%--------------------------------------------------------------------------
%
% --- Executes during object creation, after setting all properties.
function edit_GR_a1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_GR_a1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%
%--------------------------------------------------------------------------
%
function edit_GR_a2_Callback(hObject, eventdata, handles)
% hObject    handle to edit_GR_a2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
a1 = str2double(get(handles.edit_GR_a1, 'string'));
a2 = str2double(get(handles.edit_GR_a2, 'string'));
a3 = str2double(get(handles.edit_GR_a3, 'string'));
if isnan(a2)
    set(handles.edit_GR_a2, 'String', 100);
    errordlg('Input must be a number','Error');
end
if a2 <= a1
    set(handles.edit_GR_a1, 'String', -200);
    set(handles.edit_GR_a2, 'String', 100);
    errordlg('Minimum value input must be smaller than Maximum value input','Error');
end
%--------------------------------------------------------------------------
%
% --- Executes during object creation, after setting all properties.
function edit_GR_a2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_GR_a2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%
%--------------------------------------------------------------------------
%
function edit_GR_a3_Callback(hObject, eventdata, handles)
% hObject    handle to edit_GR_a3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
a3 = str2double(get(handles.edit_GR_a3, 'string'));
if isnan(a3)
    set(handles.edit_GR_a3, 'String', 10);
    errordlg('Input must be a positive integer','Error');
end
if a3<= 0
    set(handles.edit_GR_a3, 'String', 10);
    errordlg('Input must be a positive integer','Error');
end
    
%
%--------------------------------------------------------------------------
%
% --- Executes during object creation, after setting all properties.
function edit_GR_a3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_GR_a3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%
%--------------------------------------------------------------------------
%
function edit_FREQ_a1_Callback(hObject, eventdata, handles)
% hObject    handle to edit_FREQ_a1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
a1 = str2double(get(handles.edit_FREQ_a1, 'string'));
a2 = str2double(get(handles.edit_FREQ_a2, 'string'));
a3 = str2double(get(handles.edit_FREQ_a3, 'string'));
if isnan(a1)
    set(handles.edit_FREQ_a1, 'String', 0);
    errordlg('Input must be a number','Error');
end
if a2 <= a1
    set(handles.edit_FREQ_a1, 'String', 0);
    set(handles.edit_FREQ_a2, 'String', 1000);
    errordlg('Minimum value input must be smaller than Maximum value input','Error');
end
%
%--------------------------------------------------------------------------
%   
% --- Executes during object creation, after setting all properties.
function edit_FREQ_a1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_FREQ_a1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%
%--------------------------------------------------------------------------
%
function edit_FREQ_a2_Callback(hObject, eventdata, handles)
% hObject    handle to edit_FREQ_a2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
a1 = str2double(get(handles.edit_FREQ_a1, 'string'));
a2 = str2double(get(handles.edit_FREQ_a2, 'string'));
a3 = str2double(get(handles.edit_FREQ_a3, 'string'));
if isnan(a2)
    set(handles.edit_FREQ_a2, 'String', 1000);
    errordlg('Input must be a number','Error');
end
if a2 <= a1
    set(handles.edit_FREQ_a1, 'String', 0);
    set(handles.edit_FREQ_a2, 'String', 1000);
    errordlg('Minimum value input must be smaller than Maximum value input','Error');
end
%
%--------------------------------------------------------------------------
%
% --- Executes during object creation, after setting all properties.
function edit_FREQ_a2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_FREQ_a2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%
%--------------------------------------------------------------------------
%
function edit_FREQ_a3_Callback(hObject, eventdata, handles)
% hObject    handle to edit_FREQ_a3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
a3 = str2double(get(handles.edit_FREQ_a3, 'string'));
if isnan(a3)
    set(handles.edit_FREQ_a3, 'String', 10);
    errordlg('Input must be a positive integer','Error');
end
if a3<= 0
    set(handles.edit_FREQ_a3, 'String', 10);
    errordlg('Input must be a positive integer','Error');
end

%
%--------------------------------------------------------------------------
%
% --- Executes during object creation, after setting all properties.
function edit_FREQ_a3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_FREQ_a3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%
%--------------------------------------------------------------------------
%
% --- Executes when user attempts to close figure.
function figure_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(hObject);
%
%------------------------------end-----------------------------------------