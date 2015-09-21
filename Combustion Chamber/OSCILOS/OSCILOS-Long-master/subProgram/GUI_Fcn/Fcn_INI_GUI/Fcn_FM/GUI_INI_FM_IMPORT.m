function varargout = GUI_INI_FM_IMPORT(varargin)
% GUI_INI_FM_IMPORT MATLAB code for GUI_INI_FM_IMPORT.fig
%      GUI_INI_FM_IMPORT, by itself, creates a new GUI_INI_FM_IMPORT or raises the existing
%      singleton*.
%
%      H = GUI_INI_FM_IMPORT returns the handle to a new GUI_INI_FM_IMPORT or the handle to
%      the existing singleton*.
%
%      GUI_INI_FM_IMPORT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_INI_FM_IMPORT.M with the given input arguments.
%
%      GUI_INI_FM_IMPORT('Property','Value',...) creates a new GUI_INI_FM_IMPORT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_INI_FM_IMPORT_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_INI_FM_IMPORT_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_INI_FM_IMPORT

% Last Modified by GUIDE v2.5 09-Feb-2015 10:35:31

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_INI_FM_IMPORT_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_INI_FM_IMPORT_OutputFcn, ...
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
%-------------------------------------------------------------------------
%
% --- Executes just before GUI_INI_FM_IMPORT is made visible.
function GUI_INI_FM_IMPORT_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_INI_FM_IMPORT (see VARARGIN)
%--------------------------------------------------------------------------
indexEdit = 0;
switch indexEdit 
    case 0
        %--------------------------------------------------------------------------
        dontOpen = false;
        mainGuiInput = find(strcmp(varargin, 'GUI_INI_FMEXP'));
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
            GUI_INI_FM_IMPORT_Initialization(hObject, eventdata, handles)
        end
        handles = guidata(hObject);  
        guidata(hObject, handles);
        handles.output = hObject;
        guidata(hObject, handles);
        if dontOpen
           disp('-----------------------------------------------------');
           disp('This is a subprogram. It cannot be run independently.') 
           disp('Please load the program "GUI_INI_FMEXP'' from the ')
           disp('parent directory!')
           disp('-----------------------------------------------------');
        else
%            uiwait(hObject);
        end
    case 1
        global HP
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
        handles.indexApp = 0;
        HP.FMEXP.indexIMPORT = 1;
        assignin('base','HP',HP);                   % save the current information to the works
        guidata(hObject, handles);  
        GUI_INI_FM_IMPORT_Initialization(hObject, eventdata, handles)
        handles = guidata(hObject);  
        guidata(hObject, handles);
        handles.output = hObject;
        guidata(hObject, handles);
end
%
%-------------------------------------------------------------------------
%
function GUI_INI_FM_IMPORT_Initialization(varargin)
hObject = varargin{1};
handles = guidata(hObject);    
% positions reconfiguration
set(0, 'units', 'points');
screenSize  = get(0, 'ScreenSize');                     % get the screen size
sW          = handles.sW;                               % screen width
sH          = handles.sH ;                              % screen height
FigW=sW.*4/5;                                           % window width
FigH=sH*2/3;                                             % window height
set(handles.figure,     'units', 'points',...
                        'position',[(screenSize(3)-FigW)./2 (screenSize(4)-FigH)./2 FigW FigH],...
                        'name','Import flame transfer function for a velocity ratio (u''/u_mean) and corresponding fitting',...
                        'color',handles.bgcolor{3});
%--------------------------------------------------------------------------
% pannel axes
set(handles.uipanel_Axes,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[FigW*10.5/20 FigH*2.0/20 FigW*9/20 FigH*17.75/20],...
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
                        'Fontunits','points',...
                        'position',[pW*2/10 pH*4.5/10 pW*6/10 pH*3.0/10],...
                        'fontsize',handles.FontSize(1),...
                        'color',handles.bgcolor{1},...
                        'box','on');  
set(handles.axes2,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*2/10 pH*1.5/10 pW*6/10 pH*3.0/10],...
                        'fontsize',handles.FontSize(1),...
                        'color',handles.bgcolor{1},...
                        'box','on');       
set(handles.pop_plot,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*0.5/10 pH*8.5/10 pW*5.0/10 pH*1/10],...
                        'fontsize',handles.FontSize(2),...
                        'string',{  'Transfer function';...
                                    'Poles and zeros'},...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','left',...
                        'enable','on',...
                        'value',1);  
set(handles.text_fRange,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*0.5/10 pH*7.75/10 pW*5/10 pH*0.6/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Plot frequency range: [Hz]',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'horizontalalignment','left');                         
set(handles.edit_fRange,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*6/10 pH*7.875/10 pW*3/10 pH*0.6/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','[1 400]',...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','right',...
                        'Enable','on');
guidata(hObject, handles);

%-----------------------------------------
% pannels import data
set(handles.uipanel_Import,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[FigW*0.5/20 FigH*14.75/20 FigW*9.5/20 FigH*5/20],...
                        'Title','Load gain and phase data from a mat file',...
                        'visible','on',...
                        'highlightcolor',handles.bgcolor{3},...
                        'borderwidth',1,...
                        'fontsize',handles.FontSize(2),...
                        'backgroundcolor',handles.bgcolor{3}); 
pannelsize=get(handles.uipanel_Import,'position');
pW=pannelsize(3);
pH=pannelsize(4);
set(handles.text_uRatio,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*0.5/10 pH*5.25/10 pW*6/10 pH*2.0/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Velocity ratio: [-]',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'horizontalalignment','left');                         
set(handles.edit_uRatio,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*6.5/10 pH*5.25/10 pW*3/10 pH*2.5/10],...
                        'fontsize',handles.FontSize(2),...
                        'string',0,...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','right',...
                        'Enable','on');
currentFolder = pwd;
set(handles.edit_dir,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*0.5/10 pH*1.5/10 pW*5.5/10 pH*2.5/10],...
                        'fontsize',handles.FontSize(2),...
                        'string',currentFolder,...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','left');                   
set(handles.pb_Import,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*6.5/10 pH*1.5/10 pW*3/10 pH*2.5/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Import data...',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'horizontalalignment','right',...
                        'Enable','on'); 
%-----------------------------------------
% pannel sys source pannel
set(handles.uipanel_source,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[FigW*0.5/20 FigH*11.5/20 FigW*9.5/20 FigH*3.0/20],...
                        'Title','Fitted transfer function source',...
                        'visible','on',...
                        'highlightcolor',handles.bgcolor{3},...
                        'borderwidth',1,...
                        'fontsize',handles.FontSize(2),...
                        'backgroundcolor',handles.bgcolor{3}); 
pannelsize=get(handles.uipanel_source,'position');
pW=pannelsize(3);
pH=pannelsize(4);                    
set(handles.text_source,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*0.5/10 pH*2/10 pW*5/10 pH*3/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Source:',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'horizontalalignment','left');                         
set(handles.pop_source,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*4.5/10 pH*2/10 pW*5.0/10 pH*4/10],...
                        'fontsize',handles.FontSize(2),...
                        'string',{  'fitfrd fitting';...
                                    'Load from file'},...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','left',...
                        'enable','on',...
                        'value',1);  
%-----------------------------------------
% pannel fitting setting 
set(handles.uipanel_FIT,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[FigW*0.5/20 FigH*2.0/20 FigW*9.5/20 FigH*9.25/20],...
                        'Title','Command fitfrd configuration...',...
                        'visible','on',...
                        'highlightcolor',handles.bgcolor{3},...
                        'borderwidth',1,...
                        'fontsize',handles.FontSize(2),...
                        'backgroundcolor',handles.bgcolor{3}); 
pannelsize=get(handles.uipanel_FIT,'position');
pW=pannelsize(3);
pH=pannelsize(4);
set(handles.text_FIT_a1,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*0.5/10 pH*7.5/10 pW*6/10 pH*1/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Min. value of frequency range: [-]',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'horizontalalignment','left');                         
set(handles.edit_FIT_a1,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*6.5/10 pH*7.5/10 pW*3/10 pH*1.25/10],...
                        'fontsize',handles.FontSize(2),...
                        'string',0,...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','right',...
                        'Enable','on');
set(handles.text_FIT_a2,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*0.5/10 pH*5.75/10 pW*6/10 pH*1/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Max. value of frequency range: [-]',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'horizontalalignment','left');                         
set(handles.edit_FIT_a2,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*6.5/10 pH*5.75/10 pW*3/10 pH*1.25/10],...
                        'fontsize',handles.FontSize(2),...
                        'string',1000,...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','right',...
                        'Enable','on');
set(handles.text_FIT_a3,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*0.5/10 pH*4/10 pW*6/10 pH*1/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Time delay correction: [ms]',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'horizontalalignment','left');                         
set(handles.edit_FIT_a3,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*6.5/10 pH*4/10 pW*3/10 pH*1.25/10],...
                        'fontsize',handles.FontSize(2),...
                        'string',num2str(-4),...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','right',...
                        'Enable','on');
set(handles.text_FIT_a4,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*0.5/10 pH*2.25/10 pW*6/10 pH*1/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','State dimension (fitting order): [-]',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'horizontalalignment','left');                         
set(handles.edit_FIT_a4,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*6.5/10 pH*2.25/10 pW*3/10 pH*1.25/10],...
                        'fontsize',handles.FontSize(2),...
                        'string',num2str(6),...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','right',...
                        'Enable','on');
set(handles.text_FIT_a5,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*0.5/10 pH*0.0/10 pW*6/10 pH*2/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Relative degree of denominator compared to numerator: [-]',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'horizontalalignment','left');                         
set(handles.edit_FIT_a5,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*6.5/10 pH*0.5/10 pW*3/10 pH*1.25/10],...
                        'fontsize',handles.FontSize(2),...
                        'string',num2str(0),...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','right',...
                        'Enable','on');
set(handles.pb_sys_load,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*0.5/10 pH*1/10 pW*9/10 pH*1.5/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Load the transfer function models and time delay ',...
                        'backgroundcolor',handles.bgcolor{3});
%----------------------------------------
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
set(handles.pb_fitfrd,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*1/10 pH*1/10 pW*2.0/10 pH*7/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Fit and plot',...
                        'backgroundcolor',handles.bgcolor{3});
set(handles.pb_SaveFitting,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*4/10 pH*1/10 pW*2.0/10 pH*7/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Save fitting',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'enable','off');
set(handles.pb_SaveFig,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*4/10 pH*1/10 pW*2.0/10 pH*7/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Save figure',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'enable','off',...
                        'visible','off');
set(handles.pb_Cancel,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*7/10 pH*1/10 pW*2.0/10 pH*7/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Cancel',...
                        'backgroundcolor',handles.bgcolor{3});
% ----------------------------                                 
guidata(hObject, handles);
GUI_INI_FM_IMPORT_edit_update(hObject)
handles = guidata(hObject); 
guidata(hObject, handles);
%
%-------------------------------------------------------------------------
%
function GUI_INI_FM_IMPORT_edit_update(varargin)
hObject = varargin{1};
handles = guidata(hObject); 
handles.ObjVisible_FIT          = findobj('-regexp','Tag','FIT');
handles.ObjEditEnable_FIT       = findobj('-regexp','Tag','edit_FIT');
handles.ObjVisible_source       = findobj('-regexp','Tag','_source');
set(handles.ObjVisible_FIT,         'visible','on')
set(handles.ObjEditEnable_FIT,      'enable', 'on')
global HP
guidata(hObject, handles);
switch HP.FMEXP.indexIMPORT 
    case 1
        set(handles.pb_Import,              'enable','on'); 
    case 2  % edit selected FTF
        handles.FMFit   = HP.FMEXP.FTF{HP.FMEXP.indexModify};        
        fmin            = handles.FMFit.freq_band(1);
        fmax            = handles.FMFit.freq_band(2);
        set(handles.edit_uRatio,    'string', num2str(handles.FMFit.uRatio));
        set(handles.edit_FIT_a1,    'string', num2str(fmin));
        set(handles.edit_FIT_a2,    'string', num2str(fmax));
        set(handles.edit_FIT_a3,    'string', num2str(handles.FMFit.tau_correction.*1000));                         
        set(handles.edit_FIT_a4,    'string', num2str(handles.FMFit.fitting_level));
        set(handles.edit_FIT_a5,    'string', num2str(handles.FMFit.RD));
        set(handles.pb_Import,      'Enable','off');
        try
            set(handles.pop_source,     'value',handles.FMFit.popSourceValue);
        catch
            set(handles.pop_source,     'value',1);
        end
        set(handles.pb_fitfrd,      'Enable','on');
        set(handles.pb_SaveFitting, 'Enable','on');
        % -------------
        guidata(hObject, handles);
        Fcn_GUI_INI_FM_IMPORT_Plot(hObject)
        handles = guidata(hObject);  
end   
Fcn_sys_source_pop_update(hObject);
%
%-------------------------------------------------------------------------
%
function Fcn_sys_source_pop_update(varargin)
hObject = varargin{1};
handles = guidata(hObject);
popValue = get(handles.pop_source,'value');
switch popValue
    case 1                
        set(handles.edit_FIT_a3,        'enable','on');
        set(handles.text_FIT_a4,        'visible','on'); 
        set(handles.text_FIT_a5,        'visible','on'); 
        set(handles.edit_FIT_a4,        'visible','on'); 
        set(handles.edit_FIT_a5,        'visible','on');
        set(handles.pb_sys_load,        'visible','off'); 
        set(handles.pb_fitfrd,         'enable','on');
    case 2
        set(handles.edit_FIT_a3,        'enable','off');
        set(handles.text_FIT_a4,        'visible','off'); 
        set(handles.text_FIT_a5,        'visible','off'); 
        set(handles.edit_FIT_a4,        'visible','off'); 
        set(handles.edit_FIT_a5,        'visible','off');
        set(handles.pb_sys_load,        'visible','on');
        set(handles.pb_fitfrd,         'enable','off');
end


%
%-------------------------------------------------------------------------
%

% --- Outputs from this function are returned to the command line.
function varargout = GUI_INI_FM_IMPORT_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
varargout{1} = handles.output;
end
%
%-------------------------------------------------------------------------
%
% --- Executes on button press in pb_Import.
function pb_Import_Callback(hObject, eventdata, handles)
% hObject    handle to pb_Import (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles       = guidata(hObject);
currentFolder = pwd;
[filename, pathname] = uigetfile({  '*.mat'}, ...
                                    'Pick a file',...
                                     currentFolder);
if filename~=0
    addpath(pathname)               % add directory to search path
    fDataTemp = load(filename);
handles.FMFit.Gain_exp  = fDataTemp.Gain;
handles.FMFit.Phase_exp = fDataTemp.Phase;
guidata(hObject, handles);
%
Fcn_GUI_INI_FM_IMPORT_Plot(hObject)
%
handles = guidata(hObject);
fmin = 0;
fmax = ceil(max(max(handles.FMFit.Gain_exp(:,1)),...
                max(handles.FMFit.Phase_exp(:,1)))./100).*100;
set(handles.edit_FIT_a1,      'string', num2str(fmin));
set(handles.edit_FIT_a2,      'string', num2str(fmax));
set(handles.pb_SaveFitting, 'enable','off')
set(handles.pb_SaveFig,     'enable','off')
end
% Update handles structure
guidata(hObject, handles);
%
%-------------------------------------------------------------------------
%
% --- Executes on button press in pb_sys_load.
function pb_sys_load_Callback(hObject, eventdata, handles)
% hObject    handle to pb_sys_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles       = guidata(hObject);
currentFolder = pwd;
[filename, pathname] = uigetfile({  '*.mat'}, ...
                                    'Pick a file',...
                                     currentFolder);
if filename~=0
    addpath(pathname)               % add directory to search path
    load(filename);
    handles.sysFit = sys;
    set(handles.edit_FIT_a3, 'string', num2str(td.*1e3));
    set(handles.pb_fitfrd,         'enable','on');
    guidata(hObject, handles);
end
%
%-------------------------------------------------------------------------
%
function Fcn_GUI_INI_FM_IMPORT_Plot(varargin)
hObject     = varargin{1};
handles     = guidata(hObject);
hAxes1      = handles.axes1;
hAxes2      = handles.axes2;
fontSize1   = handles.FontSize(1);
fontSize2   = handles.FontSize(2);
colorPlot   = {'r','k'};
cla(hAxes1,'reset')
axes(hAxes1)
hold on
plot(hAxes1,  handles.FMFit.Gain_exp(:,1),...
              handles.FMFit.Gain_exp(:,2),'s','color',colorPlot{1},...
              'Linewidth',2,'markersize',6);
try
    plot(hAxes1,handles.FMFit.xfit,...
                abs(handles.FMFit.yfit),'-','color',colorPlot{2},'Linewidth',2); 
end
ymax        = ceil(max(max(handles.FMFit.Gain_exp(:,2))).*10)./10;
ymin        = floor(min(min(handles.FMFit.Gain_exp(:,2))).*10)./10;
ylimitUD    = [ymin-0.25*(ymax-ymin) ymax+0.25*(ymax-ymin)];
ytickUD     = linspace(ylimitUD(1),ylimitUD(2),7);
xmin        = 0;
xmax        = ceil( max(max(handles.FMFit.Gain_exp(:,1)),...
                    max(handles.FMFit.Gain_exp(:,1)))./100).*100;
for ss=1:length(ytickUD)
    yticklabelUD{ss}=num2str(ytickUD(ss));
end
yticklabelUD{1}='';
set(hAxes1,'YColor','k','Box','on','ygrid','on','xgrid','on');
set(hAxes1,'FontName','Helvetica','FontSize',fontSize1,'LineWidth',1)
xlabel(hAxes1,'','Color','k','Interpreter','LaTex','FontSize',fontSize1);
ylabel(hAxes1,'Gain [-]','Color','k','Interpreter','LaTex','FontSize',fontSize1)
set(hAxes1,'xlim',[xmin xmax],'xtick',xmin:100:xmax,'xticklabel',{},...
'YAxisLocation','left','Color','w');
set(hAxes1,'ylim',ylimitUD,'yTick',ytickUD,'yticklabel',yticklabelUD)
hold off
%--------------------------------
cla(hAxes2,'reset')
axes(hAxes2)
hold on
plot(hAxes2,    handles.FMFit.Phase_exp(:,1),...
                handles.FMFit.Phase_exp(:,2)./pi,'s','color',colorPlot{1},...
                'Linewidth',2,'markersize',6);
try
    plot(hAxes2,handles.FMFit.xfit,...
                unwrap(angle(handles.FMFit.yfit),1.9*pi)./pi,'-',...
                'color',colorPlot{2},'Linewidth',2); 
end
ymax        = 0;
ymin        = round(min(handles.FMFit.Phase_exp(:,2)./pi));
yTickset    = ymin:1:ymax;
set(hAxes2,'YColor','k','Box','on','ygrid','on','xgrid','on');
set(hAxes2,'FontName','Helvetica','FontSize',fontSize1,'LineWidth',1)
xlabel(hAxes2,'$f$ [Hz]','Color','k','Interpreter','LaTex','FontSize',fontSize1);
ylabel(hAxes2,'Phase/$\pi$ [-]','Color','k','Interpreter','LaTex','FontSize',fontSize1)
set(hAxes2,'xlim',get(hAxes1,'xlim'),'xTick',get(hAxes1,'xTick'),...
'YAxisLocation','left','Color','w');
set(hAxes2,'ylim',[ymin ymax],'yTick',yTickset)
hold off     
%
guidata(hObject,handles);
%
%--------------------------------------------------------------------------
%
% --- Executes on button press in pb_Cancel.
function pb_Cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pb_Cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
delete(handles.figure);
%
%--------------------------------------------------------------------------
%
% --- Executes on button press in pb_SaveFig.
function pb_SaveFig_Callback(varargin)
hObject = varargin{1};
handles = guidata(hObject);
Fig = figure;
set(Fig,        'units','points')
posFig = get(handles.figure,'position');
copyobj(handles.axes1, Fig);
copyobj(handles.axes2, Fig);

hAxes = get(Fig,'children');
set(hAxes(1),       'units','points',...
                    'position',[80 60 200 150],...
                    'ActivePositionProperty','position')
set(hAxes(2),       'units','points',...
                    'position',[80 210 200 150],...
                    'ActivePositionProperty','position')
pos1=get(hAxes(1),'position');
pos2=get(hAxes(2),'position');
posAxesOuter = [0 0 310 400];
set(Fig,        'units','points',...
                'position', [posFig(1)+0.5*posFig(3)-0.5*posAxesOuter(3),...
                            posFig(2)+0.5*posFig(4)-0.5*posAxesOuter(4),...
                            posAxesOuter(3:4)])       
%
%--------------------------------------------------------------------------
%
% --- Executes on button press in pb_fitfrd.
function pb_fitfrd_Callback(varargin)
hObject                             = varargin{1};
handles                             = guidata(hObject);
fmin                                = str2double(get(handles.edit_FIT_a1,'string'));
fmax                                = str2double(get(handles.edit_FIT_a2,'string'));
handles.FMFit.uRatio                = str2double(get(handles.edit_uRatio,'string'));
handles.FMFit.freq_band             = [fmin, fmax];
handles.FMFit.tau_correction        = str2double(get(handles.edit_FIT_a3,'string'))./1000;
handles.FMFit.fitting_level         = str2double(get(handles.edit_FIT_a4,'string'));
handles.FMFit.RD                    = str2double(get(handles.edit_FIT_a5,'string'));
handles.FMFit.popSourceValue        = get(handles.pop_source,'value');
switch handles.FMFit.popSourceValue
    case 1
[   handles.FMFit.F,...
    handles.FMFit.den,...
    handles.FMFit.num ] = Fcn_fitfrd(   handles.FMFit.Gain_exp,...
                                            handles.FMFit.Phase_exp,...
                                            handles.FMFit.freq_band,...
                                            handles.FMFit.tau_correction,...
                                            handles.FMFit.fitting_level,...
                                            handles.FMFit.RD,...
                                            [],...
                                            handles.FMFit.popSourceValue);
    case 2
[   handles.FMFit.F,...
    handles.FMFit.den,...
    handles.FMFit.num ] = Fcn_fitfrd(   handles.FMFit.Gain_exp,...
                                            handles.FMFit.Phase_exp,...
                                            handles.FMFit.freq_band,...
                                            handles.FMFit.tau_correction,...
                                            handles.FMFit.fitting_level,...
                                            handles.FMFit.RD,...
                                            handles.sysFit,...
                                            handles.FMFit.popSourceValue);
end
        
% used for plots
handles.FMFit.xfit      =   fmin:1:fmax;
try
    fRange = str2num(get(handles.edit_fRange,'String'));
    handles.FMFit.xfit = fRange(1):1:fRange(2);
end

s                       =   2*pi*handles.FMFit.xfit*1i;
handles.FMFit.yfit      =   polyval(handles.FMFit.num,s)./...
                            polyval(handles.FMFit.den,s).*exp(-handles.FMFit.tau_correction.*s);
% Update handles structure
guidata(hObject, handles);
%
Fcn_GUI_INI_FM_IMPORT_Plot(hObject)
%
handles = guidata(hObject);
set(handles.pb_SaveFitting, 'enable','on')
set(handles.pb_SaveFig,     'enable','on')
guidata(hObject,handles);
%
%--------------------------------------------------------------------------
%
% --- Executes on button press in pb_SaveFitting.
function pb_SaveFitting_Callback(varargin)
hObject = varargin{1};
handles = guidata(hObject);
global HP
set(handles.pb_SaveFig, 'enable','on')
switch HP.FMEXP.indexIMPORT 
    case 1  % add new FTF
        nCount                      = HP.FMEXP.nFTF + 1;
        HP.FMEXP.nFTF              = nCount;
        HP.FMEXP.FTF{nCount}       = handles.FMFit;
        HP.FMEXP.uRatio(nCount)    = handles.FMFit.uRatio; 
        assignin('base','HP',HP);                  
        main = handles.MainGUI;
        if(ishandle(main))
            mainHandles             = guidata(main);
            changeMain              = mainHandles.listbox_EXP;
            String_Listbox          = get(changeMain,'string');
            nLength                 = size(String_Listbox);
            String_Listbox{nLength(1)+1}...
                = ['Velocity ratio u''/u_mean = ' num2str(handles.FMFit.uRatio)];
            set(changeMain,'string',String_Listbox,'value',nCount);
        end   
    case 2  % edit selected FTF
        nCount                      = HP.FMEXP.indexModify;
        HP.FMEXP.FTF{nCount}       = handles.FMFit;
        HP.FMEXP.uRatio(nCount)    = handles.FMFit.uRatio; 
        assignin('base','HP',HP);                  
        main = handles.MainGUI;
        if(ishandle(main))
            mainHandles             = guidata(main);
            changeMain              = mainHandles.listbox_EXP;
            String_Listbox          = get(changeMain,'string');
            String_Listbox{nCount}...
                = ['Velocity ratio u''/u_mean = ' num2str(handles.FMFit.uRatio)];
            set(changeMain,'string',String_Listbox,'value',nCount);
        end 
end  
assignin('base','HP',HP);                   % save the current information to the works
delete(handles.figure);
%
%--------------------------------------------------------------------------
%
% --- Executes when user attempts to close figure.
function figure_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(hObject);
clear TEMP
%
%--------------------------------------------------------------------------
%
function [F,den,num]=Fcn_fitfrd(FGain, FPhase, freq_band, tau_correction,fitting_level,RD,sys,flag)
% This function is used to get the transfer function of the unflagged open
% pipe situation
% freq_band should be [fmin fmax]
% tau_correction is used to better fitting
fmin        = freq_band(1);
fmax        = freq_band(2);
F.Freq      = fmin:1:fmax;
F.Gain      = abs(interp1(FGain(:,1)',FGain(:,2)',F.Freq,'linear','extrap'));
F.Phase     = interp1(FPhase(:,1)',FPhase(:,2)',F.Freq,'linear','extrap');
F.FTF       = F.Gain.*exp(1i.*F.Phase);
F.FTF_corr  = F.FTF.*exp(1i.*2*pi*F.Freq.*tau_correction);
%-------------------
switch flag
    case 1
        F_frd       = frd(F.FTF_corr,2*pi*F.Freq);
        F_fitfrd    = fitfrd(F_frd,fitting_level,RD);
        F_fitfrd_TF = tf(F_fitfrd);
        [num,den]   = tfdata(F_fitfrd_TF,'v');
    case 2   % in case the sys is from an external file
        [num,den]   = tfdata(sys,'v');
end


%
%--------------------------------------------------------------------------
%
function edit_uRatio_Callback(hObject, eventdata, handles)
% hObject    handle to edit_uRatio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_uRatio as text
%        str2double(get(hObject,'String')) returns contents of edit_uRatio as a double
uRatio = str2double(get(handles.edit_uRatio,'string'));
if isnan(uRatio)
    set(handles.edit_uRatio, 'String', 0);
    errordlg('Input must be a number','Error');
end
if uRatio<0
    set(handles.edit_uRatio, 'String', 0);
    errordlg('Input must be a positive number','Error');
end
%
%--------------------------------------------------------------------------
%
% --- Executes during object creation, after setting all properties.
function edit_uRatio_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_uRatio (see GCBO)
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
function edit_dir_Callback(hObject, eventdata, handles)
% hObject    handle to edit_dir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_dir as text
%        str2double(get(hObject,'String')) returns contents of edit_dir as a double
%
%--------------------------------------------------------------------------
%
% --- Executes during object creation, after setting all properties.
function edit_dir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_dir (see GCBO)
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
function edit_FIT_a1_Callback(hObject, eventdata, handles)
% hObject    handle to edit_FIT_a1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

fmin      = str2double(get(handles.edit_FIT_a1,'string'));
fmax      = str2double(get(handles.edit_FIT_a2,'string'));
if isnan(fmin)
    set(hObject, 'String', 0);
    errordlg('Input must be a nonnegative number','Error');
end
if fmin<=0
    set(hObject, 'String', 0);
    errordlg('Input must be a nonnegative number','Error');
end
if fmax <= fmin
    set(hObject, 'String', 0);
    errordlg('Minimum value input must be smaller than the Maximum value input','Error');
end
set(handles.pb_SaveFitting, 'enable','off')
set(handles.pb_SaveFig,     'enable','off')
%
%--------------------------------------------------------------------------
%
% --- Executes during object creation, after setting all properties.
function edit_FIT_a1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_FIT_a1 (see GCBO)
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
function edit_FIT_a2_Callback(hObject, eventdata, handles)
% hObject    handle to edit_FIT_a2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_FIT_a2 as text
%        str2double(get(hObject,'String')) returns contents of edit_FIT_a2 as a double
fmin      = str2double(get(handles.edit_FIT_a1,'string'));
fmax      = str2double(get(handles.edit_FIT_a2,'string'));
fDelta = 500;
if isnan(fmax)
    set(hObject, 'String', num2str(fmin+fDelta));
    errordlg('Input must be a nonnegative number','Error');
end
if fmax <= fmin
    set(hObject, 'String', num2str(fmin+fDelta));
    errordlg('Minimum value input must be smaller than the Maximum value input','Error');
end
set(handles.pb_SaveFitting, 'enable','off')
set(handles.pb_SaveFig,     'enable','off')
%
%--------------------------------------------------------------------------
%
% --- Executes during object creation, after setting all properties.
function edit_FIT_a2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_FIT_a2 (see GCBO)
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
function edit_FIT_a3_Callback(hObject, eventdata, handles)
% hObject    handle to edit_FIT_a3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_FIT_a3 as text
%        str2double(get(hObject,'String')) returns contents of edit_FIT_a3 as a double
datEdit = str2double(get(hObject, 'String'));
if isnan(datEdit)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end
set(handles.pb_SaveFitting, 'enable','off')
set(handles.pb_SaveFig,     'enable','off')
%
%--------------------------------------------------------------------------
%
% --- Executes during object creation, after setting all properties.
function edit_FIT_a3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_FIT_a3 (see GCBO)
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
function edit_FIT_a4_Callback(hObject, eventdata, handles)
% hObject    handle to edit_FIT_a4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_FIT_a4 as text
%        str2double(get(hObject,'String')) returns contents of edit_FIT_a4 as a double
datEdit         = str2double(get(hObject, 'String'));
if isnan(datEdit)
    set(hObject, 'String', 2);
    errordlg('Input must be a number','Error');
end
if datEdit <=0
    set(hObject, 'String', 2);
    errordlg('Input must be a positive integer','Error');
end
if rem(datEdit ,1)~=0
    set(hObject, 'String', num2str(round(datEdit)));
    errordlg('Input must be a positive integer','Error');
end
if datEdit >50
    set(hObject, 'String', 2);
    errordlg('Too large input number! Try a smaller one','Error');
end
set(handles.pb_SaveFitting, 'enable','off')
set(handles.pb_SaveFig,     'enable','off')
%
%--------------------------------------------------------------------------
%
% --- Executes during object creation, after setting all properties.
function edit_FIT_a4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_FIT_a4 (see GCBO)
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
function edit_FIT_a5_Callback(hObject, eventdata, handles)
% hObject    handle to edit_FIT_a5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fitting_level         = str2double(get(handles.edit_FIT_a4,'string'));
datEdit = str2double(get(hObject, 'String'));
if isnan(datEdit)
    set(hObject, 'String', 0);
    errordlg('Input must be a integer','Error');
end
if datEdit <0
    set(hObject, 'String', 0);
    errordlg('Input must be a nonnegative integer','Error');
end
if rem(datEdit ,1)~=0
    set(hObject, 'String', num2str(round(datEdit)));
    errordlg('Input must be a nonnegative integer','Error');
end
if datEdit >=fitting_level
    set(hObject, 'String', 0);
    errordlg('The relative degree must be smaller than the fitting order','Error');
end
set(handles.pb_SaveFitting, 'enable','off')
set(handles.pb_SaveFig,     'enable','off')
%
%--------------------------------------------------------------------------
%
% --- Executes during object creation, after setting all properties.
function edit_FIT_a5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_FIT_a5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%
%----------------------------end-------------------------------------------


% --------------------------------------------------------------------
function uipushtool1_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Fig = figure;
set(Fig,        'units','points')
posFig = get(handles.figure,'position');
copyobj(handles.axes1, Fig);
copyobj(handles.axes2, Fig);

hAxes = get(Fig,'children');
set(hAxes(1),       'units','points',...
                    'position',[80 60 200 150],...
                    'ActivePositionProperty','position')
set(hAxes(2),       'units','points',...
                    'position',[80 210 200 150],...
                    'ActivePositionProperty','position')
pos1=get(hAxes(1),'position');
pos2=get(hAxes(2),'position');
posAxesOuter = [0 0 310 400];
set(Fig,        'units','points',...
                'position', [posFig(1)+0.5*posFig(3)-0.5*posAxesOuter(3),...
                            posFig(2)+0.5*posFig(4)-0.5*posAxesOuter(4),...
                            posAxesOuter(3:4)])   


% --- Executes on selection change in pop_source.
function pop_source_Callback(hObject, eventdata, handles)
% hObject    handle to pop_source (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pop_source contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_source
% pop_type=get(hObject,'Value');
% switch pop_type
%     case 1
%         set(handles.edit_FIT_a3,        'enable','on');
%         set(handles.text_FIT_a4,        'visible','on'); 
%         set(handles.text_FIT_a5,        'visible','on'); 
%         set(handles.edit_FIT_a4,        'visible','on'); 
%         set(handles.edit_FIT_a5,        'visible','on');
%         set(handles.pb_sys_load,        'visible','off'); 
%         set(handles.pb_fitfrd,         'enable','on');
%     case 2
%         set(handles.edit_FIT_a3,        'enable','off');
%         set(handles.text_FIT_a4,        'visible','off'); 
%         set(handles.text_FIT_a5,        'visible','off'); 
%         set(handles.edit_FIT_a4,        'visible','off'); 
%         set(handles.edit_FIT_a5,        'visible','off');
%         set(handles.pb_sys_load,        'visible','on');
%         set(handles.pb_fitfrd,         'enable','off');
% end
Fcn_sys_source_pop_update(hObject);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function pop_source_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_source (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





% --- Executes on selection change in pop_plot.
function pop_plot_Callback(hObject, eventdata, handles)
% hObject    handle to pop_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pop_plot contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_plot
pop_type=get(hObject,'Value');
pannelsize=get(handles.uipanel_Axes,'position');
pW=pannelsize(3);
pH=pannelsize(4);   
switch pop_type
    case 1
        set(handles.axes1,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*2/10 pH*5.0/10 pW*6/10 pH*3.5/10],...
                        'fontsize',handles.FontSize(1),...
                        'color',handles.bgcolor{1},...
                        'box','on','visible','on');  
        set(handles.axes2,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*2/10 pH*1.5/10 pW*6/10 pH*3.5/10],...
                        'fontsize',handles.FontSize(1),...
                        'color',handles.bgcolor{1},...
                        'box','on','visible','on'); 
    case 2
        set(handles.axes1,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*2/10 pH*1.5/10 pW*7/10 pH*7/10],...
                        'fontsize',handles.FontSize(1),...
                        'color',handles.bgcolor{1},...
                        'box','on','visible','on');  
        set(handles.axes2,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*2/10 pH*1.5/10 pW*7/10 pH*3.5/10],...
                        'fontsize',handles.FontSize(1),...
                        'color',handles.bgcolor{1},...
                        'box','on','visible','off'); 
 end
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function pop_plot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_fRange_Callback(hObject, eventdata, handles)
% hObject    handle to edit_fRange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_fRange as text
%        str2double(get(hObject,'String')) returns contents of edit_fRange as a double


% --- Executes during object creation, after setting all properties.
function edit_fRange_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_fRange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
