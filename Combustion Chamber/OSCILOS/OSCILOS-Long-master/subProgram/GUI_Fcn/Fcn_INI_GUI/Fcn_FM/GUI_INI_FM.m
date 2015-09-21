function varargout = GUI_INI_FM(varargin)
% GUI_INI_FM MATLAB code for GUI_INI_FM.fig
%      GUI_INI_FM, by itself, creates a new GUI_INI_FM or raises the existing
%      singleton*.
%
%      H = GUI_INI_FM returns the handle to a new GUI_INI_FM or the handle to
%      the existing singleton*.
%
%      GUI_INI_FM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_INI_FM.M with the given input arguments.
%
%      GUI_INI_FM('Property','Value',...) creates a new GUI_INI_FM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_INI_FM_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_INI_FM_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_INI_FM

% Last Modified by GUIDE v2.5 15-Dec-2014 23:39:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_INI_FM_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_INI_FM_OutputFcn, ...
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


% --- Executes just before GUI_INI_FM is made visible.
function GUI_INI_FM_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_INI_FM (see VARARGIN)
global CI
mainGuiInput = 0;
% handles of main GUI
handles.MainGUI = varargin{mainGuiInput+1};
% Obtain handles using GUIDATA with the caller's handle 
mainHandles = guidata(handles.MainGUI);
% background colors
handles.bgcolor=mainHandles.bgcolor;
% fontsize
handles.FontSize=mainHandles.FontSize;
%
handles.sW = mainHandles.sW;
handles.sH = mainHandles.sH;
handles.indexApp = 0;
% --------------------------
%
handles.HP_num  = varargin{2};              % the index of unsteady heat source
handles.indexFM = varargin{3};              % the index of flame model
handles.HP = CI.FM.HP{handles.HP_num};  
guidata(hObject, handles);
%
% --------------------------
% Initialization
GUI_INI_FM_Initialization(hObject, eventdata, handles)
%         end
guidata(hObject, handles);
handles.output = hObject;
guidata(hObject, handles);

% 
%-------------------------------------------------
%
function GUI_INI_FM_Initialization(varargin)
hObject = varargin{1};
handles = guidata(hObject);    
set(0, 'units', 'points');
screenSize  = get(0, 'ScreenSize');                 % get the screen size
sW          = handles.sW;                           % screen width
sH          = handles.sH ;                          % screen height
FigW=sW*4/5;                                        % window width
FigH=sH*5/8;                                        % window height
set(handles.figure,     'units', 'points',...
                        'position',[(screenSize(3)-FigW)./2 (screenSize(4)-FigH)./2 FigW FigH],...
                        'name','Flame model configurations',...
                        'color',handles.bgcolor{3});
%--------------------------------------------------------------------------
% pannel axes
set(handles.uipanel_Axes,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[FigW*10.5/20 FigH*2.25/20 FigW*9/20 FigH*17.5/20],...
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
                        'position',[pW*2/10 pH*5.0/10 pW*6/10 pH*3.5/10],...
                        'fontsize',handles.FontSize(1),...
                        'color',handles.bgcolor{1},...
                        'box','on');  
set(handles.axes2,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*2/10 pH*1.5/10 pW*6/10 pH*3.5/10],...
                        'fontsize',handles.FontSize(1),...
                        'color',handles.bgcolor{1},...
                        'box','on');  
set(handles.pop_Plot,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*0.5/10 pH*8.5/10 pW*8/10 pH*1.0/10],...
                        'fontsize',handles.FontSize(2),...
                        'string',{  'Flame transfer function model';...
                                    'Nonlinear model';...
                                    'Flame describing function'},...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','left',...
                        'enable','on');                    
guidata(hObject, handles);
%----------------------------------------               
% pannels FTF
set(handles.uipanel_FTF,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[FigW*0.5/20 FigH*11/20 FigW*9.5/20 FigH*8.75/20],...
                        'Title','Linear flame transfer function model',...
                        'visible','on',...
                        'highlightcolor',handles.bgcolor{3},...
                        'borderwidth',1,...
                        'fontsize',handles.FontSize(2),...
                        'backgroundcolor',handles.bgcolor{3}); 
pannelsize = get(handles.uipanel_FTF,'position');
pW = pannelsize(3);
pH = pannelsize(4);                      
set(handles.pop_FTF_type,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*0.5/10 pH*8/10 pW*8.5/10 pH*1/10],...
                        'fontsize',handles.FontSize(2),...
                        'string',{  'n-tau model';...
                                    '1st-order low pass filter';...
                                    '2nd-order low pass filter';...
                                    'Transfer function model'},...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','left',...
                        'enable','on',...
                        'value',handles.HP.FTF.style );  
set(handles.text_FTF_a1,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*0.5/10 pH*5.75/10 pW*5/10 pH*1/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','af: [-]',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'horizontalalignment','left');                         
set(handles.edit_FTF_a1,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*6/10 pH*5.75/10 pW*3/10 pH*1.25/10],...
                        'fontsize',handles.FontSize(2),...
                        'string',1,...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','right',...
                        'Enable','off');

set(handles.text_FTF_a2,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*0.5/10 pH*4.0/10 pW*5/10 pH*1.0/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','fc: [Hz]',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'horizontalalignment','left');                   
set(handles.edit_FTF_a2,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*6/10 pH*4/10 pW*3/10 pH*1.25/10],...
                        'fontsize',handles.FontSize(2),...
                        'string',100,...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','right',...
                        'Enable','off');  
set(handles.text_FTF_a3,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*0.5/10 pH*2.25/10 pW*5/10 pH*1.0/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Damping ratio: [-]',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'horizontalalignment','left');                   
set(handles.edit_FTF_a3,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*6/10 pH*2.25/10 pW*3/10 pH*1.25/10],...
                        'fontsize',handles.FontSize(2),...
                        'string',0.5,...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','right',...
                        'Enable','off'); 
set(handles.text_FTF_a4,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*0.5/10 pH*0.5/10 pW*5/10 pH*1.0/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','tau_f: [ms]',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'horizontalalignment','left');                   
set(handles.edit_FTF_a4,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*6/10 pH*0.5/10 pW*3/10 pH*1.25/10],...
                        'fontsize',handles.FontSize(2),...
                        'string',3,...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','right',...
                        'Enable','off');  
%----------------------------------------
% pannels FDF
set(handles.uipanel_FDF,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[FigW*0.5/20 FigH*2.25/20 FigW*9.5/20 FigH*8.5/20],...
                        'Title','Nonlinear model',...
                        'visible','on',...
                        'highlightcolor',handles.bgcolor{3},...
                        'borderwidth',1,...
                        'fontsize',handles.FontSize(2),...
                        'backgroundcolor',handles.bgcolor{3}); 
pannelsize=get(handles.uipanel_FDF,'position');
pW=pannelsize(3);
pH=pannelsize(4);                      
set(handles.pop_NL_type,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*0.5/10 pH*8/10 pW*8.5/10 pH*1/10],...
                        'fontsize',handles.FontSize(2),...
                        'string',{  'None';...
                                    'Dowling''s model (Stow and Dowling JEGTP 2009)';...
                                    'J.Li and A.S.Morgans JSV'},...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','left',...
                        'enable','on',...
                        'value',handles.HP.NL.style);                      
set(handles.text_NL_a1,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*0.5/10 pH*5/10 pW*5/10 pH*1.0/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','alpha: [-]',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'horizontalalignment','left');                         
set(handles.edit_NL_a1,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*6/10 pH*5/10 pW*3/10 pH*1.25/10],...
                        'fontsize',handles.FontSize(2),...
                        'string',0.8,...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','right',...
                        'Enable','off');

set(handles.text_NL_a2,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*0.5/10 pH*3/10 pW*5/10 pH*1.0/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','beta: [-]',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'horizontalalignment','left');                   
set(handles.edit_NL_a2,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*6/10 pH*3/10 pW*3/10 pH*1.25/10],...
                        'fontsize',handles.FontSize(2),...
                        'string',30,...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','right',...
                        'Enable','off');  
set(handles.text_NL_a3,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*0.5/10 pH*1/10 pW*5/10 pH*1/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','taufN: [ms]',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'horizontalalignment','left');                   
set(handles.edit_NL_a3,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*6/10 pH*1/10 pW*3/10 pH*1.25/10],...
                        'fontsize',handles.FontSize(2),...
                        'string',0,...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','right',...
                        'Enable','off');  
%----------------------------------------
%
% pannel AOC                   
set(handles.uipanel_AOC,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[FigW*0.5/20 FigH*0/20 FigW*19/20 FigH*1.8/20],...
                        'Title','',...
                        'visible','on',...
                        'highlightcolor',handles.bgcolor{3},...
                        'borderwidth',1,...
                        'fontsize',handles.FontSize(2),...
                        'backgroundcolor',handles.bgcolor{3}); 
pannelsize=get(handles.uipanel_AOC,'position');                    
pW=pannelsize(3);
pH=pannelsize(4);                
set(handles.pb_Apply,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*1/10 pH*2/10 pW*2.0/10 pH*6/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Plot figure',...
                        'backgroundcolor',handles.bgcolor{3});
set(handles.pb_SaveFig,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*4/10 pH*2/10 pW*2.0/10 pH*6/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Save figure',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'enable','off',...
                        'visible','off');
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
%---------------------------------------
handles.ObjVisible_FTF      = findobj('-regexp','Tag','FTF');
handles.ObjVisible_NL       = findobj('-regexp','Tag','NL');
handles.objVisible_edit_NL  = findobj('-regexp','Tag','edit_NL');
handles.objVisible_text_NL  = findobj('-regexp','Tag','text_NL');
% default visible settings
set(handles.ObjVisible_FTF,         'visible','on')
set(handles.ObjVisible_NL,          'visible','on')
%
handles.ObjEditEnable_FTF   = findobj('-regexp','Tag','edit_FTF');
handles.ObjEditEnable_NL    = findobj('-regexp','Tag','edit_NL');
% default enable settings
set(handles.ObjEditEnable_FTF,      'Enable','on')
set(handles.ObjEditEnable_NL,       'Enable','on')

%---------------------------------------
guidata(hObject, handles);
Fcn_GUI_INI_FM_FTF_POP_update(hObject)
handles = guidata(hObject);
guidata(hObject, handles);
switch handles.indexFM
    case 1
        set(handles.uipanel_FDF,          'visible','off')
    case 2
        set(handles.uipanel_FDF,          'visible','on')
end
guidata(hObject, handles);
        Fcn_GUI_INI_FM_NL_POP_update(hObject)
        handles = guidata(hObject);
guidata(hObject, handles);
% 
% ------------------------------------------------------------------------
%
function Fcn_GUI_INI_FM_FTF_POP_update(varargin)
hObject         = varargin{1};
handles         = guidata(hObject);
pop_FTF_type    = get(handles.pop_FTF_type,'Value');
%
switch pop_FTF_type
    case 1
        set(handles.edit_FTF_a1,                'visible','on',...
                                                'string', num2str(handles.HP.FTF.af));
        set(handles.text_FTF_a1,                'visible','on',...
                                                'string','af [-]'); 
        set(handles.edit_FTF_a2,                'visible','off');
        set(handles.text_FTF_a2,                'visible','off');
        set(handles.edit_FTF_a3,                'visible','off');
        set(handles.text_FTF_a3,                'visible','off');
        set(handles.edit_FTF_a4,                'visible','on',...
                                                'string', num2str(handles.HP.FTF.tauf.*1e3));
        set(handles.text_FTF_a4,                'visible','on',...
                                                'string','tauf [ms]');
    case 2
        set(handles.edit_FTF_a1,                'visible','on',...
                                                'string', num2str(handles.HP.FTF.af));
        set(handles.text_FTF_a1,                'visible','on',...
                                                'string','af [-]'); 
        set(handles.edit_FTF_a2,                'visible','on',...
                                                'string', num2str(handles.HP.FTF.fc));
        set(handles.text_FTF_a2,                'visible','on',...
                                                'string','fc [Hz]'); 
        set(handles.edit_FTF_a3,                'visible','off');
        set(handles.text_FTF_a3,                'visible','off');
        set(handles.edit_FTF_a4,                'visible','on',...
                                                'string', num2str(handles.HP.FTF.tauf.*1e3));
        set(handles.text_FTF_a4,                'visible','on',...
                                                'string','tauf [ms]'); 
    case 3
        set(handles.edit_FTF_a1,                'visible','on',...
                                                'string', num2str(handles.HP.FTF.af));
        set(handles.text_FTF_a1,                'visible','on',...
                                                'string','af [-]'); 
        set(handles.edit_FTF_a2,                'visible','on',...
                                                'string', num2str(handles.HP.FTF.fc));
        set(handles.text_FTF_a2,                'visible','on',...
                                                'string','fc [Hz]'); 
        set(handles.edit_FTF_a3,                'visible','on',...
                                                'string', num2str(handles.HP.FTF.xi));
        set(handles.text_FTF_a3,                'visible','on',...
                                                'string','xi [-]');
        set(handles.edit_FTF_a4,                'visible','on',...
                                                'string', num2str(handles.HP.FTF.tauf.*1e3));
        set(handles.text_FTF_a4,                'visible','on',...
                                                'string','tauf [ms]');  
    case 4
        handles.indexApp = Fcn_GUI_INI_FM_helpdlg(handles.indexApp);
        strNum = ['[' num2str(handles.HP.FTF.num) ']'];
        strDen = ['[' num2str(handles.HP.FTF.den) ']'];
        set(handles.edit_FTF_a1,                'visible','on',...
                                                'string', strNum);
        set(handles.text_FTF_a1,                'visible','on',...
                                                'string', 'Numerator(s):');
        set(handles.edit_FTF_a2,                'visible','on',...
                                                'string', strDen);
        set(handles.text_FTF_a2,                'visible','on',...
                                                'string', 'Denominator(s):');
        set(handles.edit_FTF_a3,                'visible','off');
        set(handles.text_FTF_a3,                'visible','off');
        set(handles.edit_FTF_a4,                'visible','on',...
                                                'string', num2str(handles.HP.FTF.tauf.*1e3));
        set(handles.text_FTF_a4,                'visible','on',...
                                                'string','tauf [ms]'); 
end
guidata(hObject, handles);
%
% ------------------------------------------------------------------------
%
function Fcn_GUI_INI_FM_NL_POP_update(varargin)
hObject         = varargin{1};
handles         = guidata(hObject);
pop_NL_type    = get(handles.pop_NL_type,'Value');
%
switch pop_NL_type
    case 1
        set(handles.pop_Plot,                   'enable','off');
        set(handles.text_NL_a1,                 'visible','off');
        set(handles.text_NL_a2,                 'visible','off');
        set(handles.text_NL_a3,                 'visible','off');
        set(handles.edit_NL_a1,                 'visible','off');
        set(handles.edit_NL_a2,                 'visible','off');
        set(handles.edit_NL_a3,                 'visible','off');
    case 2
        set(handles.pop_Plot,                  'enable','on');
        set(handles.text_NL_a2,                 'visible','off');
        set(handles.text_NL_a3,                 'visible','off');
        set(handles.edit_NL_a2,                 'visible','off');
        set(handles.edit_NL_a3,                 'visible','off');
        set(handles.edit_NL_a1,                'visible','on',...
                                               'string', num2str(handles.HP.NL.Model2.alpha));
        set(handles.text_NL_a1,                'visible','on',...
                                               'string','alpha [-]'); 
    case 3
        set(handles.pop_Plot,                  'enable','on');
        set(handles.edit_NL_a1,                'visible','on',...
                                               'string', num2str(handles.HP.NL.Model3.alpha));
        set(handles.text_NL_a1,                'visible','on',...
                                               'string','alpha [-]'); 
        set(handles.edit_NL_a2,                'visible','on',...
                                               'string', num2str(handles.HP.NL.Model3.beta));
        set(handles.text_NL_a2,                'visible','on',...
                                               'string','beta [-]'); 
        set(handles.edit_NL_a3,                'visible','on',...
                                               'string', num2str(handles.HP.NL.Model3.taufN.*1e3));
        set(handles.text_NL_a3,                'visible','on',...
                                               'string','taufN [ms]'); 
end
guidata(hObject, handles);
%
function Fcn_GUI_INI_FM_FTF_Value_update(varargin)
hObject = varargin{1};
handles = guidata(hObject);
% global CI
handles.HP.FTF.style = get(handles.pop_FTF_type,'Value');
%
switch handles.HP.FTF.style
    case 1      
        handles.HP.FTF.af        = str2double(get(handles.edit_FTF_a1,'String'));
        handles.HP.FTF.tauf      = str2double(get(handles.edit_FTF_a4,'String')).*1e-3;
        handles.HP.FTF.num       = handles.HP.FTF.af;
        handles.HP.FTF.den       = 1;
    case 2      
        handles.HP.FTF.af        = str2double(get(handles.edit_FTF_a1,'String'));
        handles.HP.FTF.tauf      = str2double(get(handles.edit_FTF_a4,'String')).*1e-3;
        handles.HP.FTF.fc        = str2double(get(handles.edit_FTF_a2,'String'));
        handles.HP.FTF.omegac    = 2*pi*handles.HP.FTF.fc;
        handles.HP.FTF.num       = handles.HP.FTF.af.*handles.HP.FTF.omegac;
        handles.HP.FTF.den       = [1 handles.HP.FTF.omegac];
    case 3      %
        handles.HP.FTF.af        = str2double(get(handles.edit_FTF_a1,'String'));
        handles.HP.FTF.tauf      = str2double(get(handles.edit_FTF_a4,'String')).*1e-3;
        handles.HP.FTF.fc        = str2double(get(handles.edit_FTF_a2,'String'));
        handles.HP.FTF.xi        = str2double(get(handles.edit_FTF_a3,'String'));
        handles.HP.FTF.omegac    = 2*pi*handles.HP.FTF.fc;
        handles.HP.FTF.num       = handles.HP.FTF.af.*handles.HP.FTF.omegac.^2;
        handles.HP.FTF.den       = [1 2*handles.HP.FTF.xi.*handles.HP.FTF.omegac handles.HP.FTF.omegac.^2];
        %
     case 4   
        handles.HP.FTF.num       = str2num(get(handles.edit_FTF_a1,'String'));
        handles.HP.FTF.den       = str2num(get(handles.edit_FTF_a2,'String'));
        handles.HP.FTF.tauf      = str2num(get(handles.edit_FTF_a4, 'String')).*1e-3; 
end
% assignin('base','CI',CI);                   % save the current information to the workspace
guidata(hObject, handles);
% ------------------------------------------------------------------------
function Fcn_GUI_INI_FM_NL_Value_update(varargin)
hObject = varargin{1};
handles = guidata(hObject);
% global CI
handles.HP.NL.style = get(handles.pop_NL_type,'Value');
%
switch handles.HP.NL.style
    case 1
    case 2
        handles.HP.NL.Model2.alpha       = str2double(get(handles.edit_NL_a1,'String'));
        DqRatio                          = 1e-3;
        qRatioMax                        = 1;
        [handles.HP.NL.Model2.qRatioLinear,handles.HP.NL.Model2.Lf]...
       = Fcn_GUI_INI_FM_Nonlinear_model_Dowling(handles.HP.NL.Model2.alpha,...
                                                DqRatio,qRatioMax);
    case 3      
        handles.HP.NL.Model3.alpha       = str2double(get(handles.edit_NL_a1,'String'));
        handles.HP.NL.Model3.beta        = str2double(get(handles.edit_NL_a2,'String'));
        handles.HP.NL.Model3.taufN       = str2double(get(handles.edit_NL_a3,'String')).*1e-3;
        DuRatio                     = 1e-3;
        uRatioMax                   = 5;
        [handles.HP.NL.Model3.uRatio,handles.HP.NL.Model3.Lf]...
      = Fcn_GUI_INI_FM_Nonlinear_model_J_Li_A_Morgas(   handles.HP.NL.Model3.alpha,...
                                                        handles.HP.NL.Model3.beta,...
                                                        DuRatio,uRatioMax);
end
% assignin('base','CI',CI);                   % save the current information to the workspace
guidata(hObject, handles);
% ------------------------------------------------------------------------
%
function indexApp = Fcn_GUI_INI_FM_helpdlg(indexApp)
if indexApp < 2
string={...
'The numerator coefficient can be a vector or matrix expression.';...
'The denominator coefficient must be a vector.';...
'The output width equals the number of rows in the numerator coefficient.';
'You should specify the coefficients in descending order of powers of s.'};
helpdlg(string,'Set transfer function coefficients')
end
indexApp = indexApp + 1;
%
% ------------------------------------------------------------------------
%
% --- Executes during object creation, after setting all properties.
function pop_FTF_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_FTF_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes on selection change in pop_FTF_type.
function pop_FTF_type_Callback(hObject, eventdata, handles)
% hObject    handle to pop_FTF_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.pb_SaveFig,'enable','off');
Fcn_GUI_INI_FM_FTF_POP_update(hObject)
handles = guidata(hObject);
guidata(hObject, handles);

%--------------------------------------------------------------------------
% --- Executes on selection change in pop_NL_type.
function pop_NL_type_Callback(hObject, eventdata, handles)
% hObject    handle to pop_NL_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = cellstr(get(hObject,'String')) returns pop_NL_type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_NL_type
set(handles.pb_SaveFig,'enable','off');
Fcn_GUI_INI_FM_NL_POP_update(hObject)
handles = guidata(hObject);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function pop_NL_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_NL_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%
% --- Executes on button press in pb_OK.
function pb_OK_Callback(hObject, eventdata, handles)
% hObject    handle to pb_OK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Fcn_GUI_INI_FM_FTF_Value_update(hObject, eventdata, handles)
handles = guidata(hObject);
guidata(hObject, handles);
Fcn_GUI_INI_FM_NL_Value_update(hObject, eventdata, handles)
handles = guidata(hObject);
guidata(hObject, handles);
Fcn_GUI_INI_FM_Plot(hObject, eventdata, handles)
handles = guidata(hObject);
handles.HP.IsRun = 1;    % index == 1 to show that this program has ever been run
guidata(hObject, handles);
global CI
CI.FM.HP{handles.HP_num}        = handles.HP;
CI.FM.indexFM(handles.HP_num)   = handles.indexFM;
assignin('base','CI',CI); 
handles = guidata(hObject);
guidata(hObject, handles);
Fcn_GUI_INI_FM_Update_Data(hObject, eventdata, handles);  % put it before % delete
delete(handles.figure);
%
% --- Executes on button press in pb_Apply.
function pb_Apply_Callback(hObject, eventdata, handles)
% hObject    handle to pb_Apply (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Fcn_GUI_INI_FM_FTF_Value_update(hObject, eventdata, handles)
Fcn_GUI_INI_FM_NL_Value_update(hObject, eventdata, handles)
Fcn_GUI_INI_FM_Plot(hObject, eventdata, handles)
handles = guidata(hObject);
set(handles.pb_SaveFig,'enable','on');
guidata(hObject, handles);
%
% --- Outputs from this function are returned to the command line.
function varargout = GUI_INI_FM_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
varargout{1} = handles.output;
end
%
%-------------------------------------------------------------------------
function Fcn_GUI_INI_FM_Plot(varargin)
hObject     = varargin{1};
handles     = guidata(hObject);
global CI
hAxes1      = handles.axes1;
hAxes2      = handles.axes2;
fontSize1   = handles.FontSize(1);
fontSize2   = handles.FontSize(2);
fs          = 0:5:1000;              
s           = 2*pi*fs*1i;
pop_PLot    = get(handles.pop_Plot,'Value');
%------------------------------------
FTF = polyval(handles.HP.FTF.num,s)./polyval(handles.HP.FTF.den,s).*exp(-s.*handles.HP.FTF.tauf);
uRatio = 0:0.1:1;
switch handles.HP.NL.style       
    case 1

    case 2
        for ss = 1:length(fs)
            for kk = 1 :length(uRatio)
                qRatioLinear    = abs(FTF(ss)).*uRatio(kk);
                LfModel2(kk,ss) = interp1(  handles.HP.NL.Model2.qRatioLinear,...
                                            handles.HP.NL.Model2.Lf,...
                                            qRatioLinear,'linear','extrap');
                FDF(kk,ss)      = LfModel2(kk,ss).*FTF(ss);
            end
        end 
    case 3      %  J.Li and A.S.Morgans JSV
        LfModel3        = interp1(  handles.HP.NL.Model3.uRatio,...
                                    handles.HP.NL.Model3.Lf,...
                                    uRatio,'linear','extrap');
        taufNSp         = handles.HP.NL.Model3.taufN.*(1-LfModel3);
        for kk = 1 :length(uRatio) 
        FDF(kk,:)       = FTF.*LfModel3(kk).*exp(-s.*taufNSp(kk));
        end
end
assignin('base','CI',CI);                   
%--------------------------------------------------------------------------
switch pop_PLot
    case 1
        colorPlot={'b','g'};
        cla(hAxes1,'reset')
        axes(hAxes1)
        hold on
        hLine11     = plot(hAxes1,fs,abs(FTF),'-','color',colorPlot{1},...
                        'Linewidth',2);
        set(hAxes1,'YColor','k','Box','on','ygrid','on','xgrid','on');
        set(hAxes1,'FontName','Helvetica','FontSize',fontSize1,'LineWidth',1)
        xlabel(hAxes1,'','Color','k','Interpreter','LaTex','FontSize',fontSize1);
        ylabel(hAxes1,'Gain [-]','Color','k','Interpreter','LaTex',...
            'FontSize',fontSize1)
        set(hAxes1,'xlim',[0 1000],'xTick',200:200:1000,'xticklabel',{},...
            'YAxisLocation','left','Color','w');
        hold off
        %--------------------------------
        cla(hAxes2,'reset')
        axes(hAxes2)
        hold on
        hLine21     = plot(hAxes2,fs,unwrap(angle(FTF),1.9*pi)./pi,'-',...
                        'color', colorPlot{1},'Linewidth',2);
        set(hAxes2,'YColor','k','Box','on','ygrid','on','xgrid','on');
        set(hAxes2,'FontName','Helvetica','FontSize',fontSize1,'LineWidth',1)
        xlabel(hAxes2,'$f$ [Hz]','Color','k','Interpreter','LaTex',...
            'FontSize',fontSize1);
        ylabel(hAxes2,'Phase/$\pi$ [-]','Color','k','Interpreter','LaTex',...
            'FontSize',fontSize1)
        set(hAxes2,'xlim',get(hAxes1,'xlim'),'xTick',get(hAxes1,'xTick'),...
        'YAxisLocation','left','Color','w');
        hold off
        try
            cbh = findobj( 0, 'tag', 'Colorbar' );
            delete(cbh)
        catch
        end
    case 2
        colorPlot={'r','g'};
        cla(hAxes1,'reset')
        axes(hAxes1)
        hold on
        set(hAxes1,'YColor','k','Box','on','ygrid','on','xgrid','on');
        set(hAxes1,'FontName','Helvetica','FontSize',fontSize1,'LineWidth',1)
        xlabel(hAxes1,'','Color','k','Interpreter','LaTex','FontSize',fontSize1);
        set(hAxes1,'xticklabel',{},'YAxisLocation','left','Color','w');
        switch handles.HP.NL.style   
            case 1
                errordlg('Nonlinearity is not accounted for!!!','Error');
            case 2
            plot(hAxes1,handles.HP.NL.Model2.qRatioLinear,...
                        handles.HP.NL.Model2.Lf.*handles.HP.NL.Model2.qRatioLinear,...
                        '-','color',colorPlot{1},...
                        'Linewidth',2);
            ylabel(hAxes1,'$\hat{\dot{q}}/\bar{\dot{q}}_{Nonlinear}$ [-]','Color','k',...
                'Interpreter','LaTex','FontSize',fontSize1)
            case 3
            plot(hAxes1,handles.HP.NL.Model3.uRatio,...
                        handles.HP.NL.Model3.Lf.*handles.HP.NL.Model3.uRatio,...
                        '-','color',colorPlot{1},'Linewidth',2);
            ylabel(hAxes1,'$\hat{\dot{q}}/\bar{\dot{q}}(f = 0)$ [-]','Color','k',...
                'Interpreter','LaTex','FontSize',fontSize1)
            set(hAxes1,'xlim',[0 1],'xtick',0:0.2:1)
        end
        hold off  
        %--------------------------------
        cla(hAxes2,'reset')
        axes(hAxes2)
        hold on
        set(hAxes2,'YColor','k','Box','on','ygrid','on','xgrid','on');
        set(hAxes2,'FontName','Helvetica','FontSize',fontSize1,'LineWidth',1)
        set(hAxes2,'YAxisLocation','left','Color','w');
        switch handles.HP.NL.style   
            case 1
                errordlg('Nonlinearity is not accounted for!!!','Error');
            case 2
                plot(hAxes2,handles.HP.NL.Model2.qRatioLinear,...
                            handles.HP.NL.Model2.Lf,...
                        '-','color',colorPlot{1},...
                        'Linewidth',2);
                xlabel(hAxes2,'$\hat{\dot{q}}/\bar{\dot{q}}_{Linear}$ [-]','Color','k','Interpreter',...
                'LaTex','FontSize',fontSize1);
                ylabel(hAxes2,'$L_f$ [-]','Color','k','Interpreter','LaTex',...
                'FontSize',fontSize1)
            case 3
                plot(hAxes2,handles.HP.NL.Model3.uRatio,...
                            handles.HP.NL.Model3.Lf,...
                            '-','color',colorPlot{1},'Linewidth',2);
                xlabel(hAxes2,'$\hat{u}_u/\bar{u}_u$ [-]','Color','k','Interpreter',...
                'LaTex','FontSize',fontSize1);
                ylabel(hAxes2,'$L_f$ [-]','Color','k','Interpreter','LaTex',...
                'FontSize',fontSize1)
        end
        set(hAxes2,'xlim',get(hAxes1,'xlim'),'xTick',get(hAxes1,'xTick'),...
            'YAxisLocation','left','Color','w');
        hold off
        try
            cbh = findobj( 0, 'tag', 'Colorbar' );
            delete(cbh)
        catch
        end
    case 3
        try
            cbh = findobj( 0, 'tag', 'Colorbar' );
            delete(cbh)
        catch
        end
        color_type      = Fcn_color_set(length(uRatio)+2);    % color set
        cla(hAxes1,'reset')
        axes(hAxes1)
        switch handles.HP.NL.style   
            case 1
                errordlg('Nonlinearity is not accounted for!!!','Error');
            otherwise
                hold on
                for ss=1:length(uRatio)
                    hLine11 = plot(hAxes1,fs,abs(FDF(ss,:)),...
                        '-','color',color_type(ss,1:3),'Linewidth',2);
                end
        end
        set(hAxes1,'YColor','k','Box','on','ygrid','on','xgrid','on');
        set(hAxes1,'FontName','Helvetica','FontSize',fontSize1,'LineWidth',1)
        xlabel(hAxes1,'','Color','k','Interpreter','LaTex','FontSize',fontSize1);
        ylabel(hAxes1,'Gain [-]','Color','k','Interpreter','LaTex',...
            'FontSize',fontSize1)
        set(hAxes1,'xlim',[0 1000],'xTick',200:200:1000,'xticklabel',{},...
        'YAxisLocation','left','Color','w');
        hold off
        %--------------------------------
        cla(hAxes2,'reset')
        axes(hAxes2)
        hold on
        switch handles.HP.NL.style   
            case 1
                errordlg('Nonlinearity is not accounted for!!!','Error');
            otherwise
                for ss=1:length(uRatio)
                    plot(hAxes2,fs,(unwrap(angle(FDF(ss,:)),1.9*pi))./pi,...
                        '-','color', color_type(ss,1:3),'Linewidth',2);
                end
        end
            set(hAxes2,'YColor','k','Box','on','ygrid','on','xgrid','on');
            set(hAxes2,'FontName','Helvetica','FontSize',fontSize1,'LineWidth',1)
            xlabel(hAxes2,'$f$ [Hz]','Color','k','Interpreter','LaTex',...
                'FontSize',fontSize1);
            ylabel(hAxes2,'Phase/$\pi$ [-]','Color','k','Interpreter','LaTex',...
                'FontSize',fontSize1)
            set(hAxes2,'xlim',get(hAxes1,'xlim'),'xTick',get(hAxes1,'xTick'),...
            'YAxisLocation','left','Color','w');
        hold off
        Pos_hAxes2      = get(hAxes2,'position');
        colormap(hot);
        z_range         = uRatio;
        hcb             = colorbar;
        set(hcb,'parent',handles.uipanel_Axes)
        hcb_ylim        = get(hcb,'ylim');
        hcb_ytick_conv  = 0:max(z_range)/5:max(z_range);
        hcb_ytick       = hcb_ytick_conv.*...
            ((max(hcb_ylim)-min(hcb_ylim))./(max(z_range)-min(z_range)));
        set(hcb,'ytick',hcb_ytick,'yticklabel',hcb_ytick_conv);
        set(hcb,'Fontsize',fontSize1,'box','on','Unit','points')
        set(hcb,'position',[Pos_hAxes2(1)+Pos_hAxes2(3),...
                            Pos_hAxes2(2),...
                            Pos_hAxes2(3)./20,...
                            Pos_hAxes2(4).*2]);
        set(hAxes2,'position',Pos_hAxes2)
        handles.hColorbar.style     = 'hot';
        handles.hColorbar.ylimit    = [0,max(z_range)];
        handles.hColorbar.ytick     = hcb_ytick;
        handles.hColorbar.ytick_conv= hcb_ytick_conv;
        %------------------------
        xt_pos = min((get(hAxes2,'xlim')))+1.25*(max((get(hAxes2,'xlim')))...
            -min((get(hAxes2,'xlim'))));
        yt_pos = min(get(hAxes1,'ylim'));
        hcbText=text(xt_pos,yt_pos,'$\hat{u}_u^{}/\bar{u}_u^{}$ [-]',...
            'FontSize',fontSize1, 'interpreter','latex','rotation',90,...
            'HorizontalAlignment','center');
        set(hcbText,'parent',hAxes2)
end
%--------------------------------
guidata(hObject,handles);
assignin('base','CI',CI);                   % save the current information to the workspace
%
% --- Update the data when clicking 'OK' or 'Apply'
function Fcn_GUI_INI_FM_Update_Data(hObject, eventdata, handles)
handles = guidata(hObject);
global CI
main = handles.MainGUI;
if(ishandle(main))
    mainHandles = guidata(main);
    changeMain = mainHandles.uitable1;
    table_cell = get(changeMain, 'data');
    table_cell{handles.HP_num,1}= CI.FM.ModelType{handles.indexFM};
    table_cell{handles.HP_num,2}= 'Y';
    set(changeMain,'data',table_cell);  % change the table 
    % -----------
    N = length(CI.FM.indexFM);
    for ss = 1:N
        HP = CI.FM.HP{ss};
        isRunSp(ss) = HP.IsRun;
    end
    if isRunSp == 1
        changeMain2 = mainHandles.pb_OK;
        set(changeMain2, 'enable','on');
    end            
end
guidata(hObject, handles);
%
% --- Executes on button press in pb_SaveFig.
function pb_SaveFig_Callback(varargin)
% global CI
hObject = varargin{1};
handles = guidata(hObject);
pop_PLot            = get(handles.pop_Plot,'Value');
pop_FTF_type        = get(handles.pop_FTF_type,'Value');
Fig                 = figure;
set(Fig,'units','points')
posFig              = get(handles.figure,'position');
copyobj(handles.axes1, Fig);
copyobj(handles.axes2, Fig);

hAxes               = get(Fig,'children');
set(hAxes(1),       'units','points',...
                    'position',[60 60 200 150],...
                    'ActivePositionProperty','position')
set(hAxes(2),       'units','points',...
                    'position',[60 210 200 150],...
                    'ActivePositionProperty','position')
pos1                = get(hAxes(1),'position');
pos2                = get(hAxes(2),'position');
switch pop_PLot
    case 3
        try
            posAxesOuter = [0 0 340 400];
            colormap(hot);
            hcb = colorbar;
            set(hcb,'parent',Fig)
            set(hcb,'ytick',handles.hColorbar.ytick,...
                    'yticklabel',handles.hColorbar.ytick_conv);
            set(hcb,'Fontsize',handles.FontSize(1),'box','on','Unit','points')
            set(hcb,'ActivePositionProperty','position');
            set(hcb,'position',[260,60,10,300]);
            set(hAxes(1),       'position',pos1)
            set(hAxes(2),       'position',pos2) 
            set(hcb,'position',[260,60,10,300]);
        catch
        end
    otherwise
        posAxesOuter = [0 0 300 400];
end
set(Fig,        'units','points',...
                'position', [posFig(1)+0.5*posFig(3)-0.5*posAxesOuter(3),...
                            posFig(2)+0.5*posFig(4)-0.5*posAxesOuter(4),...
                            posAxesOuter(3:4)])       

% --- Executes on button press in pb_Cancel.
function pb_Cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pb_Cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(handles.figure);
%
% --- Executes when user attempts to close figure.
function figure_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(hObject);
%
% --- Executes on selection change in pop_Plot.
function pop_Plot_Callback(hObject, eventdata, handles)
% hObject    handle to pop_Plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.pb_SaveFig,'enable','off');
%
% --- Executes during object creation, after setting all properties.
function pop_Plot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_Plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%
function edit_NL_a1_Callback(hObject, eventdata, handles)
% hObject    handle to edit_NL_a1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_NL_a1 as text
%        str2double(get(hObject,'String')) returns contents of edit_NL_a1 as a double
datEdit = str2double(get(hObject, 'String'));
if isnan(datEdit)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end
%
% --- Executes during object creation, after setting all properties.
function edit_NL_a1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_NL_a1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%
function edit_NL_a2_Callback(hObject, eventdata, handles)
% hObject    handle to edit_NL_a2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_NL_a2 as text
%        str2double(get(hObject,'String')) returns contents of edit_NL_a2 as a double
datEdit = str2double(get(hObject, 'String'));
if isnan(datEdit)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end
%
% --- Executes during object creation, after setting all properties.
function edit_NL_a2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_NL_a2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%
function edit_NL_a3_Callback(hObject, eventdata, handles)
% hObject    handle to edit_NL_a3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_NL_a3 as text
%        str2double(get(hObject,'String')) returns contents of edit_NL_a3 as a double
datEdit = str2double(get(hObject, 'String'));
if isnan(datEdit)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end

% --- Executes during object creation, after setting all properties.
function edit_NL_a3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_NL_a3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%
function edit_FTF_a1_Callback(hObject, eventdata, handles)
% hObject    handle to edit_FTF_a1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_FTF_a1 as text
%        str2double(get(hObject,'String')) returns contents of edit_FTF_a1 as a double
datEdit = str2num(get(hObject, 'String'));
if isnan(datEdit)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end
%
% --- Executes during object creation, after setting all properties.
function edit_FTF_a1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_FTF_a1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_FTF_a2_Callback(hObject, eventdata, handles)
% hObject    handle to edit_FTF_a2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_FTF_a2 as text
%        str2double(get(hObject,'String')) returns contents of edit_FTF_a2 as a double
datEdit = str2num(get(hObject, 'String'));
if isnan(datEdit)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end

% --- Executes during object creation, after setting all properties.
function edit_FTF_a2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_FTF_a2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_FTF_a3_Callback(hObject, eventdata, handles)
% hObject    handle to edit_FTF_a3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_FTF_a3 as text
%        str2double(get(hObject,'String')) returns contents of edit_FTF_a3 as a double
datEdit = str2num(get(hObject, 'String'));
if isnan(datEdit)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end

% --- Executes during object creation, after setting all properties.
function edit_FTF_a3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_FTF_a3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_FTF_a4_Callback(hObject, eventdata, handles)
% hObject    handle to edit_FTF_a4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_FTF_a4 as text
%        str2double(get(hObject,'String')) returns contents of edit_FTF_a4 as a double
datEdit = str2num(get(hObject, 'String'));
if isnan(datEdit)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end

% --- Executes during object creation, after setting all properties.
function edit_FTF_a4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_FTF_a4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%--------------------------------------------------------------------------
%
function color_type=Fcn_color_set(color_number)
cmap = colormap(hot);
n_sample=round(linspace(1,64,color_number));
for ss=1:length(n_sample)
    color_type(ss,1:3)=cmap(n_sample(ss),1:3);
end
%---------------------------end-------------------------------------------


% --------------------------------------------------------------------
function uipushtool1_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% hObject = varargin{1};
% handles = guidata(hObject);
pop_PLot            = get(handles.pop_Plot,'Value');
pop_FTF_type        = get(handles.pop_FTF_type,'Value');
Fig                 = figure;
set(Fig,'units','points')
posFig              = get(handles.figure,'position');
copyobj(handles.axes1, Fig);
copyobj(handles.axes2, Fig);

hAxes               = get(Fig,'children');
set(hAxes(1),       'units','points',...
                    'position',[60 60 200 150],...
                    'ActivePositionProperty','position')
set(hAxes(2),       'units','points',...
                    'position',[60 210 200 150],...
                    'ActivePositionProperty','position')
pos1                = get(hAxes(1),'position');
pos2                = get(hAxes(2),'position');
switch pop_PLot
    case 3
        try
            posAxesOuter = [0 0 340 400];
            colormap(hot);
            hcb = colorbar;
            set(hcb,'parent',Fig)
            set(hcb,'ytick',handles.hColorbar.ytick,...
                    'yticklabel',handles.hColorbar.ytick_conv);
            set(hcb,'Fontsize',handles.FontSize(1),'box','on','Unit','points')
            set(hcb,'ActivePositionProperty','position');
            set(hcb,'position',[260,60,10,300]);
            set(hAxes(1),       'position',pos1)
            set(hAxes(2),       'position',pos2) 
            set(hcb,'position',[260,60,10,300]);
        catch
        end
    otherwise
        posAxesOuter = [0 0 300 400];
end
set(Fig,        'units','points',...
                'position', [posFig(1)+0.5*posFig(3)-0.5*posAxesOuter(3),...
                            posFig(2)+0.5*posFig(4)-0.5*posAxesOuter(4),...
                            posAxesOuter(3:4)])       
