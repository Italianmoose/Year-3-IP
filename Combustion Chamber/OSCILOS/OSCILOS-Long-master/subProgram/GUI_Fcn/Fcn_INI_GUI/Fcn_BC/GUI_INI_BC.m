function varargout = GUI_INI_BC(varargin)
% GUI_INI_BC MATLAB code for GUI_INI_BC.fig
%      GUI_INI_BC, by itself, creates a new GUI_INI_BC or raises the existing
%      singleton*.
%
%      H = GUI_INI_BC returns the handle to a new GUI_INI_BC or the handle to
%      the existing singleton*.
%
%      GUI_INI_BC('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_INI_BC.M with the given input arguments.
%
%      GUI_INI_BC('Property','Value',...) creates a new GUI_INI_BC or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_INI_BC_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_INI_BC_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_INI_BC

% Last Modified by GUIDE v2.5 19-Jan-2015 11:28:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_INI_BC_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_INI_BC_OutputFcn, ...
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
              

% --- Executes just before GUI_INI_BC is made visible.
function GUI_INI_BC_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_INI_BC (see VARARGIN)
indexEdit = 0;
switch indexEdit 
    case 0
        %--------------------------------------------------------------------------
        dontOpen = false;
        mainGuiInput = find(strcmp(varargin, 'OSCILOS_long'));
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
            handles.indexApp = 0;
            % Update handles structure
            guidata(hObject, handles);
            % Initialization
            GUI_INI_BC_Initialization(hObject, eventdata, handles)
        end
        guidata(hObject, handles);
        handles.output = hObject;
        guidata(hObject, handles);
        if dontOpen
           disp('-----------------------------------------------------');
           disp('This is a subprogram. It cannot be run independently.') 
           disp('Please load the program "OSCILOS_long'' from the ')
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
        assignin('base','CI',CI);                   % save the current information to the works
        handles.indexApp = 0;
        guidata(hObject, handles);  
        GUI_INI_BC_Initialization(hObject, eventdata, handles)
        handles.output = hObject;
        guidata(hObject, handles);
end

function GUI_INI_BC_Initialization(varargin)
hObject = varargin{1};
handles = guidata(hObject);
global CI        
% --------------------------------------
try
    indexExist_CI_BC = any(strcmp('BC',fieldnames(CI)));
    if indexExist_CI_BC == 0;
        CI.BC.StyleInlet    = 2;
        CI.BC.StyleOutlet   = 1;
        CI.BC.num1          = 1;
        CI.BC.den1          = 1;
        CI.BC.tau_d1        = 0;
        CI.BC.num2          = -1;
        CI.BC.den2          = 1;
        CI.BC.tau_d2        = 0;
        CI.BC.A1            = 1;
        CI.BC.Phi1          = 0;
        CI.BC.A2            = 1;
        CI.BC.Phi2          = 0;
        CI.BC.ET.pop_type_model             = 1;
        CI.BC.ET.Dispersion.Delta_tauCs     = 0;
        CI.BC.ET.Dissipation.k              = 0; 
    end
catch
end

assignin('base','CI',CI);                   % save the current information to the workspace
%--------------------------------------
% positions reconfiguration
set(0, 'units', 'points');
screenSize  = get(0, 'ScreenSize');                             % get the screen size
sW          = handles.sW;                                       % screen width
sH          = handles.sH;                                       % screen height
FigW        = sW.*4/5;                                          % window width
FigH        = sH.*5/9;                                          % window height
set(handles.figure,     'units', 'points',...
                        'position',[(screenSize(3)-FigW)./2 (screenSize(4)-FigH)./2 FigW FigH],...
                        'name','Boundary condition configurations',...
                        'color',handles.bgcolor{3});
%----------------------------------------
% pannel axes
set(handles.uipanel_axes,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[FigW*10.5/20 FigH*2.0/20 FigW*9/20 FigH*17.75/20],...
                        'Title','Transfer function of reflection coefficient',...
                        'visible','on',...
                        'highlightcolor',handles.bgcolor{3},...
                        'borderwidth',1,...
                        'fontsize',handles.FontSize(2),...
                        'backgroundcolor',handles.bgcolor{3});  
pannelsize=get(handles.uipanel_axes,'position');
pW=pannelsize(3);
pH=pannelsize(4);                
set(handles.axes_Gain,  'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*2/10 pH*5.0/10 pW*7/10 pH*3.5/10],...
                        'fontsize',handles.FontSize(1),...
                        'color',handles.bgcolor{1},...
                        'box','on');  
set(handles.axes_Phase, 'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*2/10 pH*1.5/10 pW*7/10 pH*3.5/10],...
                        'fontsize',handles.FontSize(1),...
                        'color',handles.bgcolor{1},...
                        'box','on');  
set(handles.pop_Plot,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*0.5/10 pH*8.5/10 pW*3/10 pH*1.0/10],...
                        'fontsize',handles.FontSize(2),...
                        'string',{'Inlet';'Outlet'},...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','left',...
                        'enable','on'); 
                   
guidata(hObject, handles);
%----------------------------------------
% pannels Inlet
set(handles.uipanel_Inlet,   'units', 'points',...
                        'Fontunits','points',...
                        'position',[FigW*0.5/20 FigH*11.0/20 FigW*9.5/20 FigH*8.75/20],...
                        'Title','Inlet configuration',...
                        'visible','on',...
                        'highlightcolor',handles.bgcolor{3},...
                        'borderwidth',1,...
                        'fontsize',handles.FontSize(2),...
                        'backgroundcolor',handles.bgcolor{3}); 
pannelsize=get(handles.uipanel_Inlet,'position');
pW=pannelsize(3);
pH=pannelsize(4);   
set(handles.text_type_Inlet,      'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*0.5/10 pH*7.0/10 pW*5/10 pH*1.2/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Type:',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'horizontalalignment','left',...
                        'visible','off');                         
set(handles.pop_type_Inlet,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*0.5/10 pH*7.5/10 pW*8.5/10 pH*1.2/10],...
                        'fontsize',handles.FontSize(2),...
                        'string',{  'Open end';...
                                    'Closed end';...
                                    'Choked';...
                                    'User defined (Amp. and T.D.)...';...
                                    'User defined (Amp. and Phase)...';...
                                    'Transfer function model from numerators and denominators';...
                                    'Transfer function model loaded from external mat file'},...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','left',...
                        'enable','on',...
                        'value',CI.BC.StyleInlet);  
set(handles.text_Amp_Inlet,      'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*0.5/10 pH*4.75/10 pW*5/10 pH*1.2/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Amplitude [-]:',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'horizontalalignment','left',...
                        'visible','on');                   
set(handles.edit_Amp_Inlet,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*6/10 pH*5/10 pW*3/10 pH*1.2/10],...
                        'fontsize',handles.FontSize(2),...
                        'string',1,...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','right',...
                        'visible','on');  
set(handles.text_td_Inlet,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*0.5/10 pH*2.75/10 pW*5/10 pH*1.2/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Time delay [ms]:',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'horizontalalignment','left',...
                        'visible','on');                   
set(handles.edit_td_Inlet,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*6/10 pH*3/10 pW*3/10 pH*1.2/10],...
                        'fontsize',handles.FontSize(2),...
                        'string',0,...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','right',...
                        'visible','on');  
set(handles.text_add_Inlet,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*0.5/10 pH*0.75/10 pW*5/10 pH*1.2/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Time delay [ms]:',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'horizontalalignment','left',...
                        'visible','on');                   
set(handles.edit_add_Inlet,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*6/10 pH*1.0/10 pW*3/10 pH*1.2/10],...
                        'fontsize',handles.FontSize(2),...
                        'string',0,...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','right',...
                        'visible','on');  
set(handles.pb_inlet,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*0.5/10 pH*0.75/10 pW*8.5/10 pH*1.5/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Load the transfer function model and time delay',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'horizontalalignment','right',...
                        'Enable','on'); 
%----------------------------------------
% pannels Outlet
set(handles.uipanel_Outlet,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[FigW*0.5/20 FigH*2.0/20 FigW*9.5/20 FigH*8.75/20],...
                        'Title','Outlet configuration',...
                        'visible','on',...
                        'highlightcolor',handles.bgcolor{3},...
                        'borderwidth',1,...
                        'fontsize',handles.FontSize(2),...
                        'backgroundcolor',handles.bgcolor{3});  
pannelsize=get(handles.uipanel_Outlet,'position');
pW=pannelsize(3);
pH=pannelsize(4);   
set(handles.text_type_Outlet,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*0.5/10 pH*6.5/10 pW*5/10 pH*1.2/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Type:',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'horizontalalignment','left',...
                        'visible','off');  
set(handles.pop_type_Outlet,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*0.5/10 pH*7.5/10 pW*8.5/10 pH*1.2/10],...
                        'fontsize',handles.FontSize(2),...
                        'string',{  'Open end';...
                                    'Closed end';...
                                    'Choked';...
                                    'User defined (Amp. and T.D.)...';...
                                    'User defined (Amp. and Phase)...';...
                                    'Transfer function model from numerators and denominators';...
                                    'Transfer function model loaded from external mat file'},...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','left',...
                        'Enable','on',...
                        'value',CI.BC.StyleOutlet);  
set(handles.text_Amp_Outlet,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*0.5/10 pH*4.75/10 pW*5/10 pH*1.2/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Amplitude [-]:',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'horizontalalignment','left',...
                        'visible','on');                   
set(handles.edit_Amp_Outlet,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*6/10 pH*5/10 pW*3/10 pH*1.2/10],...
                        'fontsize',handles.FontSize(2),...
                        'string',1,...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','right',...
                        'visible','on');  
set(handles.text_td_Outlet,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*0.5/10 pH*2.75/10 pW*5/10 pH*1.2/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Time delay [ms]:',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'horizontalalignment','left',....
                        'visible','on');                   
set(handles.edit_td_Outlet,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*6/10 pH*3/10 pW*3/10 pH*1.2/10],...
                        'fontsize',handles.FontSize(2),...
                        'string',0,...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','right',....
                        'visible','on');  
set(handles.text_add_Outlet,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*0.5/10 pH*0.75/10 pW*5/10 pH*1.2/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Time delay [ms]:',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'horizontalalignment','left',...
                        'visible','on');                   
set(handles.edit_add_Outlet,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*6/10 pH*1.0/10 pW*3/10 pH*1.2/10],...
                        'fontsize',handles.FontSize(2),...
                        'string',0,...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','right',...
                        'visible','on'); 
set(handles.pb_outlet,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*0.5/10 pH*0.75/10 pW*8.5/10 pH*1.5/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Load the transfer function model and time delay',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'horizontalalignment','right',...
                        'Enable','on'); 
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
set(handles.pb_Apply,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*1/10 pH*1/10 pW*2.0/10 pH*7/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Plot',...
                        'backgroundcolor',handles.bgcolor{3});
set(handles.pb_SaveFig,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*4/10 pH*1/10 pW*2.0/10 pH*7/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Save figure',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'enable','off',...
                        'visible','off'); 
set(handles.pb_OK,       'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*4/10 pH*1/10 pW*2.0/10 pH*7/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','OK',...
                        'backgroundcolor',handles.bgcolor{3});
set(handles.pb_Cancel,....
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*7/10 pH*1/10 pW*2.0/10 pH*7/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Cancel',...
                        'backgroundcolor',handles.bgcolor{3});
%
assignin('base','CI',CI);                   % save the current information to the workspace
guidata(hObject, handles);
Fcn_GUI_INI_BC_POP1_update(hObject) 
handles=guidata(hObject);
guidata(hObject, handles);
Fcn_GUI_INI_BC_POP2_update(hObject) 
handles=guidata(hObject);
guidata(hObject, handles);
Fcn_GUI_INI_BC_PLOT_BC_TF(hObject)
handles=guidata(hObject);
guidata(hObject, handles);

% 

function Fcn_GUI_INI_BC_POP1_update(varargin)
hObject = varargin{1};
handles = guidata(hObject);
global CI
pop_type_Inlet      = get(handles.pop_type_Inlet,'Value');
%
switch pop_type_Inlet
    case 4
            set(handles.edit_Amp_Inlet,             'visible','on',...
                                                    'string', num2str(CI.BC.A1));
            set(handles.text_Amp_Inlet,             'visible','on',...
                                                    'string','Amplitude [-]'); 
            set(handles.edit_td_Inlet,              'visible','on',...
                                                    'string', num2str(CI.BC.tau_d1*1e3));
            set(handles.text_td_Inlet,              'visible','on',...
                                                    'string','Time delay [ms]:');
            set(handles.edit_add_Inlet,             'visible','off');
            set(handles.text_add_Inlet,             'visible','off');
            set(handles.pb_inlet,                   'visible','off');
    case 5
            set(handles.edit_Amp_Inlet,             'visible','on',...
                                                    'string', num2str(CI.BC.A1));
            set(handles.text_Amp_Inlet,             'visible','on',...
                                                    'string','Amplitude [-]'); 
            set(handles.edit_td_Inlet,              'visible','on',...
                                                    'string', num2str(CI.BC.Phi1/pi*180));
            set(handles.text_td_Inlet,              'visible','on',...
                                                    'string','Phase [deg]:');
            set(handles.edit_add_Inlet,             'visible','off');
            set(handles.text_add_Inlet,             'visible','off');
            set(handles.pb_inlet,                   'visible','off');
    case 6
        handles.indexApp = Fcn_GUI_INI_BC_helpdlg(handles.indexApp);
        strNum = ['[' num2str(CI.BC.num1) ']'];
        strDen = ['[' num2str(CI.BC.den1) ']'];
            set(handles.edit_Amp_Inlet,             'visible','on',...
                                                    'string', strNum);
            set(handles.text_Amp_Inlet,             'visible','on',...
                                                    'string', 'Numerator(s):');
            set(handles.edit_td_Inlet,              'visible','on',...
                                                    'string', strDen);
            set(handles.text_td_Inlet,              'visible','on',...
                                                    'string', 'Denominator(s):');
            set(handles.edit_add_Inlet,             'visible','on',...
                                                    'string', num2str(CI.BC.tau_d1.*1e3));
            set(handles.text_add_Inlet,             'visible','on',...
                                                    'string','Time delay [ms]:');
            set(handles.pb_inlet,                   'visible','off');
    case 7
            set(handles.edit_Amp_Inlet,             'visible','off');
            set(handles.text_Amp_Inlet,             'visible','off'); 
            set(handles.edit_td_Inlet,              'visible','on',...
                                                    'string', num2str(CI.BC.tau_d1.*1e3));
            set(handles.text_td_Inlet,              'visible','on',...
                                                    'string','Time delay [ms]:');
            set(handles.edit_add_Inlet,             'visible','off');
            set(handles.text_add_Inlet,             'visible','off');
            set(handles.pb_inlet,                   'visible','on');
            
    otherwise
            set(handles.edit_Amp_Inlet,             'visible','off');
            set(handles.text_Amp_Inlet,             'visible','off'); 
            set(handles.edit_td_Inlet,              'visible','off');
            set(handles.text_td_Inlet,              'visible','off');   
            set(handles.edit_add_Inlet,             'visible','off');
            set(handles.text_add_Inlet,             'visible','off');
            set(handles.pb_inlet,                   'visible','off');
end
guidata(hObject, handles);
% ----------------------
function Fcn_GUI_INI_BC_POP2_update(varargin)
hObject = varargin{1};
handles = guidata(hObject);
global CI
pop_type_Outlet     = get(handles.pop_type_Outlet,'Value');
switch pop_type_Outlet
    case 4
            set(handles.edit_Amp_Outlet,            'visible','on',...
                                                    'string', num2str(CI.BC.A2));
            set(handles.text_Amp_Outlet,            'visible','on',...
                                                    'string','Amplitude [-]'); 
            set(handles.edit_td_Outlet,             'visible','on',...
                                                    'string', num2str(CI.BC.tau_d2*1e3));
            set(handles.text_td_Outlet,             'visible','on',...
                                                    'string','Time delay [ms]:');
            set(handles.edit_add_Outlet,            'visible','off');
            set(handles.text_add_Outlet,            'visible','off');
            set(handles.pb_outlet,                   'visible','off');
    case 5
            set(handles.edit_Amp_Outlet,            'visible','on',...
                                                    'string', num2str(CI.BC.A2));
            set(handles.text_Amp_Outlet,            'visible','on',...
                                                    'string','Amplitude [-]'); 
            set(handles.edit_td_Outlet,             'visible','on',...
                                                    'string', num2str(CI.BC.Phi2/pi*180));
            set(handles.text_td_Outlet,             'visible','on',...
                                                    'string','Phase [deg]:');
            set(handles.edit_add_Outlet,            'visible','off');
            set(handles.text_add_Outlet,            'visible','off');
            set(handles.pb_outlet,                   'visible','off');
    case 6
            handles.indexApp = Fcn_GUI_INI_BC_helpdlg(handles.indexApp);
            strNum = ['[' num2str(CI.BC.num2) ']'];
            strDen = ['[' num2str(CI.BC.den2) ']'];
            set(handles.edit_Amp_Outlet,            'visible','on',...
                                                    'string', strNum);
            set(handles.text_Amp_Outlet,            'visible','on',...
                                                    'string', 'Numerator(s):');
            set(handles.edit_td_Outlet,             'visible','on',...
                                                    'string', strDen);
            set(handles.text_td_Outlet,             'visible','on',...
                                                    'string', 'Denominator(s):');
            set(handles.edit_add_Outlet,            'visible','on',...
                                                    'string', num2str(CI.BC.tau_d2.*1e3));
            set(handles.text_add_Outlet,            'visible','on',...
                                                    'string','Time delay [ms]:');
            set(handles.pb_outlet,                   'visible','off');
    case 3
            set(handles.edit_Amp_Outlet,             'visible','off');
            set(handles.text_Amp_Outlet,             'visible','off'); 
            set(handles.edit_td_Outlet,              'visible','off');
            set(handles.text_td_Outlet,              'visible','off');   
            set(handles.edit_add_Outlet,             'visible','off');
            set(handles.text_add_Outlet,             'visible','off');
            
            % Construct a questdlg with three options
            choice = questdlg('Include indirect noise from entropy waves?', ...
                'Indirect noise','Yes','No','No');
            % Handle response
            switch choice
                case 'Yes'
                    GUI_INI_BC_Entropy('GUI_INI_BC', handles.figure);
                case 'No'
                    CI.BC.ET.pop_type_model             = 1;
                    CI.BC.ET.Dispersion.Delta_tauCs     = 0;
                    CI.BC.ET.Dissipation.k              = 0; 
                    assignin('base','CI',CI);                   % save the current information to the workspace
            end
            set(handles.pb_outlet,                   'visible','off');
    case 7
            set(handles.edit_Amp_Outlet,             'visible','off');
            set(handles.text_Amp_Outlet,             'visible','off'); 
            set(handles.edit_td_Outlet,              'visible','on',...
                                                    'string', num2str(CI.BC.tau_d2.*1e3));
            set(handles.text_td_Outlet,              'visible','on',...
                                                    'string','Time delay [ms]:');
            set(handles.edit_add_Outlet,             'visible','off');
            set(handles.text_add_Outlet,             'visible','off');
            set(handles.pb_outlet,                   'visible','on');
    otherwise
            set(handles.edit_Amp_Outlet,             'visible','off');
            set(handles.text_Amp_Outlet,             'visible','off'); 
            set(handles.edit_td_Outlet,              'visible','off');
            set(handles.text_td_Outlet,              'visible','off');   
            set(handles.edit_add_Outlet,             'visible','off');
            set(handles.text_add_Outlet,             'visible','off');
            set(handles.pb_outlet,                   'visible','off');
end
guidata(hObject, handles);

function indexApp = Fcn_GUI_INI_BC_helpdlg(indexApp)
if indexApp < 2
string={...
'The numerator coefficient can be a vector or matrix expression.';...
'The denominator coefficient must be a vector.';...
'The output width equals the number of rows in the numerator coefficient.';
'You should specify the coefficients in descending order of powers of s.'};
helpdlg(string,'Set transfer function coefficients')
end
indexApp = indexApp + 1;

function Fcn_GUI_INI_BC_Value_update(varargin)
hObject = varargin{1};
handles = guidata(hObject);
global CI
CI.BC.StyleInlet = get(handles.pop_type_Inlet,'Value');
%
switch CI.BC.StyleInlet
    case 1      % open end, R=-1;
        CI.BC.num1      = -1;
        CI.BC.den1      = 1;
        CI.BC.tau_d1    = 0;
    case 2      % closed end, R=1; 
        CI.BC.num1      = 1;
        CI.BC.den1      = 1;
        CI.BC.tau_d1    = 0;
    case 3      % choked
        gamma           = CI.TP.gamma(1,1);
        M1              = CI.TP.M_mean(1,1);
        TEMP            = gamma*M1/(1+(gamma-1)*M1^2);
        %
        %%% test of boundary conditions in Coh See's paper, CST 2013, added
        %%% in 2015-02-11
%         TEMP = M1;
        %%% END
        CI.BC.num1      = (1-TEMP)/(1+TEMP);
        CI.BC.den1      = 1;
        CI.BC.tau_d1    = 0;
        %
    case 4
        CI.BC.num1      = str2double(get(handles.edit_Amp_Inlet,'String'));
        CI.BC.den1      = 1;
        CI.BC.tau_d1    = str2double(get(handles.edit_td_Inlet, 'String'))./1000; 
        CI.BC.A1        = CI.BC.num1;
    case 5
        A               = str2double(get(handles.edit_Amp_Inlet,'String'));
        Phi             = str2double(get(handles.edit_td_Inlet,'String'))./180*pi;
        CI.BC.num1      = A*exp(-Phi*1i);
        CI.BC.den1      = 1;
        CI.BC.tau_d1    = 0;
        CI.BC.A1        = A;
        CI.BC.Phi1      = Phi;
    case 6
        CI.BC.num1      = str2num(get(handles.edit_Amp_Inlet,'String'));
        CI.BC.den1      = str2num(get(handles.edit_td_Inlet,'String'));
        CI.BC.tau_d1    = str2num(get(handles.edit_add_Inlet, 'String'))./1000; 
    case 7
        [CI.BC.num1,CI.BC.den1] = tfdata(handles.sysFit1,'v');
        CI.BC.tau_d1    = str2num(get(handles.edit_td_Inlet, 'String'))./1000; 
end
assignin('base','CI',CI);                   % save the current information to the workspace
%
CI.BC.StyleOutlet     = get(handles.pop_type_Outlet,'Value');
%
switch CI.BC.StyleOutlet
    case 1      % open end, R=-1;
        CI.BC.num2 = -1;
        CI.BC.den2 = 1;
        CI.BC.tau_d2 = 0;
    case 2      % closed end, R=1; 
        CI.BC.num2 = 1;
        CI.BC.den2 = 1;
        CI.BC.tau_d2 = 0;
    case 3     % choked
        gamma   = CI.TP.gamma(1,end);
        M2      = CI.TP.M_mean(1,end);
        TEMP    = (gamma-1)*M2/2;
        %
        CI.BC.num2  = (1-TEMP)/(1+TEMP);
        CI.BC.den2 = 1;
        CI.BC.tau_d2 = 0;
        %
    case 4
        CI.BC.num2      = str2double(get(handles.edit_Amp_Outlet,'String'));
        CI.BC.den2      = 1;
        CI.BC.tau_d2    = str2double(get(handles.edit_td_Outlet, 'String'))./1000;
        CI.BC.A2        = CI.BC.num2;
    case 5
        A               = str2double(get(handles.edit_Amp_Outlet,'String'));
        Phi             = str2double(get(handles.edit_td_Outlet,'String'))./180*pi;
        CI.BC.num2      = A*exp(-Phi*1i);
        CI.BC.den2      = 1;
        CI.BC.tau_d2    = 0;
        CI.BC.A2        = A;
        CI.BC.Phi2      = Phi;
    case 6
        CI.BC.num2          = str2num(get(handles.edit_Amp_Outlet,'String'));
        CI.BC.den2          = str2num(get(handles.edit_td_Outlet,'String'));
        CI.BC.tau_d2        = str2num(get(handles.edit_add_Outlet, 'String'))./1000;
    case 7
        [CI.BC.num2,CI.BC.den2] = tfdata(handles.sysFit2,'v');
        CI.BC.tau_d2    = str2num(get(handles.edit_td_Outlet, 'String'))./1000; 
end
switch CI.BC.StyleOutlet
    case 3
    otherwise
        CI.BC.ET.pop_type_model         = 1;
        CI.BC.ET.Dissipation.k          = 0;
        CI.BC.ET.Dispersion.Delta_tauCs = 0e-3;
end
assignin('base','CI',CI);                   % save the current information to the workspace        
guidata(hObject, handles);


function pop_type_Inlet_Callback(hObject, eventdata, handles)
% hObject    handle to pop_type_Inlet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.pb_SaveFig,'enable','off');
Fcn_GUI_INI_BC_POP1_update(hObject)
handles = guidata(hObject);
guidata(hObject, handles);


%--------------------------------------------------------------------------
% Outlet pannel
% --- Executes on selection change in pop_type_Outlet.
function pop_type_Outlet_Callback(hObject, eventdata, handles)
% hObject    handle to pop_type_Outlet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.pb_SaveFig,'enable','off');
Fcn_GUI_INI_BC_POP2_update(hObject)
handles = guidata(hObject);
guidata(hObject, handles);


 

% --- Outputs from this function are returned to the command line.
function varargout = GUI_INI_BC_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
try
varargout{1} = handles.output;
end

%
% --- Executes on selection change in pop_Plot.
function pop_Plot_Callback(hObject, eventdata, handles)
% hObject    handle to pop_Plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pop_Plot contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_Plot
set(handles.pb_SaveFig,'enable','off');
%
% --- Executes during object creation, after setting all properties.
function pop_Plot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_Plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%
% --- plot the transfer function of the boundary condition
% function Fcn_GUI_INI_BC_PLOT_BC_TF(hAxes1, hAxes2, hObject, handles)
function Fcn_GUI_INI_BC_PLOT_BC_TF(varargin)
hObject     = varargin{1};
Fcn_GUI_INI_BC_Value_update(hObject)
handles     = guidata(hObject);
global CI
hAxes1      = handles.axes_Gain;
hAxes2      = handles.axes_Phase;
fontSize1   = handles.FontSize(1);
fontSize2   = handles.FontSize(2);
f_sample    = 0:1000;                                       % sampling frequency
s           = 2*pi*f_sample*1i;

R1 = polyval(CI.BC.num1,s)./polyval(CI.BC.den1,s).*exp(-s*CI.BC.tau_d1);
R2 = polyval(CI.BC.num2,s)./polyval(CI.BC.den2,s).*exp(-s*CI.BC.tau_d2);

pop_type_Inlet      = get(handles.pop_type_Inlet,'Value');
pop_type_Outlet     = get(handles.pop_type_Outlet,'Value');
pop_PLot            = get(handles.pop_Plot,'Value');
%
guidata(hObject, handles)
%
switch pop_PLot
    case 1
        R = R1;
        colorPlot={'b','g'};
    case 2
        R = R2;
        colorPlot={'r','g'};
end
TextPLot={'Inlet','Outlet'};
%--------------------------------
cla(hAxes1)
axes(hAxes1)
hold on
    color_type = Fcn_color_set(5);
    hLine1=plot(hAxes1,f_sample,abs(R),'-','color',colorPlot{1},'Linewidth',2); 
    set(hAxes1,'YColor','k','Box','on','ygrid','on','xgrid','on');
    set(hAxes1,'FontName','Helvetica','FontSize',fontSize1,'LineWidth',1)
    xlabel(hAxes1,'','Color','k','Interpreter','LaTex','FontSize',fontSize1);
    ylabel(hAxes1,'Gain [-]','Color','k','Interpreter','LaTex','FontSize',fontSize1)
    set(hAxes1,'xlim',[0 1000],'xTick',200:200:1000,'xticklabel',{},...
    'YAxisLocation','left','Color','w');
    set(hAxes1,'ylim',[0 1.2],'yTick',0:0.2:1.2,'yticklabel',{'',0.2:0.2:1})
    xLim=get(hAxes1,'xlim');
    yLim=get(hAxes1,'ylim');
    text(xLim(1)+0.5*xLim(2),yLim(1)+1.1*yLim(2),TextPLot{pop_PLot},'interpreter','latex',...
        'fontsize',fontSize1,'parent',hAxes1,'horizontalAlignment','center')
hold off
%--------------------------------
cla(hAxes2)
axes(hAxes2)
hold on
    hLine2=plot(hAxes2,f_sample,unwrap(angle(R),1.9*pi)./pi,'-','color', colorPlot{1},'Linewidth',2);
    set(hAxes2,'YColor','k','Box','on','ygrid','on','xgrid','on');
    set(hAxes2,'FontName','Helvetica','FontSize',fontSize1,'LineWidth',1)
    xlabel(hAxes2,'$f$ [Hz]','Color','k','Interpreter','LaTex','FontSize',fontSize1);
    ylabel(hAxes2,'Phase/$\pi$ [-]','Color','k','Interpreter','LaTex','FontSize',fontSize1)
    set(hAxes2,'xlim',get(hAxes1,'xlim'),'xTick',get(hAxes1,'xTick'),...
    'YAxisLocation','left','Color','w');
hold off
assignin('base','CI',CI);                   % save the current information to the workspace
%--------------------------------------------------------------------------


% --- Update the data when clicking 'OK' or 'Apply'
function Fcn_GUI_INI_BC_Update_Data(hObject, eventdata, handles)
global CI
Fcn_GUI_INI_BC_PLOT_BC_TF(hObject, eventdata, handles) 
% try
main = handles.MainGUI;
% Obtain handles using GUIDATA with the caller's handle 
if(ishandle(main))
    mainHandles = guidata(main);
    % -------------------------------------
    changeMain1 = mainHandles.FREQ;      % set the Frequency pop-up menu enable
    changeMain2 = mainHandles.TD; % set the TD pop-up menu enable
    % default 
    set(changeMain1, 'Enable', 'on');
    set(changeMain2, 'Enable', 'on');
    % ------------
    if ~isempty(CI.CD.indexHP)            % with heat perturbations
        if ~isempty(find(CI.FM.indexFM == 4))   % Frequency calculation can only be run if the G-EQuation is not being used
            set(changeMain1, 'Enable', 'off');
            set(changeMain2, 'Enable', 'on');
        end
        if ~isempty(find(CI.FM.indexFM == 3))   % with loaded experimental/CFD FDF, time domain code is not avaliable
            set(changeMain1, 'Enable', 'on');
            set(changeMain2, 'Enable', 'off');
        end
    end
    % ---------------
    String_Listbox=get(mainHandles.listbox_Info,'string');
    ind=find(ismember(String_Listbox,'<HTML><FONT color="blue">Information 4:'));
    nLength=size(String_Listbox);
    if isempty(ind)
        indStart=nLength(1);
    else
        indStart=ind-1;
        for i=nLength(1):-1:indStart 
            String_Listbox(i)=[];
        end
    end
    String_Listbox{indStart+1}=['<HTML><FONT color="blue">Information 4:'];
    String_Listbox{indStart+2}=['<HTML><FONT color="blue">Boundary condition:'];
    String_Listbox{indStart+3}=['Inlet boundary condition type:' num2str(CI.BC.StyleInlet) ];
    String_Listbox{indStart+4}=['Outlet boundary condition type:' num2str(CI.BC.StyleOutlet) ];
    set(mainHandles.listbox_Info,'string',String_Listbox);
end
% end
guidata(hObject, handles);
% guidata(hObject, handles);
assignin('base','CI',CI);                   % save the current information to the workspace
%
% --- Executes on button press in pb_OK.
function pb_OK_Callback(hObject, eventdata, handles)
% hObject    handle to pb_OK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Fcn_GUI_INI_BC_Update_Data(hObject, eventdata, handles)
global CI
CI.IsRun.GUI_INI_BC = 1;
assignin('base','CI',CI); 
delete(handles.figure);

% --- Executes on button press in pb_Apply.
function pb_Apply_Callback(hObject, eventdata, handles)
% hObject    handle to pb_Apply (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Fcn_GUI_INI_BC_Update_Data(hObject, eventdata, handles)
set(handles.pb_SaveFig,'enable','on');

% --- Executes on button press in pb_Cancel.
function pb_Cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pb_Cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(handles.figure);


% --- Executes on button press in pb_SaveFig.
function pb_SaveFig_Callback(hObject, eventdata, handles)
global CI
% hObject    handle to pb_SaveFig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);
% pop_PLot = get(handles.pop_Plot,'Value');
% pop_type_Inlet=get(handles.pop_type_Inlet,'Value');
% pop_type_Outlet=get(handles.pop_type_Outlet,'Value');
Fig = figure;
copyobj(handles.axes_Gain, Fig);
copyobj(handles.axes_Phase, Fig);
set(Fig,        'units','points')
posFig = get(handles.figure,'position');
hAxes = get(Fig,'children');
set(hAxes(1),       'units','points',...
                    'position',[60 60 200 150],...
                    'ActivePositionProperty','position')
set(hAxes(2),       'units','points',...
                    'position',[60 210 200 150],...
                    'ActivePositionProperty','position')
posAxesOuter = [0 0 300 400];
set(Fig,        'units','points',...
                'position', [posFig(1)+0.5*posFig(3)-0.5*posAxesOuter(3),...
                            posFig(2)+0.5*posFig(4)-0.5*posAxesOuter(4),...
                            posAxesOuter(3:4)]) 




% --- Executes during object creation, after setting all properties.
function pop_type_Inlet_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_type_Inlet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Amp_Inlet_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Amp_Inlet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% datEdit = str2double(get(hObject, 'String'));
% if isnan(datEdit)
%     set(hObject, 'String', 0);
%     errordlg('Input must be a number','Error');
% end


% --- Executes during object creation, after setting all properties.
function edit_Amp_Inlet_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Amp_Inlet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_td_Inlet_Callback(hObject, eventdata, handles)
% hObject    handle to edit_td_Inlet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% datEdit = str2double(get(hObject, 'String'));
% if isnan(datEdit)
%     set(hObject, 'String', 0);
%     errordlg('Input must be a number','Error');
% end


% --- Executes during object creation, after setting all properties.
function edit_td_Inlet_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_td_Inlet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_add_Inlet_Callback(hObject, eventdata, handles)
% hObject    handle to edit_add_Inlet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% datEdit = str2double(get(hObject, 'String'));
% if isnan(datEdit)
%     set(hObject, 'String', 0);
%     errordlg('Input must be a number','Error');
% end

% --- Executes during object creation, after setting all properties.
function edit_add_Inlet_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_add_Inlet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function pop_type_Outlet_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_type_Outlet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%
%
function edit_Amp_Outlet_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Amp_Outlet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% datEdit = str2double(get(hObject, 'String'));
% if isnan(datEdit)
%     set(hObject, 'String', 0);
%     errordlg('Input must be a number','Error');
% end

% --- Executes during object creation, after setting all properties.
function edit_Amp_Outlet_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Amp_Outlet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%
%
function edit_td_Outlet_Callback(hObject, eventdata, handles)
% hObject    handle to edit_td_Outlet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% datEdit = str2double(get(hObject, 'String'));
% if isnan(datEdit)
%     set(hObject, 'String', 0);
%     errordlg('Input must be a number','Error');
% end

% --- Executes during object creation, after setting all properties.
function edit_td_Outlet_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_td_Outlet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_add_Outlet_Callback(hObject, eventdata, handles)
% hObject    handle to edit_add_Outlet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% datEdit = str2double(get(hObject, 'String'));
% if isnan(datEdit)
%     set(hObject, 'String', 0);
%     errordlg('Input must be a number','Error');
% end

% --- Executes during object creation, after setting all properties.
function edit_add_Outlet_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_add_Outlet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%--------------------------------------------------------------------------
% --- Executes when user attempts to close figure.
function figure_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
delete(hObject);
%--------------------------------------------------------------------------
%-- Color set
function color_type=Fcn_color_set(color_number)
cmap = colormap(hot);
n_sample=round(linspace(1,64,color_number));
for ss=1:length(n_sample)
    color_type(ss,1:3)=cmap(n_sample(ss),1:3);
end


% --------------------------------------------------------------------
function uipushtool1_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% hObject    handle to pb_SaveFig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);
% pop_PLot = get(handles.pop_Plot,'Value');
% pop_type_Inlet=get(handles.pop_type_Inlet,'Value');
% pop_type_Outlet=get(handles.pop_type_Outlet,'Value');
Fig = figure;
copyobj(handles.axes_Gain, Fig);
copyobj(handles.axes_Phase, Fig);
set(Fig,        'units','points')
posFig = get(handles.figure,'position');
hAxes = get(Fig,'children');
set(hAxes(1),       'units','points',...
                    'position',[60 60 200 150],...
                    'ActivePositionProperty','position')
set(hAxes(2),       'units','points',...
                    'position',[60 210 200 150],...
                    'ActivePositionProperty','position')
posAxesOuter = [0 0 300 400];
set(Fig,        'units','points',...
                'position', [posFig(1)+0.5*posFig(3)-0.5*posAxesOuter(3),...
                            posFig(2)+0.5*posFig(4)-0.5*posAxesOuter(4),...
                            posAxesOuter(3:4)]) 


% --- Executes on button press in pb_outlet.
function pb_outlet_Callback(hObject, eventdata, handles)
% hObject    handle to pb_outlet (see GCBO)
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
    handles.sysFit2 = sys;
    set(handles.edit_td_Outlet, 'string', num2str(td.*1e3));
    guidata(hObject, handles);
end

% --- Executes on button press in pb_inlet.
function pb_inlet_Callback(hObject, eventdata, handles)
% hObject    handle to pb_inlet (see GCBO)
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
    handles.sysFit1 = sys;
    set(handles.edit_td_Inlet, 'string', num2str(td.*1e3));
    guidata(hObject, handles);
end