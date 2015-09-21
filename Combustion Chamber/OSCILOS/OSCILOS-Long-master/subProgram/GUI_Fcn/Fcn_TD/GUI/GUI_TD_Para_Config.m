function varargout = GUI_TD_Para_Config(varargin)
% GUI_TD_PARA_CONFIG MATLAB code for GUI_TD_Para_Config.fig
%      GUI_TD_PARA_CONFIG, by itself, creates a new GUI_TD_PARA_CONFIG or raises the existing
%      singleton*.
%
%      H = GUI_TD_PARA_CONFIG returns the handle to a new GUI_TD_PARA_CONFIG or the handle to
%      the existing singleton*.
%
%      GUI_TD_PARA_CONFIG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_TD_PARA_CONFIG.M with the given input arguments.
%
%      GUI_TD_PARA_CONFIG('Property','Value',...) creates a new GUI_TD_PARA_CONFIG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_TD_Para_Config_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_TD_Para_Config_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_TD_Para_Config

% Last Modified by GUIDE v2.5 22-Jan-2015 19:21:13

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_TD_Para_Config_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_TD_Para_Config_OutputFcn, ...
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


% --- Executes just before GUI_TD_Para_Config is made visible.
function GUI_TD_Para_Config_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_TD_Para_Config (see VARARGIN)
% --------------------------------------------------------------------------
handles.indexEdit = 0;
switch handles.indexEdit 
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
            % Update handles structure
            guidata(hObject, handles);
            % Initialization
            GUI_TD_Para_Config_Initialization(hObject, eventdata, handles)
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
        addpath(genpath('./'))                                                      % add directories to search path
        handles = Fcn_GUI_default_configuration(handles);
        handles.output = hObject;
        guidata(hObject, handles);
        CI.IsRun.GUI_TD_Para_Config = 0;
        assignin('base','CI',CI);                   % save the current information to the works
        guidata(hObject, handles);  
        GUI_TD_Para_Config_Initialization(hObject, eventdata, handles)
        guidata(hObject, handles);
        handles.output = hObject;
        guidata(hObject, handles);
end
%
%
% --- Outputs from this function are returned to the command line.
function varargout = GUI_TD_Para_Config_OutputFcn(hObject, eventdata, handles) 
try
varargout{1} = handles.output;
end
%
%
function GUI_TD_Para_Config_Initialization(varargin)
hObject = varargin{1};
handles = guidata(hObject); 
Fcn_PreProcessing
global CI
set(0, 'units', 'points');
screenSize  = get(0, 'ScreenSize');                 % get the screen size
sW          = handles.sW;                           % screen width
sH          = handles.sH;                          % screen height
FigW=sW.*4/5;                                        % window width
FigH=sH.*3/4;                                        % window height

%
set(handles.figure,     'units', 'points',...
                        'position',[(screenSize(3)-FigW)./2 (screenSize(4)-FigH)./2 FigW FigH],...
                        'name','Parameters configuration',...
                        'color',handles.bgcolor{3});
%----------------------------------------
% pannel time end and tiem step
set(handles.uipanel_Time,...
                        'units', 'points',...
                        'fontUnits','points',...
                        'position',[FigW*0.5/20 FigH*14.5/20 FigW*9/20 FigH*5.25/20],...
                        'Title','Simulation time and step',...
                        'visible','on',...
                        'highlightcolor',handles.bgcolor{3},...
                        'borderwidth',1,...
                        'fontsize',handles.FontSize(2),...
                        'backgroundcolor',handles.bgcolor{3});  
pannelsize = get(handles.uipanel_Time,'position');
pW=pannelsize(3);
pH=pannelsize(4);  
set(handles.text_a1,...
                        'units', 'points',...
                        'fontUnits','points',...
                        'position',[pW*0.5/10 pH*6.0/10 pW*6/10 pH*1.5/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Minimum time delay: [s]',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'horizontalalignment','left');                         
set(handles.edit_a1,...
                        'units', 'points',...
                        'fontUnits','points',...
                        'position',[pW*7/10 pH*6/10 pW*2.5/10 pH*1.75/10],...
                        'fontsize',handles.FontSize(2),...
                        'string',1,...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','right',...
                        'Enable','off');
set(handles.text_a2,...
                        'units', 'points',...
                        'fontUnits','points',...
                        'position',[pW*0.5/10 pH*3.5/10 pW*6/10 pH*1.5/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Stop time: [s]',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'horizontalalignment','left');                         
set(handles.edit_a2,...
                        'units', 'points',...
                        'fontUnits','points',...
                        'position',[pW*7/10 pH*3.5/10 pW*2.5/10 pH*1.75/10],...
                        'fontsize',handles.FontSize(2),...
                        'string',1,...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','right',...
                        'Enable','on');
set(handles.text_a3,...
                        'units', 'points',...
                        'fontUnits','points',...
                        'position',[pW*0.5/10 pH*1/10 pW*6/10 pH*1.5/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Time step: [s]',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'horizontalalignment','left');                         
set(handles.edit_a3,...
                        'units', 'points',...
                        'fontUnits','points',...
                        'position',[pW*7/10 pH*1.0/10 pW*2.5/10 pH*1.75/10],...
                        'fontsize',handles.FontSize(2),...
                        'string',1e-4,...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','right',...
                        'Enable','on');                    
 
%----------------------------------------
% pannel speed
set(handles.uipanel_Speed,...
                        'units', 'points',...
                        'fontUnits','points',...
                        'position',[FigW*0.5/20 FigH*10.25/20 FigW*9/20 FigH*4/20],...
                        'Title','Simulation samples per loop',...
                        'visible','on',...
                        'highlightcolor',handles.bgcolor{3},...
                        'borderwidth',1,...
                        'fontsize',handles.FontSize(2),...
                        'backgroundcolor',handles.bgcolor{3});  
pannelsize = get(handles.uipanel_Speed,'position');
pW=pannelsize(3);
pH=pannelsize(4);  
set(handles.text_d1,...
                        'units', 'points',...
                        'fontUnits','points',...
                        'position',[pW*0.5/10 pH*5.0/10 pW*6/10 pH*1.5/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Maximum number of samples per loop: [-]',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'horizontalalignment','left');    
set(handles.edit_d1,...
                        'units', 'points',...
                        'fontUnits','points',...
                        'position',[pW*7/10 pH*5/10 pW*2.5/10 pH*2.25/10],...
                        'fontsize',handles.FontSize(2),...
                        'string',1,...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','right',...
                        'Enable','off');    
set(handles.text_d2,...
                        'units', 'points',...
                        'fontUnits','points',...
                        'position',[pW*0.5/10 pH*1.5/10 pW*6/10 pH*1.5/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Number of samples per loop: [-]',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'horizontalalignment','left');    
set(handles.edit_d2,...
                        'units', 'points',...
                        'fontUnits','points',...
                        'position',[pW*7/10 pH*1.5/10 pW*2.5/10 pH*2.25/10],...
                        'fontsize',handles.FontSize(2),...
                        'string',1,...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','right',...
                        'Enable','on');        
%----------------------------------------
% pannel uRatio limit
set(handles.uipanel_Limit,...
                        'units', 'points',...
                        'fontUnits','points',...
                        'position',[FigW*0.5/20 FigH*7.5/20 FigW*9/20 FigH*2.5/20],...
                        'Title','Velocity ratio upper limit ',...
                        'visible','on',...
                        'highlightcolor',handles.bgcolor{3},...
                        'borderwidth',1,...
                        'fontsize',handles.FontSize(2),...
                        'backgroundcolor',handles.bgcolor{3});  
pannelsize = get(handles.uipanel_Limit,'position');
pW=pannelsize(3);
pH=pannelsize(4);  
set(handles.text_b1,...
                        'units', 'points',...
                        'fontUnits','points',...
                        'position',[pW*0.5/10 pH*2.0/10 pW*6/10 pH*3/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','uRatio upper limit: [-]',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'horizontalalignment','left');                         
set(handles.edit_b1,...
                        'units', 'points',...
                        'fontUnits','points',...
                        'position',[pW*7/10 pH*2/10 pW*2.5/10 pH*3.5/10],...
                        'fontsize',handles.FontSize(2),...
                        'string',1,...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','right',...
                        'Enable','on');
%----------------------------------------
% pannel background noise
set(handles.uipanel_BN,...
                        'units', 'points',...
                        'fontUnits','points',...
                        'position',[FigW*0.5/20 FigH*1.75/20 FigW*9/20 FigH*5.5/20],...
                        'Title','Background noise',...
                        'visible','on',...
                        'highlightcolor',handles.bgcolor{3},...
                        'borderwidth',1,...
                        'fontsize',handles.FontSize(2),...
                        'backgroundcolor',handles.bgcolor{3});  
pannelsize = get(handles.uipanel_BN,'position');
pW=pannelsize(3);
pH=pannelsize(4);  
set(handles.text_c1,...
                        'units', 'points',...
                        'fontUnits','points',...
                        'position',[pW*0.5/10 pH*6.0/10 pW*6/10 pH*1.5/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Start time: [s]',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'horizontalalignment','left');                         
set(handles.edit_c1,...
                        'units', 'points',...
                        'fontUnits','points',...
                        'position',[pW*7/10 pH*6/10 pW*2.5/10 pH*1.75/10],...
                        'fontsize',handles.FontSize(2),...
                        'string',0,...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','right',...
                        'Enable','on');
set(handles.text_c2,...
                        'units', 'points',...
                        'fontUnits','points',...
                        'position',[pW*0.5/10 pH*3.5/10 pW*6/10 pH*1.5/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Stop time: [s]',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'horizontalalignment','left');                         
set(handles.edit_c2,...
                        'units', 'points',...
                        'fontUnits','points',...
                        'position',[pW*7/10 pH*3.5/10 pW*2.5/10 pH*1.75/10],...
                        'fontsize',handles.FontSize(2),...
                        'string',1,...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','right',...
                        'Enable','on');                    
set(handles.text_c3,...
                        'units', 'points',...
                        'fontUnits','points',...
                        'position',[pW*0.5/10 pH*1/10 pW*6/10 pH*1.5/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Noise level: [dB]',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'horizontalalignment','left');                         
set(handles.edit_c3,...
                        'units', 'points',...
                        'fontUnits','points',...
                        'position',[pW*7/10 pH*1.0/10 pW*2.5/10 pH*1.75/10],...
                        'fontsize',handles.FontSize(2),...
                        'string',-40,...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','right',...
                        'Enable','on');    
%----------------------------------------
% pannel external forcing
set(handles.uipanel_extForcing,...
                        'units', 'points',...
                        'fontUnits','points',...
                        'position',[FigW*10/20 FigH*7.5/20 FigW*9.5/20 FigH*12.25/20],...
                        'Title','External pressure perturbations ',...
                        'visible','on',...
                        'highlightcolor',handles.bgcolor{3},...
                        'borderwidth',1,...
                        'fontsize',handles.FontSize(2),...
                        'backgroundcolor',handles.bgcolor{3});  
pannelsize = get(handles.uipanel_extForcing,'position');
pW=pannelsize(3);
pH=pannelsize(4);  
set(handles.textExtForcing,...
                        'units', 'points',...
                        'fontUnits','points',...
                        'position',[pW*0.5/10 pH*8.3/10 pW*2/10 pH*0.75/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','External forcing type',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'horizontalalignment','left',...
                        'visible','off');   
set(handles.pop_extForcing_wo,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*0.5/10 pH*8.3/10 pW*9/10 pH*0.75/10],...
                        'fontsize',handles.FontSize(2),...
                        'string',{  'Without external forcing';...
                                    'With external sinusoidal pressure perturbations'},...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','left',...
                        'enable','on',...
                        'value',1); 
set(handles.text_e2,...
                        'units', 'points',...
                        'fontUnits','points',...
                        'position',[pW*0.5/10 pH*7/10 pW*4/10 pH*0.75/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Actuator location:',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'horizontalalignment','left',...
                        'visible','on');  
set(handles.pop_extForcing_pos,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*5/10 pH*7/10 pW*4.5/10 pH*0.75/10],...
                        'fontsize',handles.FontSize(2),...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','left',...
                        'enable','on',...
                        'value',1); 
str{1} = ['inlet of the combustor'];
numSection = length(CI.CD.x_sample);
if numSection > 2
   for ss = 1:numSection-2
      str{2*ss} = ['before the interface ' num2str(ss)];
      str{2*ss+1}   = ['after the interface ' num2str(ss)];
   end
end
str{2*numSection - 2} = ['outlet of the combustor'];
set(handles.pop_extForcing_pos,'string',str);
set(handles.text_e3,...
                        'units', 'points',...
                        'fontUnits','points',...
                        'position',[pW*0.5/10 pH*5.7/10 pW*6/10 pH*0.75/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','a: ',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'horizontalalignment','left');                         
set(handles.edit_e3,...
                        'units', 'points',...
                        'fontUnits','points',...
                        'position',[pW*7/10 pH*5.7/10 pW*2.5/10 pH*0.75/10],...
                        'fontsize',handles.FontSize(2),...
                        'string',100,...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','right',...
                        'Enable','on');
set(handles.text_e4,...
                        'units', 'points',...
                        'fontUnits','points',...
                        'position',[pW*0.5/10 pH*4.4/10 pW*6/10 pH*0.75/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','b: ',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'horizontalalignment','left');                         
set(handles.edit_e4,...
                        'units', 'points',...
                        'fontUnits','points',...
                        'position',[pW*7/10 pH*4.4/10 pW*2.5/10 pH*0.75/10],...
                        'fontsize',handles.FontSize(2),...
                        'string',100,...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','right',...
                        'Enable','on');
set(handles.text_e5,...
                        'units', 'points',...
                        'fontUnits','points',...
                        'position',[pW*0.5/10 pH*3.1/10 pW*6/10 pH*0.75/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Frequency: [Hz]',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'horizontalalignment','left');                         
set(handles.edit_e5,...
                        'units', 'points',...
                        'fontUnits','points',...
                        'position',[pW*7/10 pH*3.1/10 pW*2.5/10 pH*0.75/10],...
                        'fontsize',handles.FontSize(2),...
                        'string',100,...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','right',...
                        'Enable','on');
set(handles.text_e6,...
                        'units', 'points',...
                        'fontUnits','points',...
                        'position',[pW*0.5/10 pH*1.8/10 pW*6/10 pH*0.75/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Start time: [s]',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'horizontalalignment','left');                         
set(handles.edit_e6,...
                        'units', 'points',...
                        'fontUnits','points',...
                        'position',[pW*7/10 pH*1.8/10 pW*2.5/10 pH*0.75/10],...
                        'fontsize',handles.FontSize(2),...
                        'string',0,...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','right',...
                        'Enable','on');
set(handles.text_e7,...
                        'units', 'points',...
                        'fontUnits','points',...
                        'position',[pW*0.5/10 pH*0.5/10 pW*6/10 pH*0.75/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Stop time: [s]',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'horizontalalignment','left');                         
set(handles.edit_e7,...
                        'units', 'points',...
                        'fontUnits','points',...
                        'position',[pW*7/10 pH*0.5/10 pW*2.5/10 pH*0.75/10],...
                        'fontsize',handles.FontSize(2),...
                        'string',1,...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','right',...
                        'Enable','on');
%---------------------------
% pannel AOC                   
set(handles.uipanel_AOC,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[FigW*0.5/20 FigH*0/20 FigW*19/20 FigH*1.5/20],...
                        'Title','',...
                        'visible','on',...
                        'highlightcolor',handles.bgcolor{3},...
                        'borderwidth',1,...
                        'fontsize',handles.FontSize(2),...
                        'backgroundcolor',handles.bgcolor{3}); 
pannelsize=get(handles.uipanel_AOC,'position');                    
pW=pannelsize(3);
pH=pannelsize(4);                
set(handles.pb_OK,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*0.5/10 pH*2/10 pW*2.5/10 pH*6/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','OK',...
                        'backgroundcolor',handles.bgcolor{3});
set(handles.pb_Help,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*3.75/10 pH*2/10 pW*2.5/10 pH*6/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Help',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'enable','on');
set(handles.pb_Cancel,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*7.0/10 pH*2/10 pW*2.5/10 pH*6/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Cancel',...
                        'backgroundcolor',handles.bgcolor{3});
% --------------------------------
guidata(hObject, handles);
handles     = guidata(hObject);
%specific configuration for G-Equation
% if CI.FM.NL.style == 4;
%  %set(handles.uipanel_Speed,'visible','off');
%  set(handles.uipanel_Limit,'visible','off');
%   
% set(handles.edit_d2,...
%                         'string',1,...
%                         'Enable','off');        
% end
if ~isempty(CI.CD.indexHP)    % if there are heat perturbations    
   if ~isempty(find(CI.FM.indexFM == 4))  % G-equation
    % ************* need further change **************
    %set(handles.uipanel_Speed,'visible','off');
    set(handles.uipanel_Limit,'visible','off');
  
    set(handles.edit_d2,...
                        'string',1,...
                        'Enable','off');  
% ************* need further change **************
   end
end
% --------------------------------
guidata(hObject, handles);
if CI.IsRun.GUI_TD_Para_Config == 1
     Fcn_initialize_input_values(hObject)
else
    Fcn_pre_evaluation(hObject)
    handles     = guidata(hObject);
    Fcn_pop_with_without_external_perturbation_update(hObject, handles, 1)
end
assignin('base','CI',CI)
if isempty(CI.CD.indexHP)
   set(handles.uipanel_Limit, 'visible','off');  % in case without flame, this pannel is not visible 
end
guidata(hObject, handles);
%

function Fcn_pre_evaluation(varargin)
hObject     = varargin{1};
handles     = guidata(hObject);
global CI
uRatioMax   = str2double(get(handles.edit_b1,'String'));
dt          = str2double(get(handles.edit_a3,'String'));
tEndRaw     = str2double(get(handles.edit_a2,'String'));
nGapRatio   = 1;                                                % default value
Fcn_TD_INI_samples_information(dt,tEndRaw,uRatioMax,nGapRatio)
set(handles.edit_a1, 'string', num2str(CI.TD.tauMin));
set(handles.edit_d1, 'string', num2str(CI.TD.nGapMax));
set(handles.edit_c1, 'string', num2str(CI.TD.tSpTotal(1)));
set(handles.edit_c2, 'string', num2str(CI.TD.tSpTotal(end)));
set(handles.edit_e6, 'string', num2str(CI.TD.tSpTotal(1)));
set(handles.edit_e7, 'string', num2str(CI.TD.tSpTotal(end)));
guidata(hObject, handles);


function Fcn_initialize_input_values(varargin)
hObject     = varargin{1};
handles     = guidata(hObject);
global CI
set(handles.edit_a1, 'string', num2str(CI.TD.tauMin));
set(handles.edit_a2, 'string', num2str(CI.TD.tEndRaw));
set(handles.edit_a3, 'string', num2str(CI.TD.dt));
set(handles.edit_d1, 'string', num2str(CI.TD.nGapMax));
set(handles.edit_d2, 'string', num2str(CI.TD.nGap));
set(handles.edit_b1, 'string', num2str(CI.TD.uRatioMax));
set(handles.edit_c1, 'string', num2str(CI.TD.NoiseInfo.t(1)));
set(handles.edit_c2, 'string', num2str(CI.TD.NoiseInfo.t(2)));
set(handles.edit_c3, 'string', num2str(CI.TD.NoiseInfo.level));
set(handles.pop_extForcing_wo, 'value', CI.TD.ExtForceInfo.indexForcing);
set(handles.pop_extForcing_pos, 'value', CI.TD.ExtForceInfo.indexPos);
set(handles.edit_e3, 'string', num2str(CI.TD.ExtForceInfo.a));
set(handles.edit_e4, 'string', num2str(CI.TD.ExtForceInfo.b));
set(handles.edit_e5, 'string', num2str(CI.TD.ExtForceInfo.f));
set(handles.edit_e6, 'string', num2str(CI.TD.ExtForceInfo.t(1)));
set(handles.edit_e7, 'string', num2str(CI.TD.ExtForceInfo.t(2)));
guidata(hObject, handles);
handles     = guidata(hObject);
Fcn_pop_with_without_external_perturbation_update(hObject, handles, CI.TD.ExtForceInfo.indexForcing)
guidata(hObject, handles);


% --- Executes on button press in pb_OK.
function pb_OK_Callback(hObject, eventdata, handles)
% hObject    handle to pb_OK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global CI
dt                  = str2double(get(handles.edit_a3,'String'));
tEndRaw             = str2double(get(handles.edit_a2,'String'));
uRatioMax           = str2double(get(handles.edit_b1,'String'));
nGap                = str2double(get(handles.edit_d2,'String'));
nGapRatio           = nGap./CI.TD.nGapMax;
NoiseInfo.t(1)      = str2double(get(handles.edit_c1,'String'));
NoiseInfo.t(2)      = str2double(get(handles.edit_c2,'String'));
NoiseInfo.level     = str2double(get(handles.edit_c3,'String'));
ExtForceInfo.indexForcing   = get(handles.pop_extForcing_wo,'value');    
ExtForceInfo.indexPos       = get(handles.pop_extForcing_pos,'value');   
ExtForceInfo.a      = str2double(get(handles.edit_e3,'String'));  
ExtForceInfo.b      = str2double(get(handles.edit_e4,'String')); 
ExtForceInfo.f      = str2double(get(handles.edit_e5,'String')); 
ExtForceInfo.t(1)   = str2double(get(handles.edit_e6,'String'));
ExtForceInfo.t(2)   = str2double(get(handles.edit_e7,'String'));
Fcn_TD_INI(dt,tEndRaw,uRatioMax,nGapRatio,NoiseInfo,ExtForceInfo)
CI.IsRun.GUI_TD_Para_Config = 1;
assignin('base','CI',CI)
%
%
main = handles.MainGUI;                     % get the handle of OSCILOS_long
if(ishandle(main))
    mainHandles = guidata(main);            %
   set(mainHandles.TD_SIM, 'enable' , 'on') 
end
guidata(hObject, handles);
delete(handles.figure);


function Fcn_pop_with_without_external_perturbation_update(hObject, handles, PopVal)
handles.ObjEditVisible_edit_e           = findobj('-regexp','Tag','edit_e');
handles.ObjEditVisible_text_e           = findobj('-regexp','Tag','text_e');
switch PopVal
    case 1
        set(handles.ObjEditVisible_edit_e,'visible','off');
        set(handles.pop_extForcing_pos,'visible','off');
        set(handles.ObjEditVisible_text_e,'visible','off');
    case 2
        set(handles.ObjEditVisible_edit_e,'visible','on');
        set(handles.pop_extForcing_pos,'visible','on');
        set(handles.ObjEditVisible_text_e,'visible','on');
end
guidata(hObject, handles);


% --- Executes on button press in pb_Help.
function pb_Help_Callback(hObject, eventdata, handles)
% hObject    handle to pb_Help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
string={...
            'This GUI is used to configure the parameters for the simulation.';...
            '1. The time step must be smaller than the minimum time delay.'
            '2. The number of time steps should not be larger than the suggested value.'};
helpdlg(string,'Help')

% --- Executes on button press in pb_Cancel.
function pb_Cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pb_Cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(handles.figure);


% --- Executes during object creation, after setting all properties.
function edit_a1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_a1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%
function edit_a1_Callback(hObject, eventdata, handles)
% hObject    handle to edit_a1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
function edit_a2_Callback(hObject, eventdata, handles)
% hObject    handle to edit_a2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
datEdit = str2double(get(hObject, 'String'));
ValDefault = 1; 
if isnan(datEdit) || ~isreal(datEdit) ||datEdit < 0
    set(hObject, 'String', ValDefault);
    errordlg('Input must be a non-negative real number','Error');
    % when the input is not a number, it is set to the default value
end
%
% --- Executes during object creation, after setting all properties.
function edit_a2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_a2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%
function edit_a3_Callback(hObject, eventdata, handles)
% hObject    handle to edit_a3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
datEdit = str2double(get(hObject, 'String'));
ValDefault = 1e-4; 
if isnan(datEdit) || ~isreal(datEdit) ||datEdit < 0
    set(hObject, 'String', ValDefault);
    errordlg('Input must be a non-negative real number','Error');
    % when the input is not a number, it is set to the default value
end
tdMin = str2double(get(handles.edit_a1,'String'));
if datEdit > tdMin
    set(hObject, 'String',tdMin);
    errordlg('Time step must be smaller than the minimum time delay','Error');
end
Fcn_pre_evaluation(hObject)
%
% --- Executes during object creation, after setting all properties.
function edit_a3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_a3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%
function edit_d1_Callback(hObject, eventdata, handles)
% hObject    handle to edit_d1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
% --- Executes during object creation, after setting all properties.
function edit_d1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_d1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%
function edit_d2_Callback(hObject, eventdata, handles)
% hObject    handle to edit_d2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
datEdit = str2double(get(hObject, 'String'));
ValDefault = 1; 
if isnan(datEdit) || ~isreal(datEdit) ||datEdit <= 0 ||rem(datEdit,1) ~= 0
    set(hObject, 'String', ValDefault);
    errordlg('Input must be a positive integer','Error');
    % when the input is not a number, it is set to the default value
end
nGapMax = str2double(get(handles.edit_d1,'String'));
if datEdit > nGapMax
    set(hObject, 'String',nGapMax);
    errordlg('This value should not be larger than the maximum value','Error');
end
%
% --- Executes during object creation, after setting all properties.
function edit_d2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_d2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%
function edit_b1_Callback(hObject, eventdata, handles)
% hObject    handle to edit_b1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
datEdit = str2double(get(hObject, 'String'));
ValDefault = 1; 
if isnan(datEdit) || ~isreal(datEdit) ||datEdit <= 0
    set(hObject, 'String', ValDefault);
    errordlg('Input must be a positive real number','Error');
    % when the input is not a number, it is set to the default value
end
%
% --- Executes during object creation, after setting all properties.
function edit_b1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_b1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%
function edit_c1_Callback(hObject, eventdata, handles)
% hObject    handle to edit_c1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
datEdit = str2double(get(hObject, 'String'));
ValDefault = 0; 
if isnan(datEdit) || ~isreal(datEdit) ||datEdit <0
    set(hObject, 'String', ValDefault);
    errordlg('Input must be a non-negative real number','Error');
    % when the input is not a number, it is set to the default value
end
% --- Executes during object creation, after setting all properties.
function edit_c1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_c1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%
function edit_c2_Callback(hObject, eventdata, handles)
% hObject    handle to edit_c2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
datEdit = str2double(get(hObject, 'String'));
ValDefault = str2double(get(handles.edit_a2, 'String')); 
if isnan(datEdit) || ~isreal(datEdit) ||datEdit <0
    set(hObject, 'String', ValDefault);
    errordlg('Input must be a non-negative real number','Error');
    % when the input is not a number, it is set to the default value
end
if datEdit > ValDefault
    set(hObject, 'String', ValDefault);
    errordlg('Input cannot be larger than the stop time','Error');
end
    
% --- Executes during object creation, after setting all properties.
function edit_c2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_c2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%
function edit_c3_Callback(hObject, eventdata, handles)
% hObject    handle to edit_c3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
datEdit = str2double(get(hObject, 'String'));
ValDefault = -40; 
if isnan(datEdit) || ~isreal(datEdit) 
    set(hObject, 'String', ValDefault);
    errordlg('Input must be a real number','Error');
    % when the input is not a number, it is set to the default value
end
% --- Executes during object creation, after setting all properties.
function edit_c3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_c3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%


% --- Executes when user attempts to close figure.
function figure_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
delete(hObject);


% --- Executes during object creation, after setting all properties.
function figure_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on selection change in pop_extForcing_wo.
function pop_extForcing_wo_Callback(varargin)
hObject     = varargin{1};
handles     = guidata(hObject);
popVal      = get(handles.pop_extForcing_wo, 'value');
Fcn_pop_with_without_external_perturbation_update(hObject, handles, popVal)


% --- Executes during object creation, after setting all properties.
function pop_extForcing_wo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_extForcing_wo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in pop_extForcing_pos.
function pop_extForcing_pos_Callback(hObject, eventdata, handles)
% hObject    handle to pop_extForcing_pos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pop_extForcing_pos contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_extForcing_pos


% --- Executes during object creation, after setting all properties.
function pop_extForcing_pos_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_extForcing_pos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_e3_Callback(hObject, eventdata, handles)
% hObject    handle to edit_e3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_e3 as text
%        str2double(get(hObject,'String')) returns contents of edit_e3 as a double
datEdit = str2double(get(hObject, 'String'));
ValDefault = 100; 
if isnan(datEdit) || ~isreal(datEdit) 
    set(hObject, 'String', ValDefault);
    errordlg('Input must be a real number','Error');
    % when the input is not a number, it is set to the default value
end

% --- Executes during object creation, after setting all properties.
function edit_e3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_e3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_e4_Callback(hObject, eventdata, handles)
% hObject    handle to edit_e4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_e4 as text
%        str2double(get(hObject,'String')) returns contents of edit_e4 as a double
datEdit = str2double(get(hObject, 'String'));
ValDefault = 100; 
if isnan(datEdit) || ~isreal(datEdit) 
    set(hObject, 'String', ValDefault);
    errordlg('Input must be a real number','Error');
    % when the input is not a number, it is set to the default value
end

% --- Executes during object creation, after setting all properties.
function edit_e4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_e4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_e5_Callback(hObject, eventdata, handles)
% hObject    handle to edit_e5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_e5 as text
%        str2double(get(hObject,'String')) returns contents of edit_e5 as a double
datEdit = str2double(get(hObject, 'String'));
ValDefault = 100; 
if isnan(datEdit) || ~isreal(datEdit) 
    set(hObject, 'String', ValDefault);
    errordlg('Input must be a real number','Error');
    % when the input is not a number, it is set to the default value
end

% --- Executes during object creation, after setting all properties.
function edit_e5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_e5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_e6_Callback(hObject, eventdata, handles)
% hObject    handle to edit_e6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_e6 as text
%        str2double(get(hObject,'String')) returns contents of edit_e6 as a double
datEdit = str2double(get(hObject, 'String'));
ValDefault = 0; 
if isnan(datEdit) || ~isreal(datEdit) 
    set(hObject, 'String', ValDefault);
    errordlg('Input must be a real number','Error');
    % when the input is not a number, it is set to the default value
end

% --- Executes during object creation, after setting all properties.
function edit_e6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_e6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_e7_Callback(hObject, eventdata, handles)
% hObject    handle to edit_e7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_e7 as text
%        str2double(get(hObject,'String')) returns contents of edit_e7 as a double
datEdit = str2double(get(hObject, 'String'));
ValDefault = 1; 
if isnan(datEdit) || ~isreal(datEdit) 
    set(hObject, 'String', ValDefault);
    errordlg('Input must be a real number','Error');
    % when the input is not a number, it is set to the default value
end

% --- Executes during object creation, after setting all properties.
function edit_e7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_e7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
