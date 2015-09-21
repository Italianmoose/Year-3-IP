function varargout = GUI_INI_GEQU(varargin)
% GUI_INI_GEQU MATLAB code for GUI_INI_GEQU.fig
%      GUI_INI_GEQU, by itself, creates a new GUI_INI_GEQU or raises the existing
%      singleton*.
%
%      H = GUI_INI_GEQU returns the handle to a new GUI_INI_GEQU or the handle to
%      the existing singleton*.
%
%      GUI_INI_GEQU('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_INI_GEQU.M with the given input arguments.
%
%      GUI_INI_GEQU('Property','Value',...) creates a new GUI_INI_GEQU or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_INI_GEQU_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_INI_GEQU_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_INI_GEQU

% Last Modified by GUIDE v2.5 26-Nov-2014 16:20:56

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_INI_GEQU_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_INI_GEQU_OutputFcn, ...
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


% --- Executes just before GUI_INI_GEQU is made visible.
function GUI_INI_GEQU_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_INI_GEQU (see VARARGIN)
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
% --------------------------
% Global initialisation : for variables if the GUI has never been called
% yet
if handles.HP.IsRun == 0
    Fcn_GUI_INI_FM_GEQU_Var_initilisation(hObject, eventdata, handles)
    handles = guidata(hObject);
    guidata(hObject, handles);
end

% Initialization
GUI_INI_GEQU_Initialization(hObject, eventdata, handles);
%         end
guidata(hObject, handles);
handles.output = hObject;
guidata(hObject, handles);
        
% UIWAIT makes GUI_INI_GEQU wait for user response (see UIRESUME)
% uiwait(handles.figure1);
function Fcn_GUI_INI_FM_GEQU_Var_initilisation(varargin)
hObject = varargin{1};
handles = guidata(hObject);
global CI
% These values serve as defaults
handles.HP.NL.style      = 4; % G-EQuation (Williams 1985)
handles.HP.GEQU.nb_points = 35; % Number of points used for discretisation along r
handles.HP.GEQU.rb = CI.CD.r_sample(CI.CD.indexHP(handles.HP_num)); % If there are multiple heat zones in the duct, this is a vector
handles.HP.GEQU.ra   = handles.HP.GEQU.rb/2; % in m, also a vector
handles.HP.GEQU.rho1 = CI.TP.rho_mean(1,max(CI.CD.indexHP(handles.HP_num) - 1,1));% This is a vector if there are multple flame. The max function is required is the flame is at the begining of the duct.
handles.HP.GEQU.U1 = CI.TP.u_mean(1,max(CI.CD.indexHP(handles.HP_num) - 1,1));% This is a vector if there are multple flame. The max function is required is the flame is at the begining of the duct.
handles.HP.GEQU.SU   = handles.HP.GEQU.U1 * 0.088; % in m/s, this is a vector if there are multple flame.
handles.HP.GEQU.time_integration = 3; % Set to adam bashforth
handles.HP.GEQU.dowst_of_heat_lengths = CI.CD.x_sample(min(CI.CD.indexHP(handles.HP_num) + 1,end)) - CI.CD.x_sample(CI.CD.indexHP(handles.HP_num));
handles.HP.GEQU.tau_f_factor = 0.42;
guidata(hObject, handles);
%
%-------------------------------------------------
%
function GUI_INI_GEQU_Initialization(varargin) % graphic initialisation
hObject = varargin{1};
handles = guidata(hObject);    
set(0, 'units', 'points');
screenSize  = get(0, 'ScreenSize');                 % get the screen size
sW          = handles.sW;                           % screen width
sH          = handles.sH ;                          % screen height
FigW=sW*4/5;                                        % window width
FigH=sH*3/4;                                        % window height
set(handles.figure,     'units', 'points',...
                        'position',[(screenSize(3)-FigW)./2 (screenSize(4)-FigH)./2 FigW FigH],...
                        'name','G-Equation configuration',...
                        'color',handles.bgcolor{3});
%--------------------------------------------------------------------------
% pannel axes
set(handles.uipanel_Axes,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[FigW*0.5/20 FigH*1.75/20 FigW*19/20 FigH*13.25/20],... % Left Bottom Width Height
                        'Title','Steady Flame Shapes',...
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
                        'position',[pW*1/10 pH*1.5/10 pW*8/10 pH*7.5/10],...
                        'fontsize',handles.FontSize(1),...
                        'color',handles.bgcolor{1},...
                        'box','on');  
guidata(hObject, handles);
%----------------------------------------               
% pannel GEQU
set(handles.uipanel_GEQU,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[FigW*0.5/20 FigH*15.25/20 FigW*19/20 FigH*4.5/20],...
                        'Title','Enter flame properties',...
                        'visible','on',...
                        'highlightcolor',handles.bgcolor{3},...
                        'borderwidth',1,...
                        'fontsize',handles.FontSize(2),...
                        'backgroundcolor',handles.bgcolor{3}); 
pannelsize = get(handles.uipanel_GEQU,'position');
pW = pannelsize(3);
pH = pannelsize(4);                      
set(handles.text_GEQU_a1,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*0.5/10 pH*5.25/10 pW*3.0/10 pH*1.5/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Laminar burning vel. S_u [m/s]',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'horizontalalignment','left');                         
set(handles.edit_GEQU_a1,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*3.25/10 pH*5.25/10 pW*1.5/10 pH*2/10],...
                        'fontsize',handles.FontSize(2),...
                        'string',num2str(handles.HP.GEQU.SU),...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','right',...
                        'Enable','on');

set(handles.text_GEQU_a2,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*0.5/10 pH*2.0/10 pW*3.0/10 pH*1.5/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Flame time delay factor tau_factor [-]',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'horizontalalignment','left');                   
set(handles.edit_GEQU_a2,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*3.25/10 pH*2.0/10 pW*1.5/10 pH*2/10],...
                        'fontsize',handles.FontSize(2),...
                        'string',num2str(handles.HP.GEQU.tau_f_factor) ,...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','right',...
                        'Enable','on');  
set(handles.text_GEQU_a3,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*5.25/10 pH*5.25/10 pW*3.0/10 pH*1.5/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Space discret. (n points) [-]',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'horizontalalignment','left');                   
set(handles.edit_GEQU_a3,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*8.0/10 pH*5.25/10 pW*1.5/10 pH*2/10],...
                        'fontsize',handles.FontSize(2),...
                        'string',num2str(handles.HP.GEQU.nb_points),...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','right',...
                        'Enable','on'); 
set(handles.text_GEQU_a4,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*5.25/10 pH*2.0/10 pW*3/10 pH*1.5/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Flame holder radius r_a [m]',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'horizontalalignment','left');                   
set(handles.edit_GEQU_a4,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*8.0/10 pH*2.0/10 pW*1.5/10 pH*2/10],...
                        'fontsize',handles.FontSize(2),...
                        'string',num2str(handles.HP.GEQU.ra),...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','right',...
                        'Enable','on');  
%----------------------------------------
%
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
set(handles.pb_Apply,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*0.4/10 pH*2/10 pW*2.0/10 pH*6/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Plot figure',...
                        'backgroundcolor',handles.bgcolor{3});
set(handles.pb_SaveFig,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*2.8/10 pH*2/10 pW*2.0/10 pH*6/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Save figure',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'enable','off');
set(handles.pb_OK,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*5.2/10 pH*2/10 pW*2.0/10 pH*6/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','OK',...
                        'backgroundcolor',handles.bgcolor{3});
set(handles.pb_Cancel,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*7.6/10 pH*2/10 pW*2.0/10 pH*6/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Cancel',...
                        'backgroundcolor',handles.bgcolor{3});
%---------------------------------------
handles.obj_edit_GEQU  = findobj('-regexp','Tag','edit_GEQU');
handles.obj_text_GEQU  = findobj('-regexp','Tag','text_GEQU');
% default visible settings
set(handles.obj_edit_GEQU,          'visible','on')
set(handles.obj_text_GEQU,          'visible','on')
%
%---------------------------------------
guidata(hObject, handles);
%Fcn_GUI_INI_FM_GEQU_Value_update(hObject);
handles = guidata(hObject);
guidata(hObject, handles);
%Fcn_GUI_INI_FM_Plot_GEQU_Shape(hObject);
handles = guidata(hObject);
guidata(hObject, handles);
%
% ------------------------------------------------------------------------
% 
function Fcn_GUI_INI_FM_GEQU_Value_update(varargin)
hObject = varargin{1};
handles = guidata(hObject);
%global CI

% Get data from GUI cells
handles.HP.GEQU.SU = str2num(get(handles.edit_GEQU_a1,'String')); % vector if multiple flames. str2num required here to be able to deal with vector inputs
handles.HP.GEQU.tau_f_factor = str2num(get(handles.edit_GEQU_a2,'String'));
handles.HP.GEQU.nb_points = str2num(get(handles.edit_GEQU_a3,'String')); % str2num required here to be able to deal with vector inputs
handles.HP.GEQU.ra = str2num(get(handles.edit_GEQU_a4,'String'));
% Compute important values
handles.HP.GEQU.area_ratio = 1.0 -(handles.HP.GEQU.ra./handles.HP.GEQU.rb).^2; % vector if there are multiple flames
handles.HP.GEQU.Ugs = Fcn_TD_Gequ_calc_ugutter( handles.HP.GEQU.U1,handles.HP.GEQU.area_ratio,0,0 ); % vector if there are multiple flames
handles.HP.GEQU.tau_f = handles.HP.GEQU.tau_f_factor .* handles.HP.GEQU.dowst_of_heat_lengths./handles.HP.GEQU.Ugs; % vector of time delays for everyflame in the duct

for runner = 1:length(handles.HP.GEQU.ra)
    handles.HP.GEQU.y_vec(runner,:) = linspace(handles.HP.GEQU.ra(runner),handles.HP.GEQU.rb(runner),handles.HP.GEQU.nb_points(runner)); % currently the nb of points for all flames nees to be the same
end
handles.HP.GEQU.area_ratio = 1.0 -(handles.HP.GEQU.ra./handles.HP.GEQU.rb).^2; % vector if there are multiple flames
handles.HP.GEQU.Ugs = Fcn_TD_Gequ_calc_ugutter( handles.HP.GEQU.U1,handles.HP.GEQU.area_ratio,0,0 ); % vector if there are multiple flames
handles.HP.GEQU.xi_steady       = ...
    Fcn_TD_Gequ_steady_flame( handles.HP.GEQU.Ugs,handles.HP.GEQU.SU,handles.HP.GEQU.y_vec ); % If there are multiple flame, this is a matrix (lines are flames in different sections, columns come along r)
handles.HP.GEQU.xi = handles.HP.GEQU.xi_steady; % initialise xi value
handles.HP.GEQU.bashforth_data = zeros(length(handles.HP.GEQU.rb),max(handles.HP.GEQU.nb_points),3); % Three dimentional matrix for three step adam bashforth

% Set the flame transfer function numerator and denominator to
% 1, as they are used for the Green's function initilisation,
handles.HP.FTF.num = 1;
handles.HP.FTF.den = 1;
        
guidata(hObject, handles);
%
% ------------------------------------------------------------------------
% 
% --- Executes on button press in pb_OK.
function pb_OK_Callback(hObject, eventdata, handles)
% hObject    handle to pb_OK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Fcn_GUI_INI_FM_GEQU_Value_update(hObject, eventdata, handles)
handles = guidata(hObject);
guidata(hObject, handles);
Fcn_GUI_INI_FM_Plot_GEQU_Shape(hObject, eventdata, handles)
handles = guidata(hObject);
handles.HP.IsRun = 1;    % index == 1 to show that this program has been run
guidata(hObject, handles);
global CI
CI.FM.HP{handles.HP_num}        = handles.HP;
CI.FM.indexFM(handles.HP_num)   = handles.indexFM;
assignin('base','CI',CI); 
handles = guidata(hObject);
guidata(hObject, handles);
Fcn_GUI_INI_FM_GEQU_Update_Data(hObject, eventdata, handles)
delete(handles.figure);
%
% ------------------------------------------------------------------------
% 
% --- Executes on button press in pb_Apply.
function pb_Apply_Callback(hObject, eventdata, handles)
% hObject    handle to pb_Apply (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Fcn_GUI_INI_FM_GEQU_Value_update(hObject, eventdata, handles)
Fcn_GUI_INI_FM_Plot_GEQU_Shape(hObject, eventdata, handles)
handles = guidata(hObject);
set(handles.pb_SaveFig,'enable','on');
guidata(hObject, handles);
%
% ------------------------------------------------------------------------
% 
% --- Update the data when clicking 'OK' or 'Apply'
function Fcn_GUI_INI_FM_GEQU_Update_Data(hObject, eventdata, handles)
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
% ------------------------------------------------------------------------
% 
% --- plot the combustor shape
function Fcn_GUI_INI_FM_Plot_GEQU_Shape(varargin)
hObject = varargin{1};
handles = guidata(hObject);
global CI
% 
try
x_sample    = CI.CD.x_sample;                
r_sample    = CI.CD.r_sample;
%-------------------------------------
W           = abs(x_sample(end) - x_sample(1));             % Length of the combustor
H           = 2*max(r_sample);                              % Diameter of the combustor
plot_ratio  = 1.5;                                          % Ratio of the axes limit to the combustor dimension
axes_W      = plot_ratio*W;                                 % axes x width
axes_H      = 2*plot_ratio*H;                               % axes y width        
x_min       = x_sample(1) - (axes_W-W)./2;                  % axes x min
y_min       = -axes_H./2;                                   % axes y min
%--------------------------------------
hAxes = handles.axes1;
cla(hAxes);                             % clear the current axes 
axis(hAxes);
hold on
%--------------------------------------
% plot the approximate profile of the combustor which consisting of several
% sections 
for s=1:length(x_sample)-1
    rectangle(  'Position',1000*[x_sample(s),-r_sample(s),...
                 x_sample(s+1)-x_sample(s),2*r_sample(s)],...
                'Curvature',[0,0],'LineWidth',1,'LineStyle','-');
end

%--------------------------------------
% Plot the flame holder
flame_holder_length = 0.01; % in m
 % This loop isn't erally required, as each time we are only plotting one flame. 
 % Still its good to keep the code general
for s = 1:length(CI.CD.indexHP(handles.HP_num))
rectangle(  'Position',1000*[x_sample(CI.CD.indexHP(handles.HP_num(s))) - flame_holder_length,-handles.HP.GEQU.ra(s),...
                flame_holder_length,...
                2*handles.HP.GEQU.ra(s)],...
                'Curvature',[0,0],'LineWidth',1,'LineStyle','-','facecolor','b');
end

 % This loop isn't erally required, as each time we are only plotting one flame. 
 % Still its good to keep the code general
%plot the flame sheet
for s = 1:length(CI.CD.indexHP(handles.HP_num))
plot(1000 * handles.HP.GEQU.xi(s,:) + 1000 *x_sample(CI.CD.indexHP(handles.HP_num(s))), 1000* handles.HP.GEQU.y_vec(s,:),'-r',...
    1000 * handles.HP.GEQU.xi(s,:) + 1000 *x_sample(CI.CD.indexHP(handles.HP_num(s))), -1000 * handles.HP.GEQU.y_vec(s,:),'-r')
end

set(hAxes,'xlim',1000*[x_min, x_min+axes_W]);
set(hAxes,'box','on');
grid on
set(hAxes,'ylim',1000*[y_min, y_min+axes_H]);
text(1000*(x_sample(1)-W./8), 0, 'inlet','FontSize',15,...
    'interpreter','latex','HorizontalAlignment','center');
text(1000*(x_sample(end)+W.*5/32), 0, 'outlet','FontSize',15,...
    'interpreter','latex','HorizontalAlignment','center');
text(1000*x_sample(index_flame), 1000*2*r_sample(index_flame),...
    'flame','FontSize',15,...
    'interpreter','latex','HorizontalAlignment','center');
xlabel(hAxes,'x~ [mm]','Color','k','Interpreter','LaTex');
ylabel(hAxes,'r~ [mm]','Color','k','Interpreter','LaTex');   
catch
end
%
% ------------------------------------------------------------------------
% 
% --- Executes on button press in pb_SaveFig.
function pb_SaveFig_Callback(varargin)
global CI
hObject = varargin{1};
handles = guidata(hObject);
Fig                 = figure;
set(Fig,'units','points')
posFig              = get(handles.figure,'position');
copyobj(handles.axes1, Fig);

hAxes               = get(Fig,'children');
set(hAxes(1),       'units','points',...
                    'position',[60 60 200 150],...
                    'ActivePositionProperty','position')
pos1                = get(hAxes(1),'position');

posAxesOuter = [0 0 340 400];
colormap(hot);
hcb = colorbar;
set(hcb,'parent',Fig)
set(hcb,'Fontsize',handles.FontSize(1),'box','on','Unit','points')
set(hcb,'ActivePositionProperty','position');
set(hcb,'position',[260,60,10,300]);
set(hAxes(1),       'position',pos1)
set(hcb,'position',[260,60,10,300]);


set(Fig,        'units','points',...
                'position', [posFig(1)+0.5*posFig(3)-0.5*posAxesOuter(3),...
                            posFig(2)+0.5*posFig(4)-0.5*posAxesOuter(4),...
                            posAxesOuter(3:4)])       
%
% ------------------------------------------------------------------------
% 
% --- Outputs from this function are returned to the command line.
function varargout = GUI_INI_GEQU_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
%
% ------------------------------------------------------------------------
% 
% --- Executes on button press in pb_Cancel.
function pb_Cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pb_Cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(handles.figure);
%
% ------------------------------------------------------------------------
%
% --- Executes when user attempts to close figure.
function figure_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(hObject);
