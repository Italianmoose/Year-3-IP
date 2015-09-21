function varargout = GUI_TD_Cal_JLI_AMorgans(varargin)
% GUI_TD_Cal_JLI_AMorgans MATLAB code for GUI_TD_Cal_JLI_AMorgans.fig
%      GUI_TD_Cal_JLI_AMorgans, by itself, creates a new GUI_TD_Cal_JLI_AMorgans or raises the existing
%      singleton*.
%
%      H = GUI_TD_Cal_JLI_AMorgans returns the handle to a new GUI_TD_Cal_JLI_AMorgans or the handle to
%      the existing singleton*.
%
%      GUI_TD_Cal_JLI_AMorgans('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_TD_Cal_JLI_AMorgans.M with the given input arguments.
%
%      GUI_TD_Cal_JLI_AMorgans('Property','Value',...) creates a new GUI_TD_Cal_JLI_AMorgans or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_TD_Cal_JLI_AMorgans_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_TD_Cal_JLI_AMorgans_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_TD_Cal_JLI_AMorgans

% Last Modified by GUIDE v2.5 05-Jan-2015 14:09:19

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_TD_Cal_JLI_AMorgans_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_TD_Cal_JLI_AMorgans_OutputFcn, ...
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
% -------------------------------------------------------------------------
%          
% --- Executes just before GUI_TD_Cal_JLI_AMorgans is made visible.
function GUI_TD_Cal_JLI_AMorgans_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_TD_Cal_JLI_AMorgans (see VARARGIN)
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
            mainHandles         = guidata(handles.MainGUI);
            % background colors
            handles.bgcolor     = mainHandles.bgcolor;
            % fontsize
            handles.FontSize    = mainHandles.FontSize;
            %
            handles.sW          = mainHandles.sW;
            handles.sH          = mainHandles.sH;
            handles.indexApp    = 0;
            % Update handles structure
            guidata(hObject, handles);
            handles.output = hObject;
            guidata(hObject, handles);
            % Initialization
            GUI_Pannel_Initialization(hObject, eventdata, handles)
        end
        % Update handles structure
        guidata(hObject, handles);
        if dontOpen
           disp('-----------------------------------------------------');
           disp('This is a subprogram. It cannot be run independently.') 
           disp('Please load the program "OSCILOS_long'' from the ')
           disp('parent directory!')
           disp('-----------------------------------------------------');
        end
    case 1
        global CI
        handles = Fcn_GUI_default_configuration(handles);
        handles.output = hObject;
        guidata(hObject, handles);
        CI.IsRun.GUI_TD_Cal_JLI_AMorgans = 0;
        assignin('base','CI',CI);                   % save the current information to the works
        guidata(hObject, handles);  
        GUI_Pannel_Initialization(hObject, eventdata, handles)
end
%
% -------------------------------------------------------------------------
%
function GUI_Pannel_Initialization(varargin)
hObject = varargin{1};
handles = guidata(hObject);
global CI        
CI.TD.JLI.iteration.tol = 1e-5;
CI.TD.JLI.iteration.num = 50;
CI.TD.nPeriod           = 10;
CI.TD.RatioGapPadding   = 1;
assignin('base','CI',CI);                   % save the current information to the workspace
%
set(0, 'units', 'points');
screenSize  = get(0, 'ScreenSize');                             % get the screen size
sW          = handles.sW;                                       % screen width
sH          = handles.sH;                                       % screen height
FigW        = sW.*11/20;                                          % window width
FigH        = sH.*3/5;                                          % window height
set(handles.figure,     'units', 'points',...
                        'position',[(screenSize(3)-FigW)./2 (screenSize(4)-FigH)./2 FigW FigH],...
                        'name','Calculation monitor',...
                        'color',handles.bgcolor{3});
%----------------------------------------
% pannel axes
set(handles.uipanel_axes,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[FigW*0.5/20 FigH*7.0/20 FigW*19/20 FigH*12.5/20],...
                        'Title','',...
                        'visible','on',...
                        'highlightcolor',handles.bgcolor{3},...
                        'borderwidth',1,...
                        'fontsize',handles.FontSize(2),...
                        'backgroundcolor',handles.bgcolor{3});  
pannelsize=get(handles.uipanel_axes,'position');
pW=pannelsize(3);
pH=pannelsize(4);                
set(handles.axes1,      'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*1.5/10 pH*2.0/10 pW*7.5/10 pH*7/10],...
                        'fontsize',handles.FontSize(1),...
                        'color',handles.bgcolor{1},...
                        'box','on');                     
guidata(hObject, handles);
set(handles.pop_type,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*0.5/10 pH*9/10 pW*5.5/10 pH*0.5/10],...
                        'fontsize',handles.FontSize(2),...
                        'string',{ },...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','left',...
                        'enable','on',...
                        'value',1,...
                        'visible','off');  

%----------------------------------------
% pannel input
set(handles.uipanel_Output,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[FigW*0.5/20 FigH*2.5/20 FigW*19/20 FigH*4.25/20],...
                        'Title','',...
                        'visible','on',...
                        'highlightcolor',handles.bgcolor{3},...
                        'borderwidth',1,...
                        'fontsize',handles.FontSize(2),...
                        'backgroundcolor',handles.bgcolor{3}); 
pannelsize=get(handles.uipanel_Output,'position');
pW=pannelsize(3);
pH=pannelsize(4);                       
%
set(handles.text_a1,    'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*0.5/10 pH*6/10 pW*5/10 pH*1.5/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Finished:',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'horizontalalignment','left',...
                        'visible','on');                   
set(handles.edit_a1,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*3/10 pH*6/10 pW*3/10 pH*2.0/10],...
                        'fontsize',handles.FontSize(1)+2,...
                        'string','1/10',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'foregroundcolor','r',...
                        'horizontalalignment','right',...
                        'visible','on',...
                        'enable','on');  
set(handles.text_a2,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*0.5/10 pH*2/10 pW*5/10 pH*1.5/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Calculation error:',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'horizontalalignment','left',...
                        'visible','on');                   
set(handles.edit_a2,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*3/10 pH*2/10 pW*3/10 pH*2.0/10],...
                        'fontsize',handles.FontSize(1)+2,...
                        'string',1,...
                        'backgroundcolor',handles.bgcolor{3},...
                        'foregroundcolor','r',...
                        'horizontalalignment','right',...
                        'visible','on',...
                        'enable','on'); 

%----------------------------------------
% pannel AOC                   
set(handles.uipanel_AOC,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[FigW*0.5/20 FigH*0/20 FigW*19/20 FigH*2.25/20],...
                        'Title','',...
                        'visible','on',...
                        'highlightcolor',handles.bgcolor{3},...
                        'borderwidth',1,...
                        'fontsize',handles.FontSize(2),...
                        'backgroundcolor',handles.bgcolor{3}); 
pannelsize=get(handles.uipanel_AOC,'position');                    
pW=pannelsize(3);
pH=pannelsize(4);    
set(handles.pb_Config,....
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*0.5/10 pH*2/10 pW*1.8/10 pH*6/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Configuration',...
                        'backgroundcolor',handles.bgcolor{3});
set(handles.pb_Cal,      'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*2.9/10 pH*2/10 pW*1.8/10 pH*6/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Calculate',...
                        'backgroundcolor',handles.bgcolor{3});
set(handles.pb_OK,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*5.3/10 pH*2/10 pW*1.8/10 pH*6/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','OK',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'enable','on',...
                        'visible','on'); 
set(handles.pb_Cancel,....
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*7.7/10 pH*2/10 pW*1.8/10 pH*6/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Cancel',...
                        'backgroundcolor',handles.bgcolor{3});
%
guidata(hObject, handles);  

switch CI.IsRun.GUI_TD_Cal_JLI_AMorgans
    case 1
        Fcn_Plots(hObject)
        Fcn_set_ui_enable(hObject,handles,1);
        text1 = [num2str(CI.TD.nPeriod) '/' num2str(CI.TD.nPeriod)];
        set(handles.edit_a1, 'string', text1);
        set(handles.edit_a2, 'string', num2str(CI.TD.JLI.error));
    otherwise
end

%
% -------------------------------------------------------------------------
function Fcn_main_calculation_JLI_AMorgans(varargin)
hObject     = varargin{1};
handles     = guidata(hObject);
global CI
indexFM2 = find(CI.FM.indexFM == 2);
for ss = 1:length(indexFM2)
   if CI.FM.HP{1,indexFM2(ss)}.NL.style == 3
       CI.TD.indexFMJLI = indexFM2(ss);   % the maximum number of this kind of interface is 1
   end
end
% 
%
assignin('base','CI',CI);
Fcn_TD_INI_Period_Seperation(CI.TD.nPeriod,CI.TD.RatioGapPadding)
Fcn_set_ui_enable(hObject,handles,0);
for mm = 1:CI.TD.nPeriod
    % ------------------------
    text1 = [num2str(mm) '/' num2str(CI.TD.nPeriod)];
    set(handles.edit_a1, 'string', text1);
    drawnow;
    % ------------------------
    Var(1)      = CI.TD.nPadding + (CI.TD.indexPeriod(mm,1)-1)*CI.TD.nGap + 1;  % Start of the period
    Var(2)      = CI.TD.nPadding +  CI.TD.indexPeriod(mm,2)*CI.TD.nGap;         % End of the period                                           % Length of the period
    N1          = CI.TD.numGap*CI.TD.nGap;
    if mm == CI.TD.nPeriod
        N1 = CI.TD.nTotal - Var(1) + 1;
    end
    % We only focus on the envelope from samples Var(1) to Var(1)-1+N1, 
    % the rest from Var(1) + N1 to Var(2) is the padding signal which is used
    % to avoid errors in Hilbert Transform calculation
    VarPad(1)   = CI.TD.nPadding + 1;                                           % Start of the signal to be processed by Hilbert transform
    VarPad(2)   = Var(2);                                                       % End of the signal to ...
    N           = VarPad(2) - VarPad(1) + 1;
    % ------------------------------------
    uRatioEnv2  = zeros(1,N);
    for ss = 1 : CI.TD.JLI.iteration.num
        % --------------------------
        Fcn_TD_calculation_one_period(mm)
        % --------------------------
        uRatioEnv1 = Fcn_calculation_envelope(CI.TD.uRatio(CI.TD.indexFMJLI,VarPad(1):VarPad(2)),N);       % calculate uRatio
        uRatioEnv1 = 0.5*(uRatioEnv1 + uRatioEnv2);
        [   CI.TD.Lf(Var(1):VarPad(2)),...
            CI.TD.tauf(Var(1):VarPad(2))]...
            = Fcn_flame_model_NL_JLi_AMorgans(uRatioEnv1(Var(1)-CI.TD.nPadding:end),CI.TD.indexFMJLI);       % only update the nonlinear value from Var(1) to VarPad(2)
            CI.TD.taufRem(CI.TD.indexFMJLI,Var(1):VarPad(2))  = CI.TD.tauf(CI.TD.indexFMJLI,Var(1):VarPad(2)) - CI.TD.taufMin(CI.TD.indexFMJLI);
        index = [Var(1), Var(1)+N1-1] - CI.TD.nPadding;
        error = abs(std(uRatioEnv1(index(1) : index(2)) - uRatioEnv2(index(1) : index(2)))); 
        % ---------------------
        Fcn_Plots(hObject)
        set(handles.edit_a2, 'string', num2str(error));
        drawnow
        if error < CI.TD.JLI.iteration.tol
            break;
        end
        uRatioEnv2 = uRatioEnv1;
        %----------------
    end
    % We only calculate this envelope, no others!!!
    CI.TD.uRatioEnv(Var(1):Var(2)) = uRatioEnv2((Var(1):Var(2))-CI.TD.nPadding);
end
CI.TD.JLI.error = error;
assignin('base','CI',CI);
Fcn_set_ui_enable(hObject,handles,1);

%
% -------------------------------------------------------------------------
%
function Fcn_Plots(varargin)
hObject     = varargin{1};
handles     = guidata(hObject);
hAxes1      = handles.axes1;
fontSize1   = handles.FontSize(1);
global CI
popVal      = get(handles.pop_type,'Value');
cla(hAxes1)
axes(hAxes1)
hold on
switch popVal 
    case 1
        plot(hAxes1,CI.TD.tSpTotal,CI.TD.uRatio(CI.TD.indexFMJLI,:),'-k','linewidth',0.5);
end
set(hAxes1,'YColor','k','Box','on','ygrid','on','xgrid','on');
set(hAxes1,'FontName','Helvetica','FontSize',fontSize1,'LineWidth',1)
xlabel(hAxes1,'Time [ms]','Color','k','Interpreter','LaTex','FontSize',fontSize1);
ylabel(hAxes1,'uRatio [-]','Color','k','Interpreter','LaTex','FontSize',fontSize1)
hold off


%
% -------------------------------------------------------------------------
%
function pop_type_Callback(hObject, eventdata, handles)
% hObject    handle to pop_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
% -------------------------------------------------------------------------
%
% --- Outputs from this function are returned to the command line.
function varargout = GUI_TD_Cal_JLI_AMorgans_OutputFcn(hObject, eventdata, handles) 
try
varargout{1} = handles.output;
end
%
% -------------------------------------------------------------------------
%
% --- Executes on button press in pb_Cal.
function pb_Cal_Callback(hObject, eventdata, handles)
% hObject    handle to pb_Cal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Fcn_main_calculation_JLI_AMorgans(hObject)
%
% -------------------------------------------------------------------------
%
% --- Executes on button press in pb_Cancel.
function pb_Cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pb_Cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(handles.figure);
%
% -------------------------------------------------------------------------
%
% --- Executes on button press in pb_OK.
function pb_OK_Callback(hObject, eventdata, handles)
handles         = guidata(hObject);
global CI
CI.IsRun.GUI_TD_Cal_JLI_AMorgans = 1;
assignin('base','CI',CI)
%
%
main = handles.MainGUI;                     % get the handle of OSCILOS_long
if(ishandle(main))
    mainHandles = guidata(main);            %
   set(mainHandles.TD_Plots, 'enable' , 'on') 
end
guidata(hObject, handles);
delete(handles.figure);

%
% -------------------------------------------------------------------------
%
% --- Executes during object creation, after setting all properties.
function pop_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%
% -------------------------------------------------------------------------
%
function edit_a1_Callback(hObject, eventdata, handles)
% hObject    handle to edit_a1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
% -------------------------------------------------------------------------
%
% --- Executes during object creation, after setting all properties.
function edit_a1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_a1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%
% -------------------------------------------------------------------------
%
function edit_a2_Callback(hObject, eventdata, handles)
% hObject    handle to edit_a2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function edit_a2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_a2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes when user attempts to close figure.
function figure_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(hObject);
%
% --- Executes on button press in pb_Config.
function pb_Config_Callback(hObject, eventdata, handles)
% hObject    handle to pb_Config (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global CI
prompt      = {'Termination tolerance:','Maximum loops:','Number of periods','Padding ratio'};
dlg_title   = 'Iteration configuration';
num_lines   = 1;
def         = {'1e-5','50','10','1'};
answer      = inputdlg(prompt,dlg_title,num_lines,def);
if ~isempty(answer)
    CI.TD.JLI.iteration.tol = str2double(answer{1});
    CI.TD.JLI.iteration.num = str2double(answer{2});
    CI.TD.nPeriod           = str2double(answer{3});
    CI.TD.RatioGapPadding   = str2double(answer{4});
else
    CI.TD.JLI.iteration.tol = str2double(def{1});
    CI.TD.JLI.iteration.num = str2double(def{2});
    CI.TD.nPeriod           = str2double(def{3});
    CI.TD.RatioGapPadding   = str2double(def{4});
end    
assignin('base','CI',CI);

function Fcn_set_ui_enable(hObject,handles,indexValiable)
switch indexValiable
    case 0
        set(handles.pb_Config,      'enable','off');
        set(handles.pb_OK,          'enable','off');
        set(handles.pb_Cal,         'enable','off');
        set(handles.pb_Cancel,      'enable','off');
    case 1
        set(handles.pb_Config,      'enable','on');
        set(handles.pb_OK,          'enable','on');
        set(handles.pb_Cal,         'enable','on');
        set(handles.pb_Cancel,      'enable','on');
end
% --- Executes during object creation, after setting all properties.
function figure_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --------------------------------------------------------------------
function uipushtool1_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles         = guidata(hObject);
Fig             = figure;
copyobj(handles.axes1, Fig);
set(Fig,        'units','points')
posFig          = get(handles.figure,'position');
hAxes           = get(Fig,'children');
set(hAxes(1),       'units','points',...
                    'position',[60 60 400 400],...
                    'ActivePositionProperty','position')
posAxesOuter    = [0 0 500 500];
set(Fig,        'units','points',...
                'position', [   posFig(1)+0.5*posFig(3)-0.5*posAxesOuter(3),...
                                posFig(2)+0.5*posFig(4)-0.5*posAxesOuter(4),...
                                posAxesOuter(3:4)]) 
