function varargout = GUI_INI_FM_Sel(varargin)
% GUI_INI_FM_Sel MATLAB code for GUI_INI_FM_Sel.fig
%      GUI_INI_FM_Sel, by itself, creates a new GUI_INI_FM_Sel or raises the existing
%      singleton*.
%
%      H = GUI_INI_FM_Sel returns the handle to a new GUI_INI_FM_Sel or the handle to
%      the existing singleton*.
%
%      GUI_INI_FM_Sel('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_INI_FM_Sel.M with the given input arguments.
%
%      GUI_INI_FM_Sel('Property','Value',...) creates a new GUI_INI_FM_Sel or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_INI_FM_Sel_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_INI_FM_Sel_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_INI_FM_Sel

% Last Modified by GUIDE v2.5 11-Dec-2014 11:50:06

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_INI_FM_Sel_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_INI_FM_Sel_OutputFcn, ...
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
% --- Executes just before GUI_INI_FM_Sel is made visible.
function GUI_INI_FM_Sel_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_INI_FM_Sel (see VARARGIN)
global CI
indexEdit = 0;
switch indexEdit 
    case 0
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
    case 1
        CI.IsRun.GUI_INI_FM_Sel = 0;
        assignin('base','CI',CI);
        handles = Fcn_GUI_default_configuration(handles);
end
guidata(hObject, handles);
handles.output = hObject;
guidata(hObject, handles);
GUI_Initialization(hObject, eventdata, handles)

%
% --- Executes during object creation, after setting all properties.
function figure_CreateFcn(hObject, eventdata, handles) %#ok<*DEFNU,*INUSD> 
%
% -------------------------------------------------------------------------
%
function GUI_Initialization(varargin)
hObject = varargin{1};
handles = guidata(hObject);
global CI
CI.FM.ModelType = {     'Linear FTF model';...
                        'Nonlinear FDF model';...
                        'Experimental/CFD fitted FDF';...
                        'Fully non-linear G-Equation model'};
assignin('base','CI',CI);
%--------------------------------------
% positions reconfiguration
set(0, 'units', 'points');
screenSize  = get(0, 'ScreenSize');                 % get the screen size
sW          = handles.sW;                           % screen width
sH          = handles.sH ;                          % screen height
FigW=sW.*1/2;                                        % window width
FigH=sH.*3/4;                                        % window height
set(handles.figure,     'units', 'points',...
                        'position',[(screenSize(3)-FigW)./2 (screenSize(4)-FigH)./2 FigW FigH],...
                        'name','Choose (flame) model for each unsteady heat source',...
                        'color',handles.bgcolor{3});

%----------------------------------------
set(handles.uipanel_Msg,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[FigW*0.5/20 FigH*12.75/20 FigW*19/20 FigH*6.75/20],...
                        'Title','',...
                        'visible','on',...
                        'highlightcolor',handles.bgcolor{3},...
                        'borderwidth',1,...
                        'fontsize',handles.FontSize(2),...
                        'backgroundcolor',handles.bgcolor{3});
pannelsize=get(handles.uipanel_Msg,'position');
pW=pannelsize(3);
pH=pannelsize(4); 
msg={   '<HTML><FONT color="blue">This GUI is used to choose (flame) model for each unsteady heat source;';...
        '<HTML><FONT color="blue">Four choices can be prescribed:';...
        '<HTML><FONT color="blue">1. A linear (flame) transfer function model;';...
        '<HTML><FONT color="blue">2. A simple nonlinear (flame) describing function model,';...
        '<HTML><FONT color="blue">which is assumed can be decoupled as a nonlinear model and a';...
        '<HTML><FONT color="blue">linear (flame) transfer function (FTF) model.';...
        '<HTML><FONT color="blue">3. Experimental/CFD (flame) transfer functions(loaded from a file) for';...
        '<HTML><FONT color="blue">different velocities, and fit the FTF data with state-space models.';...
        '<HTML><FONT color="blue">4. The fully non-linear G-Equation model (Williams 1988)'};
set(handles.listbox1,...
                        'units', 'points',...
                        'fontunits','points',...
                        'position',[pW*0/20 pH*0/20 pW*20/20 pH*20/20],...
                        'fontsize',handles.FontSize(1),...
                        'string',msg,...
                        'backgroundcolor',handles.bgcolor{4},...
                        'value',1);  
%----------------------------------------
set(handles.uipanel_pop,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[FigW*0.5/20 FigH*6.5/20 FigW*19/20 FigH*6.0/20],...
                        'Title','',...
                        'visible','on',...
                        'highlightcolor',handles.bgcolor{3},...
                        'borderwidth',1,...
                        'fontsize',handles.FontSize(2),...
                        'backgroundcolor',handles.bgcolor{3});
pannelsize=get(handles.uipanel_pop,'position');
pW=pannelsize(3);
pH=pannelsize(4); 
set(handles.text1,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*0.5/10 pH*7/10 pW*4/10 pH*1.5/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Choose heat source:',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'horizontalalignment','left');   
StringPop1  = { 'Unsteady Heat source 1'};
set(handles.pop1,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*4.5/10 pH*7/10 pW*5/10 pH*2/10],...
                        'fontsize',handles.FontSize(2),...
                        'string',StringPop1,...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','left',...
                        'enable','on',...
                        'value',1);  
set(handles.text2,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*0.5/10 pH*4/10 pW*4/10 pH*1.5/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Choose (Flame) model:',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'horizontalalignment','left');
set(handles.pop2,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*4.5/10 pH*4/10 pW*5.0/10 pH*2/10],...
                        'fontsize',handles.FontSize(2),...
                        'string',CI.FM.ModelType,...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','left',...
                        'enable','on',...
                        'value',1);
set(handles.text3,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*0.5/10 pH*1/10 pW*5/10 pH*1.5/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Edit selected flame model :',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'horizontalalignment','left');
set(handles.pb_Edit,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*7/10 pH*1/10 pW*2.5/10 pH*2/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Edit',...
                        'backgroundcolor',handles.bgcolor{3});
%----------------------------------------
set(handles.uipanel_table,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[FigW*0.5/20 FigH*2/20 FigW*19/20 FigH*4.25/20],...
                        'Title','',...
                        'visible','on',...
                        'highlightcolor',handles.bgcolor{3},...
                        'borderwidth',1,...
                        'fontsize',handles.FontSize(2),...
                        'backgroundcolor',handles.bgcolor{3});
pannelsize=get(handles.uipanel_table,'position');
pW=pannelsize(3);
pH=pannelsize(4); 

%--------------
rnames = {'Heat source 1'};
cnames = {'(Flame) model type', 'Edited?'};
table_cell(1,1)={'Linear FTF model'};
table_cell(1,2)={'N'};
%
set(handles.uitable1,...
                        'units',        'points',...
                        'Fontunits',    'points',...
                        'position',     [pW*0.5/10 pH*1/10 pW*9/10 pH*8/10],...
                        'fontsize',     handles.FontSize(2),...
                        'rowName',      rnames,...
                        'columnName',   cnames,...
                        'data',         table_cell,...
                        'RearrangeableColumns','on');  
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
set(handles.pb_OK,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*1/10 pH*1.5/10 pW*3/10 pH*6/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','OK',...
                        'backgroundcolor',handles.bgcolor{3});
set(handles.pb_Cancel,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*6/10 pH*1.5/10 pW*3/10 pH*6/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Cancel',...
                        'backgroundcolor',handles.bgcolor{3}); 
guidata(hObject, handles);
Fcn_Further_Initialization(hObject);
%
% -------------------------------------------------------------------------
%
% --- Outputs from this function are returned to the command line.
function varargout = GUI_INI_FM_Sel_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Get default command line output from handles structure
% varargout{1} = [];
% delete(hObject);
try
varargout{1} = handles.output;
end

% --- Executes when user attempts to close figure.
function figure_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% uiresume(hObject);
delete(hObject);


% --- Executes on selection change in pop1.
function pop1_Callback(hObject, eventdata, handles)
% hObject    handle to pop1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Fcn_UI2Para(hObject)
% Fcn_Para2UI(hObject)
global CI
indexPop1 = get(handles.pop1, 'value');
set(handles.pop2, 'value', CI.FM.indexFM(indexPop1));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function pop1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in pop2.
function pop2_Callback(hObject, eventdata, handles)
% hObject    handle to pop2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pop2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop2




% --- Executes during object creation, after setting all properties.
function pop2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function Fcn_Further_Initialization(varargin)
hObject = varargin{1};
handles = guidata(hObject);
global CI
switch CI.IsRun.GUI_INI_FM_Sel
    case 0
        if ~isempty(CI.CD.indexHP)
            NumHP = length(CI.CD.indexHP);
            for ss = 1:NumHP
                CI.FM.indexFM(ss)       =   1;      % defult setting: linear FTF
                CI.FM.HP{ss}         =   Fcn_initialization_flame_model(CI.FM.indexFM(ss));
                rnames{ss}           =   ['Heat source '  num2str(ss)];
                table_cell(ss,1)     =   CI.FM.ModelType(1);
                table_cell(ss,2)     =   {'N'};
                StringPop1{ss}       =   ['Unsteady heat source '  num2str(ss)];
            end
            set(handles.uitable1,       'rowName',      rnames,...
                                        'data',         table_cell); 
            set(handles.pop1,...
                                        'string',StringPop1,...
                                        'value',1); 
            set(handles.pb_OK, 'enable','off');
        end
    case 1
        if ~isempty(CI.CD.indexHP)
            NumHP = length(CI.CD.indexHP);
            for ss = 1:NumHP
               rnames{ss}           =   ['Heat source '  num2str(ss)];
               table_cell(ss,1)     =   CI.FM.ModelType(CI.FM.indexFM(ss));
               HP = CI.FM.HP{ss};
               if HP.IsRun == 0
                    table_cell(ss,2) = {'N'};
               else
                    table_cell(ss,2) = {'Y'};
               end
               StringPop1{ss}       =   ['Unsteady heat source '  num2str(ss)];
            end
            set(handles.uitable1,       'rowName',      rnames,...
                                        'data',         table_cell); 
            set(handles.pop1,...
                                        'string',StringPop1,...
                                        'value',1);
            set(handles.pop2,           'value',CI.FM.indexFM(1))
            set(handles.pb_OK,          'enable','on');
        end
end
assignin('base','CI',CI);
guidata(hObject, handles);


function HP = Fcn_initialization_flame_model(indexFM)
% This function is used to initialize the flame model for a selected
% unsteady heat source
% HP means heat perturbations
HP.IsRun = 0;    % an index to show if this flame model has ever been edited
switch indexFM
    case 1
        HP.FTF.style     = 3;
        HP.NL.style      = 1;    % when this value ==1, no nonlinear model
        HP.FTF.af        = 1;
        HP.FTF.tauf      = 3e-3;
        HP.FTF.fc        = 75;
        HP.FTF.omegac    = 2*pi*HP.FTF.fc;
        HP.FTF.num       = HP.FTF.omegac;
        HP.FTF.den       = [1 HP.FTF.omegac];
        % --------------------------
        HP.FTF.xi        = 1;    % damping coefficient in the second order filter model
        % --------------------------
    case 2
        HP.FTF.style     = 3;
        HP.NL.style      = 3;
        HP.FTF.af        = 1;
        HP.FTF.tauf      = 3e-3;
        HP.FTF.fc        = 75;
        HP.FTF.omegac    = 2*pi*HP.FTF.fc;
        HP.FTF.num       = HP.FTF.omegac;
        HP.FTF.den       = [1 HP.FTF.omegac];
        % --------------------------
        HP.FTF.xi        = 1;
        % --------------------------
        % J.Li and A.Morgans's model
        HP.NL.Model3.alpha       = 0.85;
        HP.NL.Model3.beta        = 50;
        HP.NL.Model3.taufN       = 0e-3;
        % -------------------------
        % Dowling's Model
        HP.NL.Model2.alpha       = 0.3;
    case 3
        HP.FMEXP.indexIMPORT    = 1;
        HP.FMEXP.nFTF           = 0;
        HP.FMEXP.indexModify    = 0;
    case 4
        
end
%
% --- Executes on button press in pb_Edit.
function pb_Edit_Callback(hObject, eventdata, handles)
% hObject    handle to pb_Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global CI
HP_num  = get(handles.pop1, 'value');
HP      = CI.FM.HP{HP_num};
indexFM = get(handles.pop2,'value');
switch HP.IsRun
    case 0   % never run before
        CI.FM.HP{HP_num} = Fcn_initialization_flame_model(indexFM);
    case 1   % has ever been run
        if indexFM ~= CI.FM.indexFM(HP_num)
            string={'The new flame model type is different from before. Do you want to continue?'};
            selection = questdlg(string,'Different flame model',...
              'Yes','No','Yes'); 
            switch selection, 
              case 'Yes'
                  CI.FM.HP{HP_num} = Fcn_initialization_flame_model(indexFM);
              case 'No'
              return 
            end
        else
        end
end
switch indexFM
    case {1,2}
        GUI_INI_FM(handles.figure, HP_num, indexFM);
    case 3
        GUI_INI_FMEXP(handles.figure, HP_num, indexFM);
    case 4 % G-equation case
        GUI_INI_GEQU(handles.figure, HP_num, indexFM);
end
assignin('base','CI',CI);
guidata(hObject, handles);


function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double
datEdit = str2double(get(hObject, 'String'));
ValDefault = 50; 
if isnan(datEdit) || ~isreal(datEdit) ||datEdit < 0
    set(hObject, 'String', ValDefault);
    errordlg('Input must be a positive integer','Error');
    % when the input is not a number, it is set to the default value
end

% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pb_OK.
function pb_OK_Callback(hObject, eventdata, handles)
% hObject    handle to pb_OK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global CI
CI.IsRun.GUI_INI_FM_Sel = 1;
assignin('base','CI',CI);
Fcn_GUI_Update_Data(hObject, eventdata, handles)
delete(handles.figure);
% uiresume(handles.figure);

% --- Update the data when clicking 'OK' or 'Apply'
function Fcn_GUI_Update_Data(hObject, eventdata, handles)
handles = guidata(hObject);
global CI
main = handles.MainGUI;
% Obtain handles using GUIDATA with the caller's handle
if(ishandle(main))
    mainHandles = guidata(main);
    changeMain = mainHandles.INI_BC;
    set(changeMain, 'Enable', 'on');
    String_Listbox=get(mainHandles.listbox_Info,'string');
    ind=find(ismember(String_Listbox,'<HTML><FONT color="blue">Information 4:'));
    nLength=size(String_Listbox);
    if isempty(ind)
        indStart=nLength(1);
    else
    indStart=ind-1;
        for i=nLength(1):-1:indStart+1
            String_Listbox(i)=[];
        end
    end
    String_Listbox{indStart+1}=['<HTML><FONT color="blue">Information 4:'];
    String_Listbox{indStart+2}=['<HTML><FONT color="blue"> (Flame) models:'];
    String_Listbox{indStart+3}=['The unsteady heat source(s) is(are) located at:'];
    String_Listbox{indStart+4}=[num2str(CI.CD.x_sample(CI.CD.indexHP)) ' m' ];
    String_Listbox{indStart+5}=['The types of selected (flame) model(s) is(are):'];
    for ss = 1:length(CI.FM.indexFM)
        String_Listbox{indStart+5+ss}=['(' num2str(ss) ') ' CI.FM.ModelType{CI.FM.indexFM(ss)}];
    end
    set(mainHandles.listbox_Info,'string',String_Listbox);
end
guidata(hObject, handles);
% guidata(hObject, handles);
assignin('base','CI',CI); % save the current information to the works


% --- Executes on button press in pb_Cancel.
function pb_Cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pb_Cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(handles.figure);
% uiresume(handles.figure);


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
