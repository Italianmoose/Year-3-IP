function varargout = GUI_INI_CD_Split(varargin)
% GUI_INI_CD_Split MATLAB code for GUI_INI_CD_Split.fig
%      GUI_INI_CD_Split, by itself, creates a new GUI_INI_CD_Split or raises the existing
%      singleton*.
%
%      H = GUI_INI_CD_Split returns the handle to a new GUI_INI_CD_Split or the handle to
%      the existing singleton*.
%
%      GUI_INI_CD_Split('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_INI_CD_Split.M with the given input arguments.
%
%      GUI_INI_CD_Split('Property','Value',...) creates a new GUI_INI_CD_Split or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_INI_CD_Split_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_INI_CD_Split_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_INI_CD_Split

% Last Modified by GUIDE v2.5 08-Dec-2014 09:45:49

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_INI_CD_Split_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_INI_CD_Split_OutputFcn, ...
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
% --- Executes just before GUI_INI_CD_Split is made visible.
function GUI_INI_CD_Split_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_INI_CD_Split (see VARARGIN)
handles = Fcn_GUI_default_configuration(handles);
guidata(hObject, handles);
handles.output = hObject;
guidata(hObject, handles);
GUI_Initialization(hObject, eventdata, handles)
uiwait(hObject);
% end
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
%--------------------------------------
% positions reconfiguration
set(0, 'units', 'points');
screenSize  = get(0, 'ScreenSize');                 % get the screen size
sW          = handles.sW;                           % screen width
sH          = handles.sH ;                          % screen height
FigW=sW.*1/2;                                        % window width
FigH=sH.*2/5;                                        % window height
set(handles.figure,     'units', 'points',...
                        'position',[(screenSize(3)-FigW)./2 (screenSize(4)-FigH)./2 FigW FigH],...
                        'name','Set the split number',...
                        'color',handles.bgcolor{3});

%----------------------------------------
set(handles.uipanel2,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[FigW*0.5/20 FigH*12/20 FigW*19/20 FigH*7.5/20],...
                        'Title','',...
                        'visible','on',...
                        'highlightcolor',handles.bgcolor{3},...
                        'borderwidth',1,...
                        'fontsize',handles.FontSize(2),...
                        'backgroundcolor',handles.bgcolor{3});
pannelsize=get(handles.uipanel2,'position');
pW=pannelsize(3);
pH=pannelsize(4); 
% Tstring={'Some sections are set to gradually area change sections; It is necessary to split these sections to small pieces'};
% Tstring={'Some sections are set to gradually '};
msg={   '<HTML><FONT color="blue">Some sections are set to gradually area change sections;';...
        '<HTML><FONT color="blue">It is necessary to split these sections to small pieces;';...
        '<HTML><FONT color="blue">Please set the split number for every section;';...
        '<HTML><FONT color="blue">The default number for every section is 50;';...
        ''};
set(handles.listbox1,...
                        'units', 'points',...
                        'fontunits','points',...
                        'position',[pW*0/20 pH*0/20 pW*20/20 pH*20/20],...
                        'fontsize',handles.FontSize(1),...
                        'string',msg,...
                        'backgroundcolor',handles.bgcolor{4},...
                        'value',1);  
%----------------------------------------
set(handles.uipanel1,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[FigW*0.5/20 FigH*3.5/20 FigW*19/20 FigH*8/20],...
                        'Title','',...
                        'visible','on',...
                        'highlightcolor',handles.bgcolor{3},...
                        'borderwidth',1,...
                        'fontsize',handles.FontSize(2),...
                        'backgroundcolor',handles.bgcolor{3});
pannelsize=get(handles.uipanel1,'position');
pW=pannelsize(3);
pH=pannelsize(4); 
set(handles.text1,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*0.5/10 pH*6/10 pW*4/10 pH*1.5/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','From left to right ',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'horizontalalignment','left');   
% for k = 1:length(CI.CD.indexGC)
for k = 1:2
    StringSplit{k} = ['split section number:' num2str(k)]; 
end
set(handles.pop1,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*5.5/10 pH*6/10 pW*4/10 pH*2/10],...
                        'fontsize',handles.FontSize(2),...
                        'string',StringSplit,...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','left',...
                        'enable','on',...
                        'value',1);  
set(handles.text2,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*0.5/10 pH*2/10 pW*4/10 pH*1.5/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','split number [-]:',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'horizontalalignment','left');
set(handles.edit2,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*5.5/10 pH*2/10 pW*4/10 pH*2/10],...
                        'fontsize',handles.FontSize(2),...
                        'string',50,...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','right');
%----------------------------------------
% pannel AOC                    
set(handles.uipanel_AOC,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[FigW*0.5/20 FigH*0/20 FigW*19/20 FigH*3/20],...
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
Fcn_Para_Initialization
%
% -------------------------------------------------------------------------
%
% --- Outputs from this function are returned to the command line.
function varargout = GUI_INI_CD_Split_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Get default command line output from handles structure
% varargout{1} = [];
% delete(hObject);
% try
% varargout{1} = handles.output;
% end
varargout{1} = [];
delete(hObject);

% --- Executes when user attempts to close figure.
function figure_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiresume(hObject);
% delete(hObject);


% --- Executes on selection change in pop1.
function pop1_Callback(hObject, eventdata, handles)
% hObject    handle to pop1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Fcn_UI2Para(hObject)
Fcn_Para2UI(hObject)


function Fcn_Para_Initialization
global CI
N = length(CI.CD.indexGC); % number of sections to be split
CI.CD.GC.numSplitting(1:N)  = 50;       % default splitting number 
CI.CD.GC.popupIndex         = 1;        % default popup menu number
assignin('base','CI',CI);


function Fcn_Para2UI(hObject)
handles = guidata(hObject);
popNum  = get(handles.pop1,'value');
global CI
CI.CD.GC.popupIndex = popNum;
set(handles.edit2, 'string', num2str(CI.CD.GC.numSplitting(popNum)));
guidata(hObject, handles);

function Fcn_UI2Para(hObject)
handles = guidata(hObject);
global CI
CI.CD.GC.numSplitting(CI.CD.GC.popupIndex) = str2double(get(handles.edit2,'string'));
assignin('base','CI',CI);


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

% --- Executes on button press in pb_OK.
function pb_OK_Callback(hObject, eventdata, handles)
% hObject    handle to pb_OK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Fcn_UI2Para(hObject)
% delete(handles.figure);
uiresume(handles.figure);


% --- Executes on button press in pb_Cancel.
function pb_Cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pb_Cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% delete(handles.figure);
uiresume(handles.figure);


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

