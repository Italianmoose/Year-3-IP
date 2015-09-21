function varargout = GUI_INI_CD(varargin)
% GUI_INI_CD MATLAB code for GUI_INI_CD.fig
%      GUI_INI_CD, by itself, creates a new GUI_INI_CD or raises the existing
%      singleton*.
%
%      H = GUI_INI_CD returns the handle to a new GUI_INI_CD or the handle to
%      the existing singleton*.
%
%      GUI_INI_CD('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_INI_CD.M with the given input arguments.
%
%      GUI_INI_CD('Property','Value',...) creates a new GUI_INI_CD or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_INI_CD_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_INI_CD_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_INI_CD

% Last Modified by GUIDE v2.5 15-Dec-2014 22:59:50

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_INI_CD_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_INI_CD_OutputFcn, ...
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
% --- Executes just before GUI_INI_CD is made visible.
function GUI_INI_CD_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_INI_CD (see VARARGIN)
indexEdit = 0;
switch indexEdit 
    case 0
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
            handles.output = hObject;
            % Initialization
            GUI_INI_CD_Initialization(hObject, eventdata, handles)
        end
        % Update handles structure
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
        guidata(hObject, handles);
        handles.output = hObject;
        guidata(hObject, handles);
        CI.IsRun.GUI_INI_CD = 0;
        assignin('base','CI',CI);                   % save the current information to the works
        GUI_INI_CD_Initialization(hObject, eventdata, handles)
%         uiwait(hObject);
end
%
% --- Executes during object creation, after setting all properties.
function figure_CreateFcn(hObject, eventdata, handles) %#ok<*DEFNU,*INUSD> 

% -------------------------------------------------------------------------
%
function GUI_INI_CD_Initialization(varargin)
hObject = varargin{1};
handles = guidata(hObject);
global CI
%--------------------------------------
% positions reconfiguration
set(0, 'units', 'points');
screenSize  = get(0, 'ScreenSize');                 % get the screen size
sW          = handles.sW;                           % screen width
sH          = handles.sH ;                          % screen height
FigW=sW.*7/10;                                        % window width
FigH=sH.*7/10;                                        % window height
set(handles.figure,     'units', 'points',...
                        'position',[(screenSize(3)-FigW)./2 (screenSize(4)-FigH)./2 FigW FigH],...
                        'name','Combustor dimensions configurations',...
                        'color',handles.bgcolor{3});
%----------------------------------------
% pannel axes
set(handles.uipanel_axes,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[FigW*0.5/20 FigH*9.25/20 FigW*19/20 FigH*10.5/20],...
                        'Title','Combustor shape preview',...
                        'visible','on',...
                        'highlightcolor',handles.bgcolor{3},...
                        'borderwidth',1,...
                        'fontsize',handles.FontSize(2),...
                        'backgroundcolor',handles.bgcolor{3});  
pannelsize = get(handles.uipanel_axes,'position');
pW = pannelsize(3);
pH = pannelsize(4);                
set(handles.axes1,      'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*0.75/10 pH*2.0/10 pW*6.25/10 pH*6/10],...
                        'fontsize',handles.FontSize(1),...
                        'color',handles.bgcolor{1});  
%----------------------------------------
%
% default combustor: Rijke tube
% default combustor dimensions: 
% Denoting the distance along the tube axis by the vector x_sample
% Radius by r_sample
handles.RT.Default = [300 700 50]; % default dimensions of the Rijke tube
%
switch CI.IsRun.GUI_INI_CD
    case 0
    CI.CD.x_sample          = [0 0.3 1.0];
    CI.CD.r_sample          = 50./1000.*[1 1 1];
    CI.CD.SectionIndex      = [0 11 0];
    CI.CD.TubeIndex         = [0 0 0];
    CI.CD.pop_CD_type = 1;
    otherwise
end 
assignin('base','CI',CI);                   % save the current information to the works
%
%----------------------------------------
% pannels combustor style
set(handles.uipanel_CB_type,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[FigW*0.5/20 FigH*7.25/20 FigW*19/20 FigH*1.75/20],...
                        'Title','',...
                        'visible','on',...
                        'highlightcolor',handles.bgcolor{3},...
                        'borderwidth',1,...
                        'fontsize',handles.FontSize(2),...
                        'backgroundcolor',handles.bgcolor{3});
pannelsize=get(handles.uipanel_CB_type,'position');
pW=pannelsize(3);
pH=pannelsize(4); 
set(handles.text_CB_type,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*0.5/10 pH*2.5/10 pW*5/10 pH*4/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Combustor type:',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'horizontalalignment','left');                         
set(handles.pop_CB_type,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*5.5/10 pH*3.5/10 pW*3.0/10 pH*4/10],...
                        'fontsize',handles.FontSize(2),...
                        'string',{  'Rijke tube';...
                                    'Load from txt file'},...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','left',...
                        'enable','on',...
                        'value',1);  
%----------------------------------------
% pannels Rijke tube configuration
set(handles.uipanel_Rijke_DM,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[FigW*0.5/20 FigH*1.75/20 FigW*19.0/20 FigH*5.25/20],...
                        'Title','',...
                        'visible','on',...
                        'highlightcolor',handles.bgcolor{3},...
                        'borderwidth',1,...
                        'fontsize',handles.FontSize(2),...
                        'backgroundcolor',handles.bgcolor{3}); 
pannelsize=get(handles.uipanel_Rijke_DM,'position');
pW=pannelsize(3);
pH=pannelsize(4);   
set(handles.text_US,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*0.5/10 pH*7.5/10 pW*4/10 pH*1.5/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Upstream length [mm]:',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'horizontalalignment','left');
set(handles.text_DS,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*0.5/10 pH*5.25/10 pW*4/10 pH*1.5/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Downstream length [mm]:',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'horizontalalignment','left');
set(handles.text_Diameter,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*0.5/10 pH*3/10 pW*4/10 pH*1.5/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Radius [mm]:',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'horizontalalignment','left');
                    
set(handles.ed_US,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*5.5/10 pH*7.5/10 pW*2/10 pH*1.75/10],...
                        'fontsize',handles.FontSize(2),...
                        'string',handles.RT.Default(1),...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','right');
set(handles.ed_DS,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*5.5/10 pH*5.25/10 pW*2/10 pH*1.75/10],...
                        'fontsize',handles.FontSize(2),...
                        'string',handles.RT.Default(2),...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','right');
set(handles.ed_Diameter,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*5.5/10 pH*3/10 pW*2/10 pH*1.75/10],...
                        'fontsize',handles.FontSize(2),...
                        'string',handles.RT.Default(3),...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','right');
set(handles.checkbox1,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*0.5/10 pH*0.75/10 pW*4/10 pH*1.5/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','With mean heat addition',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'horizontalalignment','right',...
                        'value',1);
set(handles.checkbox2,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*5.5/10 pH*0.75/10 pW*4/10 pH*1.5/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','With heat perturbations',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'horizontalalignment','right',...
                        'value',1);
%----------------------------------------
% pannels Load data                    
set(handles.uipanel_Load,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[FigW*0.5/20 FigH*1.75/20 FigW*19/20 FigH*5.25/20],...
                        'Title',get(handles.uipanel_Rijke_DM,'Title'),...
                        'visible','off',...
                        'highlightcolor',handles.bgcolor{3},...
                        'borderwidth',1,...
                        'fontsize',handles.FontSize(2),...
                        'backgroundcolor',handles.bgcolor{3});                   
pannelsize=get(handles.uipanel_Load,'position');
pW=pannelsize(3);
pH=pannelsize(4);  
currentFolder = pwd;
currentFolder = fullfile(currentFolder,'subProgram');
set(handles.edit_dir,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*0.5/10 pH*7.5/10 pW*6/10 pH*1.75/10],...
                        'fontsize',handles.FontSize(2),...
                        'string',currentFolder,...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','left');
set(handles.pb_load,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*7/10 pH*7.375/10 pW*2.5/10 pH*2.0/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Load',...
                        'backgroundcolor',handles.bgcolor{3});
%--------------
rnames = {'x [mm]','r [mm]','Section index', 'Tube index'};
% section index:
% 0: no flame, heat addition or other damping devices
% 1: flame or heat addition
% 2: damping devices
% ... for future development
table_data_num      = cat(1,1000*CI.CD.x_sample,1000*CI.CD.r_sample,CI.CD.SectionIndex,CI.CD.TubeIndex);
table_data_cell     = num2cell(table_data_num);
%
set(handles.uitable_DM,...
                        'units',        'points',...
                        'Fontunits',    'points',...
                        'position',     [pW*0.5/10 pH*1/10 pW*9/10 pH*6/10],...
                        'fontsize',     handles.FontSize(2),...
                        'rowName',      rnames,...
                        'data',         table_data_cell,...
                        'RearrangeableColumns','on');   
%----------------------------------------
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
set(handles.pb_Plot,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*1/10 pH*1.5/10 pW*2.0/10 pH*6/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Plot',...
                        'backgroundcolor',handles.bgcolor{3});
set(handles.pb_SaveFig,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*4/10 pH*1.5/10 pW*2.0/10 pH*6/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Save figure',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'visible','off');
set(handles.pb_OK,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*4/10 pH*1.5/10 pW*2.0/10 pH*6/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','OK',...
                        'backgroundcolor',handles.bgcolor{3});
set(handles.pb_Cancel,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*7/10 pH*1.5/10 pW*2.0/10 pH*6/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Cancel',...
                        'backgroundcolor',handles.bgcolor{3}); 
guidata(hObject, handles);
% switch CI.IsRun.GUI_INI_CD
%     case 1
%         if CI.IsRun.GUI_INI_PD==1 && (CI.CD.NUM_HR+CI.CD.NUM_Liner) >0
%             % Enable all the geometry settings
%             set(handles.pop_CB_type,        'enable','off'); 
%             set(handles.ed_US,              'enable','off'); 
%             set(handles.ed_DS,              'enable','off'); 
%             set(handles.ed_Diameter,        'enable','off'); 
%             set(handles.checkbox1,          'enable','off'); 
%             set(handles.checkbox2,          'enable','off'); 
%             set(handles.pb_load,            'enable','off'); 
%             set(handles.pb_Plot,            'enable','off');
%             set(handles.pb_SaveFig,         'enable','off');
%             set(handles.pb_OK,              'enable','off');
%             set(handles.pb_Cancel,          'enable','off');
%             switch CI.CD.pop_CD_type
%                 case 1
%                 set(handles.uipanel_Rijke_DM,   'visible','on'); 
%                 set(handles.uipanel_Load,       'visible','off');    
%                 case 2
%                 set(handles.uipanel_Rijke_DM,   'visible','off'); 
%                 set(handles.uipanel_Load,       'visible','on'); 
%                 set(handles.uitable_DM,         'data',table_data_cell);
%                 set(handles.pop_CB_type,        'value',2); 
%             end
%         else
%           switch CI.CD.pop_CD_type
%             case 1
%                 set(handles.uipanel_Rijke_DM,   'visible','on'); 
%                 set(handles.uipanel_Load,       'visible','off'); 
%                 set(handles.ed_US,              'string', num2str((CI.CD.x_sample(2)-CI.CD.x_sample(1)).*1000));
%                 set(handles.ed_DS,              'string', num2str((CI.CD.x_sample(3)-CI.CD.x_sample(2)).*1000));
%                 set(handles.ed_Diameter,        'string', num2str(1000*CI.CD.r_sample(1)));
%                 set(handles.pop_CB_type,        'value',1);
%                   if isempty(CI.CD.indexHP) == 0
%                     set(handles.checkbox1,'value', 1);
%                     set(handles.checkbox2,'value', 1);
%                   end
%                   if isempty(CI.CD.indexHP) == 1 && isempty(CI.CD.indexHA) == 0
%                     set(handles.checkbox1,'value', 1);
%                     set(handles.checkbox2,'value', 0);
%                   end
%                   if isempty(CI.CD.indexHP) == 1 && isempty(CI.CD.indexHA) == 1
%                     set(handles.checkbox1,'value', 0);
%                     set(handles.checkbox2,'value', 0);
%                   end         
%             case 2
%                 set(handles.uipanel_Rijke_DM,   'visible','off'); 
%                 set(handles.uipanel_Load,       'visible','on'); 
%                 set(handles.uitable_DM,         'data',table_data_cell);
%                 set(handles.pop_CB_type,        'value',2);
%           end
%         end
%     otherwise
% end



switch CI.IsRun.GUI_INI_CD
    case 1
        switch CI.CD.pop_CD_type
            case 2
                set(handles.uipanel_Rijke_DM,   'visible','off'); 
                set(handles.uipanel_Load,       'visible','on'); 
                set(handles.uitable_DM,         'data',table_data_cell);
                set(handles.pop_CB_type,        'value',2);
            case 1     
                try
                    if CI.IsRun.GUI_INI_PD==1 && (CI.CD.NUM_HR+CI.CD.NUM_Liner) >0
                        indexPD = 1;
                    else 
                        indexPD = 0;
                    end
                catch
                    indexPD = 0;
                end
                if indexPD == 1
                    % Enable all the geometry settings
                    set(handles.pop_CB_type,        'enable','off'); 
                    set(handles.ed_US,              'enable','off'); 
                    set(handles.ed_DS,              'enable','off'); 
                    set(handles.ed_Diameter,        'enable','off'); 
                    set(handles.checkbox1,          'enable','off'); 
                    set(handles.checkbox2,          'enable','off'); 
                    set(handles.pb_load,            'enable','off'); 
                    set(handles.pb_Plot,            'enable','off');
                    set(handles.pb_SaveFig,         'enable','off');
                    set(handles.pb_OK,              'enable','off');
                    set(handles.pb_Cancel,          'enable','off');
                    set(handles.uipanel_Rijke_DM,   'visible','on'); 
                    set(handles.uipanel_Load,       'visible','off');  
                else
                    set(handles.uipanel_Rijke_DM,   'visible','on'); 
                    set(handles.uipanel_Load,       'visible','off'); 
                    set(handles.ed_US,              'string', num2str((CI.CD.x_sample(2)-CI.CD.x_sample(1)).*1000));
                    set(handles.ed_DS,              'string', num2str((CI.CD.x_sample(3)-CI.CD.x_sample(2)).*1000));
                    set(handles.ed_Diameter,        'string', num2str(1000*CI.CD.r_sample(1)));
                    set(handles.pop_CB_type,        'value',1);
                    if isempty(CI.CD.indexHP) == 0
                        set(handles.checkbox1,'value', 1);
                        set(handles.checkbox2,'value', 1);
                    end
                    if isempty(CI.CD.indexHP) == 1 && isempty(CI.CD.indexHA) == 0
                        set(handles.checkbox1,'value', 1);
                        set(handles.checkbox2,'value', 0);
                    end
                    if isempty(CI.CD.indexHP) == 1 && isempty(CI.CD.indexHA) == 1
                        set(handles.checkbox1,'value', 0);
                        set(handles.checkbox2,'value', 0);
                    end         
                end
        end
    otherwise
end
%----------------------------------------
% Update handles structure
guidata(hObject, handles);
Fcn_GUI_INI_CD_Plot_CD_Shape(hObject);        % draw the combustor shape
guidata(hObject, handles);
%
% -------------------------------------------------------------------------
%
% --- Outputs from this function are returned to the command line.
function varargout = GUI_INI_CD_OutputFcn(hObject, eventdata, handles) 
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


function ed_US_Callback(hObject, eventdata, handles)
% hObject    handle to ed_US (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_US as text
%        str2double(get(hObject,'String')) returns contents of ed_US as a double
datEdit = str2double(get(hObject, 'String'));
ValDefault = 300; 
if isnan(datEdit) || ~isreal(datEdit) ||datEdit < 0
    set(hObject, 'String', ValDefault);
    errordlg('Input must be a non-negative real number','Error');
    % when the input is not a number, it is set to the default value
end

% --- Executes during object creation, after setting all properties.
function ed_US_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_US (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ed_DS_Callback(hObject, eventdata, handles)
% hObject    handle to ed_DS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_DS as text
%        str2double(get(hObject,'String')) returns contents of ed_DS as a double
datEdit = str2double(get(hObject, 'String'));
ValDefault = 700; 
if isnan(datEdit) || ~isreal(datEdit) ||datEdit < 0
    set(hObject, 'String', ValDefault);
    errordlg('Input must be a non-negative real number','Error');
    % when the input is not a number, it is set to the default value
end
% --- Executes during object creation, after setting all properties.
function ed_DS_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_DS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function ed_Diameter_Callback(hObject, eventdata, handles)
% hObject    handle to ed_Diameter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_Diameter as text
%        str2double(get(hObject,'String')) returns contents of ed_Diameter as a double
datEdit = str2double(get(hObject, 'String'));
ValDefault = 50; 
if isnan(datEdit) || ~isreal(datEdit) ||datEdit < 0
    set(hObject, 'String', ValDefault);
    errordlg('Input must be a non-negative real number','Error');
    % when the input is not a number, it is set to the default value
end

% --- Executes during object creation, after setting all properties.
function ed_Diameter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_Diameter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of checkbox1
indexHA_mean = get(handles.checkbox1,'value');
indexHA_pert = get(handles.checkbox2,'value');
if indexHA_mean == 0 && indexHA_pert == 1
    set(handles.checkbox1,'value',0);
    set(handles.checkbox2,'value',0);
end
guidata(hObject, handles);
Fcn_GUI_INI_CD_Plot_CD_Shape(hObject);


% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of checkbox2
indexHA_mean = get(handles.checkbox1,'value');
indexHA_pert = get(handles.checkbox2,'value');
if indexHA_mean == 0 && indexHA_pert == 1
    set(handles.checkbox1,'value',1);
    set(handles.checkbox2,'value',1);
end
guidata(hObject, handles);
Fcn_GUI_INI_CD_Plot_CD_Shape(hObject);

% --- Executes on selection change in pop_CB_type.
function pop_CB_type_Callback(hObject, eventdata, handles)
% hObject    handle to pop_CB_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pop_CB_type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_CB_type
pop_CD_type=get(hObject,'Value');
switch pop_CD_type
    case 1
        set(handles.uipanel_Rijke_DM,   'visible','on'); 
        set(handles.uipanel_Load,       'visible','off'); 
    case 2
        set(handles.uipanel_Rijke_DM,   'visible','off'); 
        set(handles.uipanel_Load,       'visible','on');        
        string={...
        'The default txt file named as ''CD_example.txt'' is automatically';...
        'created in the directory ''subProgram''.'
        'The default data format is:';...
        'x [m]    r[m]    InterfaceIndex     ModuleIndex';...
        '0        0.05    0                0';...
        '0.2      0.01    0                1';...
        '0.3      0.03    0                0';...
        '...'
        '0.5      0.07    11               0';...
        '0.6      0.07    0                0';...
        '...'
        'where, x means axial position of each sectional interface,';...
        'r indicates the radius of each section,';...
        'InterfaceIndex represents the type of interface between modules:';...
        '0: a simple area change';...
        '10: with heat addition, but without heat perturbation';...
        '11: with heat addition and heat perturbations';...
        'ModuleIndex indicates the type of tube between';...
        'this and the following interface'
        '0: straight constant area duct';...
        '1: duct with linearly changing radius'};
        helpdlg(string,'Load combustor dimensions from an external txt file')
 end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function pop_CB_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_CB_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pb_load.
function pb_load_Callback(hObject, eventdata, handles)
% hObject    handle to pb_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
idExistTxtFile=exist('CD_example.txt');
if(idExistTxtFile==0)
   Fcn_CUI_INI_CD_creat_txt_file % creat one txt file including the combustor shape
end
%---------------------------------------
global CI
currentFolder = pwd;
currentFolder = fullfile(currentFolder,CI.SD.name_program);
[filename, pathname] = uigetfile({  '*.txt'}, ...
                                    'Pick a file',...
                                     currentFolder);
%---------------------------------------
if filename~=0
    addpath(pathname)               % add directory to search path
    fid=fopen(filename);
%     frewind(fid);                                     % index the first row
    C_title         = textscan(fid, '%s', 4);           % read title
    C_cell          = textscan(fid, '%f %f %f %f');        % read numeric data
    fclose(fid);
    data_num(1,:)   = C_cell{1}.*1000;
    data_num(2,:)   = C_cell{2}.*1000;
    data_num(3,:)   = C_cell{3};
    data_num(4,:)   = C_cell{4};
    data_cell       = num2cell(data_num);
    set(handles.uitable_DM,'data',data_cell);         % Update the table
    CI.CD.x_sample      = C_cell{1};                    
    CI.CD.r_sample      = C_cell{2};               
    CI.CD.SectionIndex  = C_cell{3};      
    CI.CD.TubeIndex     = C_cell{4};     
end
assignin('base','CI',CI)
%
CI.CD.indexGC   = find(CI.CD.TubeIndex == 1);       % locate the gradually area change section
CI.CD.isNoGC    = isempty(CI.CD.indexGC);
assignin('base','CI',CI)
Fcn_CD_Load_Processing
assignin('base','CI',CI)
% -----------------------
clear data_num
data_num(1,:)   = CI.CD.x_sample'.*1e3;
data_num(2,:)   = CI.CD.r_sample'.*1e3;
data_num(3,:)   = CI.CD.SectionIndex';
data_num(4,:)   = CI.CD.TubeIndex';
set(handles.uitable_DM,'data',num2cell(data_num));         % Update the table
% -----------------------
% Update handles structure
guidata(hObject, handles);
 
% --- Executes on button press in pb_OK.
function pb_OK_Callback(hObject, eventdata, handles)
% hObject    handle to pb_OK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Fcn_GUI_INI_CD_Update_Data(hObject, eventdata, handles);
global CI
handles = guidata(hObject);
guidata(hObject, handles);
Fcn_GUI_INI_CD_Plot_CD_Shape(hObject);        % draw the combustor shap

CI.IsRun.GUI_INI_CD = 1;
Fcn_GUI_INI_CD_Update_Data(hObject);
Fcn_GUI_INI_CD_Update_Main_GUI(hObject);
Fcn_Interface_location
assignin('base','CI',CI); 
delete(handles.figure);

% --- Executes on button press in pb_Plot.
function pb_Plot_Callback(hObject, eventdata, handles)
% hObject    handle to pb_Plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Fcn_GUI_INI_CD_Update_Data(hObject, eventdata, handles);
handles = guidata(hObject);
guidata(hObject, handles);
Fcn_GUI_INI_CD_Plot_CD_Shape(hObject);        % draw the combustor shape   
% uiresume(handles.figure);


% --- Executes on button press in pb_SaveFig.
function pb_SaveFig_Callback(hObject, eventdata, handles)
% hObject    handle to pb_SaveFig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
hFontsize1  = handles.FontSize(1);
Fig     = figure;
set(Fig,'units','points')
posFig = get(handles.figure,'position');
copyobj(handles.axes1, Fig);
hAxes = get(Fig,'children');
set(hAxes,      'units','points',...
                'position',[60 60 400 150])
posAxes = get(hAxes,'position');
posAxesOuter = [0 0 800 400];
set(Fig,        'units','points',...
                'position', [posFig(1)+0.5*posFig(3)-0.5*posAxesOuter(3),...
                             posFig(2)+0.5*posFig(4)-0.5*posAxesOuter(4),...
                             posAxesOuter(3:4)])
% -------------legend ---------------
newline = char(10);
legend1 = ['with mean heat addition', newline, 'and heat perturbations'];
legend2 = ['with mean heat addition', newline, 'but no heat perturbations'];
hlegend = legend(hAxes,...
                        'inlet','outlet',...
                        legend1,...
                        legend2);
set(hlegend,'fontsize',hFontsize1,'location','northeastoutside');
set(hAxes,'position',posAxes)
% hl=legend('show');
% set(hl,'interpreter','latex','Fontsize',handles.FontSize(2),'box','off','Unit','points')

% --- Executes on button press in pb_Cancel.
function pb_Cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pb_Cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(handles.figure);

% --- plot the combustor shape
function Fcn_GUI_INI_CD_Plot_CD_Shape(varargin)
hObject = varargin{1};
Fcn_GUI_INI_CD_Update_Data(hObject);
handles = guidata(hObject);
guidata(hObject, handles);
Fcn_CD_plot_with_dampers(handles.axes1,handles,1)
%
% 
% %%% Yang Dong's change, does not work! 
% %%% Start
% % --- Update the data 
% function Fcn_GUI_INI_CD_Update_Data(varargin)
% hObject = varargin{1};
% handles = guidata(hObject);
% global CI
% CI.CD.pop_CD_type=get(handles.pop_CB_type,'Value');
% switch CI.IsRun.GUI_INI_CD
%     case 0
%     switch CI.CD.pop_CD_type    
%        case 1
%         l(1)    = 0;
%         l(2)    = str2num(get(handles.ed_US,'string'));
%         l(3)    = str2num(get(handles.ed_DS,'string'));
%         r       = str2num(get(handles.ed_Diameter,'string'))./1000;
%         indexHA_mean = get(handles.checkbox1,'value');
%         indexHA_pert = get(handles.checkbox2,'value');
%         CI.CD.x_sample = cumsum(l)./1000;
%         CI.CD.r_sample = 0*CI.CD.x_sample + r;
%         if indexHA_mean == 1 && indexHA_pert == 1
%             CI.CD.SectionIndex = [0 11 0];             % with heat addition and heat perturbations
%         elseif indexHA_mean == 1 && indexHA_pert == 0
%             CI.CD.SectionIndex = [0 10 0];             % with heat addition but heat perturbations
%         elseif indexHA_mean == 0 && indexHA_pert == 0
%             CI.CD.SectionIndex = [0 0 0];              % without heat addition
%         end
%         CI.CD.TubeIndex = [0 0 0]; 
%        case 2
%         data_cell = get(handles.uitable_DM,'data');
%         CI.CD.x_sample          = cell2mat(data_cell(1,:))./1000;
%         CI.CD.r_sample          = cell2mat(data_cell(2,:))./1000;
%         CI.CD.SectionIndex      = cell2mat(data_cell(3,:));
%         CI.CD.TubeIndex         = cell2mat(data_cell(4,:));
%        otherwise
%         % Code for when there is no match.
%     end
% otherwise
% end 
% 
% %Compute the length downstream of a flame. If a flame is at the end of a section, use the min to set that to 0
% %index_flame =   find(CI.CD.SectionIndex==11);
% %CI.CD.dowst_of_heat_lengths = CI.CD.x_sample(min(index_flame + 1,end)) - CI.CD.x_sample(index_flame);
% %CI.CD.index_flame = index_flame; % This is needed for other functions
% 
% %
% assignin('base','CI',CI);                   % save the current information to the workspace
% %%% end


% --- Update the data 
function Fcn_GUI_INI_CD_Update_Data(varargin)
hObject = varargin{1};
handles = guidata(hObject);
global CI
CI.CD.pop_CD_type=get(handles.pop_CB_type,'Value');
switch CI.CD.pop_CD_type    
    case 1
        l(1)    = 0;
        l(2)    = str2num(get(handles.ed_US,'string'));
        l(3)    = str2num(get(handles.ed_DS,'string'));
        r       = str2num(get(handles.ed_Diameter,'string'))./1000;
        indexHA_mean = get(handles.checkbox1,'value');
        indexHA_pert = get(handles.checkbox2,'value');
        CI.CD.x_sample = cumsum(l)./1000;
        CI.CD.r_sample = 0*CI.CD.x_sample + r;
        if indexHA_mean == 1 && indexHA_pert == 1
            CI.CD.SectionIndex = [0 11 0];             % with heat addition and heat perturbations
        elseif indexHA_mean == 1 && indexHA_pert == 0
            CI.CD.SectionIndex = [0 10 0];             % with heat addition but heat perturbations
        elseif indexHA_mean == 0 && indexHA_pert == 0
            CI.CD.SectionIndex = [0 0 0];              % without heat addition
        end
        CI.CD.TubeIndex = [0 0 0]; 
    case 2
        data_cell = get(handles.uitable_DM,'data');
        CI.CD.x_sample          = cell2mat(data_cell(1,:))./1000;
        CI.CD.r_sample          = cell2mat(data_cell(2,:))./1000;
        CI.CD.SectionIndex      = cell2mat(data_cell(3,:));
        CI.CD.TubeIndex         = cell2mat(data_cell(4,:));
    otherwise
        % Code for when there is no match.
end
CI.CD.pop_CD_type = get(handles.pop_CB_type,'Value');

%Compute the length downstream of a flame. If a flame is at the end of a section, use the min to set that to 0
%index_flame =   find(CI.CD.SectionIndex==11);
%CI.CD.dowst_of_heat_lengths = CI.CD.x_sample(min(index_flame + 1,end)) - CI.CD.x_sample(index_flame);
%CI.CD.index_flame = index_flame; % This is needed for other functions

%
assignin('base','CI',CI);                   % save the current information to the workspace


function Fcn_GUI_INI_CD_Update_Main_GUI(varargin)
hObject = varargin{1};
handles = guidata(hObject);
main = handles.MainGUI;
global CI
% Obtain handles using GUIDATA with the caller's handle 
if(ishandle(main))
    mainHandles = guidata(main);
    %%% add a questdlg:
    selection = questdlg('Do you want to use passsive dampers?',...
      'Passive damper',...
      'Yes','No','No'); 
    switch selection 
        case 'Yes'
            changeMain = mainHandles.INI_PD;
            set(changeMain, 'Enable', 'on');
            set(changeMain, 'Visible', 'on');
            changeMain = mainHandles.INI_CD;
            set(changeMain, 'Enable', 'off');
            set(changeMain, 'Visible', 'off');
        case 'No'
            changeMain = mainHandles.INI_TP;
            set(changeMain, 'Enable', 'on');
            changeMain = mainHandles.INI_PD;
            set(changeMain, 'Visible', 'off');
    end
    String_Listbox=get(mainHandles.listbox_Info,'string');
    ind=find(ismember(String_Listbox,'<HTML><FONT color="blue">Information 1:'));
    nLength=size(String_Listbox);
    if isempty(ind)
        indStart=nLength(1);
    else
        indStart=ind-1;
        for i=nLength(1):-1:indStart+1
            String_Listbox(i)=[];
        end
    end
    String_Listbox{indStart+1}=['<HTML><FONT color="blue">Information 1:'];
    String_Listbox{indStart+2}=['<HTML><FONT color="blue">Combustor dimension:'];
    String_Listbox{indStart+3}=['The combustor has ' num2str(length(CI.CD.r_sample)-1) ' sections'];
%     String_Listbox{indStart+4}=['The flame is located after section ' num2str(CI.CD.SectionIndex_flame-1)];
% try
    String_Listbox{indStart+4}=['The x-coordinates are: '];
    String_Listbox{indStart+5}=[num2str(CI.CD.x_sample) ' m'];
    String_Listbox{indStart+6}=['The r-coordinates are: '];
    String_Listbox{indStart+7}=[num2str(1000.*CI.CD.r_sample) ' mm'];
% catch 
% end
    set(mainHandles.listbox_Info,'string',String_Listbox,'value',1);
end
guidata(hObject, handles);
%
%-----------------------------------end------------------------------------


% --------------------------------------------------------------------
function uipushtool1_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
hFontsize1  = handles.FontSize(1);
Fig     = figure;
set(Fig,'units','points')
posFig = get(handles.figure,'position');
copyobj(handles.axes1, Fig);
hAxes = get(Fig,'children');
set(hAxes,      'units','points',...
                'position',[60 60 400 150])
posAxes = get(hAxes,'position');
posAxesOuter = [0 0 800 400];
set(Fig,        'units','points',...
                'position', [posFig(1)+0.5*posFig(3)-0.5*posAxesOuter(3),...
                             posFig(2)+0.5*posFig(4)-0.5*posAxesOuter(4),...
                             posAxesOuter(3:4)])
% -------------legend ---------------
newline = char(10);
legend1 = ['with mean heat addition', newline, 'and heat perturbations'];
legend2 = ['with mean heat addition', newline, 'but no heat perturbations'];
hlegend = legend(hAxes,...
                        'inlet','outlet',...
                        legend1,...
                        legend2);
set(hlegend,'fontsize',hFontsize1,'location','northeastoutside');
set(hAxes,'position',posAxes)

 
