function varargout = GUI_INI_BC_Entropy(varargin)
% GUI_INI_BC_ENTROPY MATLAB code for GUI_INI_BC_Entropy.fig
%      GUI_INI_BC_ENTROPY, by itself, creates a new GUI_INI_BC_ENTROPY or raises the existing
%      singleton*.
%
%      H = GUI_INI_BC_ENTROPY returns the handle to a new GUI_INI_BC_ENTROPY or the handle to
%      the existing singleton*.
%
%      GUI_INI_BC_ENTROPY('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_INI_BC_ENTROPY.M with the given input arguments.
%
%      GUI_INI_BC_ENTROPY('Property','Value',...) creates a new GUI_INI_BC_ENTROPY or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_INI_BC_Entropy_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_INI_BC_Entropy_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_INI_BC_Entropy

% Last Modified by GUIDE v2.5 16-Dec-2014 08:28:53

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_INI_BC_Entropy_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_INI_BC_Entropy_OutputFcn, ...
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


% --- Executes just before GUI_INI_BC_Entropy is made visible.
function GUI_INI_BC_Entropy_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_INI_BC_Entropy (see VARARGIN)
indexEdit = 0;
switch indexEdit 
    case 0
        %--------------------------------------------------------------------------
        dontOpen = false;
        mainGuiInput = find(strcmp(varargin, 'GUI_INI_BC'));
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
            GUI_INI_BC_Entropy_Initialization(hObject, eventdata, handles)
        end
        guidata(hObject, handles);
        handles.output = hObject;
        guidata(hObject, handles);
        if dontOpen
           disp('-----------------------------------------------------');
           disp('This is a subprogram. It cannot be run independently.') 
           disp('Please load the program "GUI_INI_BC'' from the ')
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
        guidata(hObject, handles);  
        GUI_INI_BC_Entropy_Initialization(hObject, eventdata, handles)
        handles.output = hObject;
        guidata(hObject, handles);
end



function GUI_INI_BC_Entropy_Initialization(varargin)
hObject = varargin{1};
handles = guidata(hObject);
global CI
assignin('base','CI',CI);                   % save the current information to the workspace            
% guidata(hObject, handles);
% ------------------------------------------
% positions reconfiguration
set(0, 'units', 'points');
screenSize  = get(0, 'ScreenSize');                     % get the screen size
sW          = handles.sW;                               % screen width
sH          = handles.sH;                               % screen height
FigW=sW*3/5;                                                % window width
FigH=sH.*1/2;                                           % window height
set(handles.figure,     'units', 'points',...
                        'position',[(screenSize(3)-FigW)./2 (screenSize(4)-FigH)./2 FigW FigH],...
                        'name','Entropy advection model configuration',...
                        'color',handles.bgcolor{3});
%----------------------------------------
% pannel axes
set(handles.uipanel_axes,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[FigW*9.5/20 FigH*2.5/20 FigW*10/20 FigH*17.25/20],...
                        'Title','Plots',...
                        'visible','on',...
                        'highlightcolor',handles.bgcolor{3},...
                        'borderwidth',1,...
                        'fontsize',handles.FontSize(2),...
                        'backgroundcolor',handles.bgcolor{3});  
pannelsize=get(handles.uipanel_axes,'position');
pW=pannelsize(3);
pH=pannelsize(4);    
set(handles.pop_plot,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*0.5/10 pH*8.25/10 pW*7/10 pH*0.8/10],...
                        'fontsize',handles.FontSize(2),...
                        'string',{'Time domain results';...
                        'Frequency domain results'},...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','left',...
                        'Enable','on');  
set(handles.axes1,      'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*2.5/10 pH*2.0/10 pW*7/10 pH*5.5/10],...
                        'fontsize',handles.FontSize(1),...
                        'color',handles.bgcolor{1},...
                        'box','on');  
guidata(hObject, handles);
%----------------------------------------
% pannel dissipation
set(handles.uipanel_Dissipation,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[FigW*0.5/20 FigH*13.25/20 FigW*8.5/20 FigH*6.5/20],...
                        'Title','Dissipation model',...
                        'visible','on',...
                        'highlightcolor',handles.bgcolor{3},...
                        'borderwidth',1,...
                        'fontsize',handles.FontSize(2),...
                        'backgroundcolor',handles.bgcolor{3}); 
pannelsize = get(handles.uipanel_Dissipation,'position');
pW = pannelsize(3);
pH = pannelsize(4);   
%
set(handles.text_dissipation,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*0.5/10 pH*2.0/10 pW*3/10 pH*3.5/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Dissipation ratio [-]:',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'horizontalalignment','left',...
                        'visible','on');                   
set(handles.edit_dissipation,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*6/10 pH*3.5/10 pW*3.5/10 pH*2.0/10],...
                        'fontsize',handles.FontSize(2),...
                        'string',1,...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','right',...
                        'visible','on');  
%----------------------------------------
% pannels Dispersion
set(handles.uipanel_Dispersion,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[FigW*0.5/20 FigH*2.5/20 FigW*8.5/20 FigH*10.5/20],...
                        'Title','Dispersion model',...
                        'visible','on',...
                        'highlightcolor',handles.bgcolor{3},...
                        'borderwidth',1,...
                        'fontsize',handles.FontSize(2),...
                        'backgroundcolor',handles.bgcolor{3});  
pannelsize=get(handles.uipanel_Dispersion,'position');
pW=pannelsize(3);
pH=pannelsize(4);   
set(handles.text_Disp_model,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*0.5/10 pH*6.5/10 pW*5/10 pH*1.2/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Type:',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'horizontalalignment','left',...
                        'visible','off');  
set(handles.pop_Disp_model,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*0.5/10 pH*7/10 pW*9/10 pH*0.8/10],...
                        'fontsize',handles.FontSize(2),...
                        'string',{'No dispersion';...
                        'Morgans et al.(2013)';...
                        'Sattelmayer (2003)'},...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','left',...
                        'Enable','on');  
set(handles.text_Delta_tauCs,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*0.5/10 pH*1/10 pW*5/10 pH*4/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Dispersion time estimation ---- Delta tauCs [ms]:',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'horizontalalignment','left',...
                        'visible','on');                   
set(handles.edit_Delta_tauCs,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*6/10 pH*3.75/10 pW*3.5/10 pH*1.25/10],...
                        'fontsize',handles.FontSize(2),...
                        'string',0,...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','right',...
                        'visible','on',...
                        'enable','off');  
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
set(handles.pb_Apply,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*1/10 pH*2/10 pW*2.0/10 pH*6/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Plot',...
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
guidata(hObject, handles);
%
global ET
switch CI.IsRun.GUI_INI_BC_Entropy 
    case 0
        ET.pop_type_model         = 2;
        ET.Dissipation.k          = 1;
        ET.Dispersion.Delta_tauCs = 3e-3;
        set(handles.pop_Disp_model, 'value',ET.pop_type_model);
        set(handles.edit_dissipation,'string',num2str(ET.Dissipation.k));
        set(handles.edit_Delta_tauCs,'string',num2str(ET.Dispersion.Delta_tauCs.*1e3));
    case 1
        set(handles.pop_Disp_model, 'value',CI.BC.ET.pop_type_model);
        set(handles.edit_dissipation,'string',num2str(CI.BC.ET.Dissipation.k));
        set(handles.edit_Delta_tauCs,'string',num2str(CI.BC.ET.Dispersion.Delta_tauCs.*1e3));
end
%
assignin('base','ET',ET);                   % save the current information to the workspace
Fcn_GUI_INI_BC_Entropy_update_plot(hObject)
guidata(hObject, handles);

% --- plot the dispersion of  resident time of entropy
function Fcn_GUI_INI_BC_Entropy_update_plot(varargin)
hObject     = varargin{1};
handles     = guidata(hObject);
global ET
hAxes1      = handles.axes1;
fontSize1   = handles.FontSize(1);
fontSize2   = handles.FontSize(2);
% -------------
ET.pop_type_model               = get(handles.pop_Disp_model,'Value');
ET.Dispersion.Delta_tauCs       = 1e-3*str2num(get(handles.edit_Delta_tauCs,'String'));
ET.Dissipation.k                = str2num(get(handles.edit_dissipation,'String'));  
DtCs                            = ET.Dispersion.Delta_tauCs;
k                               = ET.Dissipation.k;
% --------
cla(hAxes1)
axes(hAxes1)
set(hAxes1,'YColor','k','Box','on','ygrid','on','xgrid','on');
    set(hAxes1,'FontName','Helvetica','FontSize',fontSize1,'LineWidth',1)
hold on
pop_plot = get(handles.pop_plot,'Value');
switch pop_plot
    case 1
        tMax = 2*DtCs;
        if DtCs == 0
            DtCs = eps;
            tMax = 1e-3;
        end
        N = 2001;
        tSp  = linspace(-tMax, tMax, N); 
        dt   = tSp(2) - tSp(1);
        set(hAxes1,'xlim',[-tMax tMax],'xTick',-tMax/2:tMax/2:tMax/2,'xticklabel',{'','',''},...
        'YAxisLocation','left','Color','w');
        hXlabel = xlabel(hAxes1,'Time [s]','Color','k','Interpreter','LaTex','FontSize',fontSize1);
        hYlabel = ylabel(hAxes1,'$E_C^{outlet}$ ','Color','k','Interpreter','LaTex','FontSize',fontSize1);
        set(hXlabel,'Unit','points');
        posXlabel = get(hXlabel,'position');
        set(hXlabel,'position',[83 -25 0]);
        set(hYlabel,'Unit','points');
        posYlabel = get(hYlabel,'position');
        set(hYlabel,'position',[-30 70 0]);
        switch ET.pop_type_model
            case 1
                Eout = 0*tSp;
                Eout((N+1)./2) = 1;
                Tr = 1;
                text( 0,-0.15.*Tr,'$\tau_C^{s}$',...
                    'FontSize',fontSize2,'interpreter','latex',...
                    'parent',hAxes1,'horizontalAlignment','center',...
                    'rotation',0)
                text( -tMax*1.1,Tr,'$\delta $',...
                    'FontSize',fontSize2,'interpreter','latex',...
                    'parent',hAxes1,'horizontalAlignment','center',...
                    'rotation',0)
                
                set(hAxes1,'ylim',[0 1.2],'yTick',[0 1],'yticklabel',{'0',''});
            case 2
                Eout = 1/(pi^0.5*DtCs).*exp(-(tSp./DtCs).^2);
                Tr = 1/(pi^0.5*DtCs);
                text(-tMax./2,-0.15*Tr,'$\tau_C^{s} - \Delta \tau_C^{s}$',...
                    'FontSize',fontSize2,'interpreter','latex',...
                    'parent',hAxes1,'horizontalAlignment','center',...
                    'rotation',60)
                text( tMax./2,-0.15*Tr,'$\tau_C^{s} + \Delta \tau_C^{s}$',...
                    'FontSize',fontSize2,'interpreter','latex',...
                    'parent',hAxes1,'horizontalAlignment','center',...
                    'rotation',60)
                text( 0,-0.15*Tr,'$\tau_C^{s}$',...
                    'FontSize',fontSize2,'interpreter','latex',...
                    'parent',hAxes1,'horizontalAlignment','center',...
                    'rotation',0)
                text( -tMax*1.1,Tr,'$1\Big/(\sqrt{\pi} \Delta \tau_C^s) $',...
                    'FontSize',fontSize2,'interpreter','latex',...
                    'parent',hAxes1,'horizontalAlignment','center',...
                    'rotation',0)
                set(hAxes1,'ylim',[0 1.2*Tr],'yTick',[0 Tr],'yticklabel',{'0',''});
            case 3
                Eout = 0*tSp;
                Eout(ceil((tMax - DtCs)./dt):ceil((tMax + DtCs)./dt)) = 1./(2*DtCs);
                Tr = 1./(2*DtCs);
                text(-tMax./2,-0.15*Tr,'$\tau_C^{s} - \Delta \tau_C^{s}$',...
                    'FontSize',fontSize2,'interpreter','latex',...
                    'parent',hAxes1,'horizontalAlignment','center',...
                    'rotation',60)
                text( tMax./2,-0.15*Tr,'$\tau_C^{s} + \Delta \tau_C^{s}$',...
                    'FontSize',fontSize2,'interpreter','latex',...
                    'parent',hAxes1,'horizontalAlignment','center',...
                    'rotation',60)
                text( 0,-0.15*Tr,'$\tau_C^{s}$',...
                    'FontSize',fontSize2,'interpreter','latex',...
                    'parent',hAxes1,'horizontalAlignment','center',...
                    'rotation',0)
                text( -tMax*1.1,Tr,'$1\Big/{2 \Delta \tau_C^s} $',...
                    'FontSize',fontSize2,'interpreter','latex',...
                    'parent',hAxes1,'horizontalAlignment','center',...
                    'rotation',0)
                set(hAxes1,'ylim',[0 1.2*Tr],'yTick',[0 Tr],'yticklabel',{'0',''});
        end
        hLine = plot(hAxes1,tSp,k.*Eout,'-','color','b','linewidth',2);
    case 2
        fSp = -2000:1:2000;
        s = j*2*pi*fSp;
        set(hAxes1,'xlim',[fSp(1) fSp(end)]./1e3,'xTick',linspace(fSp(1),fSp(end),5)./1e3,...
            'xticklabel',num2str(linspace(fSp(1),fSp(end),5)'./1e3),...
        'YAxisLocation','left','Color','w');
        hXlabel = xlabel(hAxes1,'Frequency [kHz]','Color','k','Interpreter','LaTex','FontSize',fontSize1);
        hYlabel = ylabel(hAxes1,'$\widetilde{\mathcal{E}}$ ','Color','k','Interpreter','LaTex','FontSize',fontSize1);
        set(hXlabel,'Unit','points');
        posXlabel = get(hXlabel,'position');
        set(hXlabel,'position',[83 -25 0]);
        set(hYlabel,'Unit','points');
        posYlabel = get(hYlabel,'position');
        set(hYlabel,'position',[-30 70 0]);
        switch ET.pop_type_model
            case 1
                Etf = fSp.*0 + 1;
                Tr = 1;
            case 2
                Etf = abs(exp((DtCs.*s).^2./4));
                Tr = 1;
            case 3
                Etf = abs(sinc(2*pi*fSp.*DtCs./pi));
                Tr = 1;
        end
        hLine = plot(hAxes1,fSp./1e3,k.*Etf,'-','color','b','linewidth',2);
        set(hAxes1,'ylim',[0 1.2*Tr],'yTick',[0 Tr],'yticklabel',{'0','1'});
end
hold off
assignin('base','ET',ET);                   % save the current information to the workspace
set(handles.pb_SaveFig,                     'enable','on'); 
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_INI_BC_Entropy_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
try
varargout{1} = handles.output;
end


% --- Executes on button press in pb_OK.
function pb_OK_Callback(hObject, eventdata, handles)
% hObject    handle to pb_OK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global CI
global ET
Fcn_GUI_INI_BC_Entropy_update_plot(hObject)
CI.BC.ET                    = ET;
clear ET
CI.IsRun.GUI_INI_BC_Entropy = 1;
delete(handles.figure);

% --- Executes on button press in pb_Apply.
function pb_Apply_Callback(hObject, eventdata, handles)
% hObject    handle to pb_Apply (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Fcn_GUI_INI_BC_Entropy_update_plot(hObject)

% --- Executes on button press in pb_Cancel.
function pb_Cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pb_Cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(handles.figure);

% --- Executes on button press in pb_SaveFig.
function pb_SaveFig_Callback(hObject, eventdata, handles)
handles=guidata(hObject);
Fig = figure;
copyobj(handles.axes1, Fig);
set(Fig,        'units','points')
posFig = get(handles.figure,'position');
hAxes = get(Fig,'children');
set(hAxes(1),       'units','points',...
                    'position',[60 60 200 150],...
                    'ActivePositionProperty','position')
posAxesOuter = [0 0 300 250];
set(Fig,        'units','points',...
                'position', [posFig(1)+0.5*posFig(3)-0.5*posAxesOuter(3),...
                            posFig(2)+0.5*posFig(4)-0.5*posAxesOuter(4),...
                            posAxesOuter(3:4)]) 

% --- Executes on selection change in pop_Disp_model.
function pop_Disp_model_Callback(hObject, eventdata, handles)
% hObject    handle to pop_Disp_model (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pop_Disp_model contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_Disp_model
pop_type_model      = get(hObject,'Value');;
switch pop_type_model
    case 1
        set(handles.edit_Delta_tauCs,               'enable','off',...
                                                    'value',0);
    otherwise
        set(handles.edit_Delta_tauCs,               'enable','on');
end

% --- Executes during object creation, after setting all properties.
function pop_Disp_model_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_Disp_model (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Delta_tauCs_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Delta_tauCs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Delta_tauCs as text
%        str2double(get(hObject,'String')) returns contents of edit_Delta_tauCs as a double
datEdit = str2double(get(hObject, 'String'));
if isnan(datEdit)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% --- Executes during object creation, after setting all properties.
function edit_Delta_tauCs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Delta_tauCs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_dissipation_Callback(hObject, eventdata, handles)
% hObject    handle to edit_dissipation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_dissipation as text
%        str2double(get(hObject,'String')) returns contents of edit_dissipation as a double
datEdit = str2double(get(hObject, 'String'));
if isnan(datEdit)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% --- Executes during object creation, after setting all properties.
function edit_dissipation_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_dissipation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object deletion, before destroying properties.
function figure_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to figure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when user attempts to close figure.
function figure_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
delete(hObject);


% --- Executes on selection change in pop_plot.
function pop_plot_Callback(hObject, eventdata, handles)
% hObject    handle to pop_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pop_plot contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_plot


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


% --------------------------------------------------------------------
function uipushtool1_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);
Fig = figure;
copyobj(handles.axes1, Fig);
set(Fig,        'units','points')
posFig = get(handles.figure,'position');
hAxes = get(Fig,'children');
set(hAxes(1),       'units','points',...
                    'position',[60 60 200 150],...
                    'ActivePositionProperty','position')
posAxesOuter = [0 0 300 250];
set(Fig,        'units','points',...
                'position', [posFig(1)+0.5*posFig(3)-0.5*posAxesOuter(3),...
                            posFig(2)+0.5*posFig(4)-0.5*posAxesOuter(4),...
                            posAxesOuter(3:4)]) 
