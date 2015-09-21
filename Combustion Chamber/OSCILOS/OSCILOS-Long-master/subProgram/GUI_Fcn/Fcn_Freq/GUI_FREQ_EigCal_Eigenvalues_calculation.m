function GUI_FREQ_EigCal_Eigenvalues_calculation(varargin)
% This function is used to calculate the eigenvalues
hObject = varargin{1};
handles = guidata(hObject);
global CI
% ------------------------split the scan domain and contour domain --------
%
GUI_FREQ_EigCal_Scan_domain_Contour_domain_splitting
%
% ------------------------Initialize some nonlinear parameters --------
%
Fcn_Para_initialization
%
% -------------------------% calculation depend on nonlinearities----------
% 
switch CI.EIG.APP_style
    case {11,12}  
        Fcn_calculation_APP_1(hObject)
    case {21,22}                             % nonlinear flame model with enforcing before the flame
        Fcn_calculation_APP_2(hObject)
    case {3}
        Fcn_calculation_APP_3(hObject)       % nonlinear system with enforcing at the inlet
end
handles = guidata(hObject);
guidata(hObject, handles);
set(handles.pop_numMode,        'enable','on');  
set(handles.pop_PlotType,       'enable','on');
set(handles.pb_AOC_Plot,        'enable','on');
set(handles.pb_AOC_SaveFig,     'enable','on');
set(handles.pb_AOC_OK,          'enable','on');
% --
eigenvalue = CI.EIG.Scan.EigValCol{1};
for k = 1:length(eigenvalue)
   StringMode{k} = ['Mode number: ' num2str(k)]; 
end
set(handles.pop_numMode,    'string',StringMode);
%
guidata(hObject, handles);
%
% -------------------------------------------------------------------------
%
function Fcn_calculation_APP_1(varargin)
% CI.EIG.APP_style = {11, 12} 
% linear
hObject = varargin{1};
handles = guidata(hObject);
global CI
hWaitBar = waitbar(0,'The calculation may take several minutes, please wait...');
% for ss =1:CI.EIG.FDF.uRatioNum 
ss = 1;
    CI.EIG.Scan.EigValCol{ss}   = Fcn_calculation_eigenvalues(1);
    CI.EIG.Cont.ValCol{ss}      = Fcn_calculation_contour(1);
waitbar(ss/CI.EIG.FDF.uRatioNum);
drawnow
% end
close(hWaitBar)
% ---------------update the table---------------
data_num(:,1)   = abs(imag(CI.EIG.Scan.EigValCol{1})./2./pi);
data_num(:,2)   = real(CI.EIG.Scan.EigValCol{1});
data_cell       = num2cell(data_num);
set(handles.uitable,'data',data_cell);         % Update the table
% ---------------update the table---------------
% ---------------update popup plot -------------
 set(handles.pop_PlotType,...
                'string',{'Map of eigenvalues';'Modeshape'},...
                'enable','on'); 
% ---------------update popup plot -------------
assignin('base','CI',CI);
guidata(hObject, handles);
% -------------------------------------------------------------------------
%
function Fcn_calculation_APP_2(varargin)
% CI.EIG.APP_style = {21, 22} 
% nonlinear
hObject = varargin{1};
handles = guidata(hObject);
global CI
global FDF
global HP
% -----------------------get flame describing function --------------------
%
if  CI.EIG.APP_style == 21     % From flame model, CI.FM.indexFM == 2
    % use nonlinear flame model
    uRatioMin           = str2double(get(handles.edit_uRatio_min,'string'));
    uRatioMax           = str2double(get(handles.edit_uRatio_max,'string'));
    uRatioNum           = str2double(get(handles.edit_uRatio_SampNum,'string'));
    CI.EIG.FDF.uRatioSp = linspace(uRatioMin,uRatioMax,uRatioNum);
    HP = CI.FM.HP{CI.FM.indexMainHPinHp};
    % nonlinear model style: HP.NL.style
    switch HP.NL.style
    case 1
        for ss = 1:uRatioNum
            CI.EIG.FDF.num{ss}   = HP.FTF.num;
            CI.EIG.FDF.den{ss}   = HP.FTF.den;
            CI.EIG.FDF.tauf(ss)  = HP.FTF.tauf;     
        end  
    case 2
        for ss = 1:uRatioNum
            CI.EIG.FDF.num{ss}   = HP.FTF.num;
            CI.EIG.FDF.den{ss}   = HP.FTF.den;
            CI.EIG.FDF.tauf(ss)  = HP.FTF.tauf;     
            % something must be done to include the nonlinear in the Fcn_DetEqn.m
        end  
    case 3      
        CI.EIG.FDF.Lf       = interp1(  HP.NL.Model3.uRatio,...
                                        HP.NL.Model3.Lf,...
                                        CI.EIG.FDF.uRatioSp,'linear','extrap');         % Nonlinear ratio
        CI.EIG.FDF.taufNSp  = HP.NL.Model3.taufN.*(1-CI.EIG.FDF.Lf);                % Nonlinear time delay
        for ss = 1:uRatioNum
            CI.EIG.FDF.num{ss}   = CI.EIG.FDF.Lf(ss).*HP.FTF.num;
            CI.EIG.FDF.den{ss}   = HP.FTF.den;
            CI.EIG.FDF.tauf(ss)  = HP.FTF.tauf + CI.EIG.FDF.taufNSp(ss);     % Time delay of the FDF model
        end    
    end
elseif CI.EIG.APP_style == 22  % From experimental FDF
    
    HP = CI.FM.HP{CI.FM.indexMainHPinHp};
    CI.EIG.FDF.uRatioSp             = HP.FMEXP.uRatio;
    uRatioNum            = length(CI.EIG.FDF.uRatioSp);
    for ss = 1:uRatioNum
        FTF = HP.FMEXP.FTF{ss};
        CI.EIG.FDF.num{ss}  = FTF.num;
        CI.EIG.FDF.den{ss}  = FTF.den;
        CI.EIG.FDF.tauf(ss) = FTF.tau_correction;  % Be careful: we always use a = a0 exp(-a*tau)
    end
end
CI.EIG.FDF.uRatioNum = uRatioNum;
assignin('base','CI',CI);
assignin('base','HP',HP);
%
% ---------------calculate the eigenvalues---------------
%
hWaitBar = waitbar(0,'The calculation may take several minutes, please wait...');
for ss =1:CI.EIG.FDF.uRatioNum
    FDF.uRatio  = CI.EIG.FDF.uRatioSp(ss);
    FDF.num     = CI.EIG.FDF.num{ss};
    FDF.den     = CI.EIG.FDF.den{ss};
    FDF.tauf    = CI.EIG.FDF.tauf(ss);
    assignin('base','FDF',FDF);
    CI.EIG.Scan.EigValCol{ss}   = Fcn_calculation_eigenvalues(2);
    CI.EIG.Cont.ValCol{ss}      = Fcn_calculation_contour(2);
    waitbar(ss/CI.EIG.FDF.uRatioNum);
    drawnow
end
close(hWaitBar)
assignin('base','CI',CI);
% ---------------update the table---------------
data_num(:,1)   = abs(imag(CI.EIG.Scan.EigValCol{1})./2./pi);
data_num(:,2)   = real(CI.EIG.Scan.EigValCol{1});
data_cell       = num2cell(data_num);
set(handles.uitable,'data',data_cell);         % Update the table
% ---------------update the table---------------
%
% ---------------update the slider---------------
if CI.EIG.FDF.uRatioNum == 1
    set(handles.slider_uRatio,'visible','off')
else
set(handles.slider_uRatio,...
                        'Enable','on',...
                        'min',1,...
                        'max',CI.EIG.FDF.uRatioNum,...
                        'value',1,...
                        'SliderStep',[1/(CI.EIG.FDF.uRatioNum-1), 1/(CI.EIG.FDF.uRatioNum-1)]);
end
set(handles.edit_uRatio,'string',num2str(CI.EIG.FDF.uRatioSp(1)));
% ---------------update the slider---------------
% ---------------update popup plot -------------
set(handles.pop_PlotType,...
                'string',{'Map of eigenvalues';'Modeshape'; 'Evolution of eigenvalue with velocity ratio'},...
                'enable','on'); 
% ---------------update popup plot -------------
guidata(hObject, handles);
%
%--------------------------------------------------------------------------

function Fcn_calculation_APP_3(varargin)
% CI.EIG.APP_style = {3} 
% nonlinear
hObject = varargin{1};
handles = guidata(hObject);
global CI
global FDF
% -----------------------get flame describing function --------------------
    % use nonlinear flame model
    pRatioMin           = str2double(get(handles.edit_pRatio_min,'string'));
    pRatioMax           = str2double(get(handles.edit_pRatio_max,'string'));
    pRatioNum           = str2double(get(handles.edit_pRatio_SampNum,'string'));
    CI.EIG.FDF.pRatioSp = linspace(pRatioMin,pRatioMax,pRatioNum);
    CI.EIG.FDF.pRatioNum = pRatioNum;
    CI.EIG.PR.R1abs     = abs(polyval(CI.BC.num1,0)./polyval(CI.BC.den1,0));
    CI.EIG.PR.A1minusSp = CI.EIG.FDF.pRatioSp.*10^(-5).*CI.TP.p_mean(1,1)./(1+CI.EIG.PR.R1abs);
assignin('base','CI',CI);
%
% ---------------calculate the eigenvalues---------------
%
hWaitBar = waitbar(0,'The calculation may take several minutes, please wait...');
for ss =1:CI.EIG.FDF.pRatioNum
    FDF.pRatio  = CI.EIG.FDF.pRatioSp(ss);
%     FDF.num     = CI.EIG.FDF.num{ss};
%     FDF.den     = CI.EIG.FDF.den{ss};
%     FDF.tauf    = CI.EIG.FDF.tauf(ss);
    CI.EIG.PR.A1minus_this_step=CI.EIG.PR.A1minusSp(ss);
    assignin('base','CI',CI);
    assignin('base','FDF',FDF);
    CI.EIG.Scan.EigValCol{ss}   = Fcn_calculation_eigenvalues(3);
    CI.EIG.Cont.ValCol{ss}      = Fcn_calculation_contour(3);
    waitbar(ss/CI.EIG.FDF.pRatioNum);
    drawnow
end
close(hWaitBar)
assignin('base','CI',CI);
% ---------------update the table---------------
data_num(:,1)   = abs(imag(CI.EIG.Scan.EigValCol{1})./2./pi);
data_num(:,2)   = real(CI.EIG.Scan.EigValCol{1});
data_cell       = num2cell(data_num);
set(handles.uitable,'data',data_cell);         % Update the table
% ---------------update the table---------------
%
% ---------------update the slider---------------
if CI.EIG.FDF.pRatioNum == 1
    set(handles.slider_pRatio,'visible','off')
else
set(handles.slider_pRatio,...
                        'Enable','on',...
                        'min',1,...
                        'max',CI.EIG.FDF.pRatioNum,...
                        'value',1,...
                        'SliderStep',[1/(CI.EIG.FDF.pRatioNum-1), 1/(CI.EIG.FDF.pRatioNum-1)]);
end
set(handles.edit_pRatio,'string',num2str(CI.EIG.FDF.pRatioSp(1)));
% ---------------update the slider---------------
% ---------------update popup plot -------------
set(handles.pop_PlotType,...
                'string',{'Map of eigenvalues';'Modeshape'; 'Evolution of eigenvalue with pressure ratio'},...
                'enable','on'); 
% ---------------update popup plot -------------
guidata(hObject, handles);
%
% -------------------------------------------------------------------------
%
function Fcn_Para_initialization
global CI
CI.EIG.FDF.uRatioNum    = 1;
CI.EIG.FDF.uRatioSp(1)  = 0;
CI.EIG.FDF.pRatioNum    = 1;
CI.EIG.FDF.pRatioSp(1)  = 0;
CI.EIG.FDF.num{1}       = [];
CI.EIG.FDF.den{1}       = [];
CI.EIG.FDF.tauf(1)      = 0;
assignin('base','CI',CI);
%
% -------------------------------------------------------------------------
%
