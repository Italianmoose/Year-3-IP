function Fcn_GUI_INI_TP_Main_calculation(varargin)
% This function is used to calculate the mean properties in every section
% first created: 2014-04-01
% last modified: 2015-06-03
% author: Jingxuan LI (jingxuan.li@imperial.ac.uk)
hObject     = varargin{1};
handles     = guidata(hObject);
% ----------------------------
global CI
CI.TP.index_gamma = get(handles.pop_gamma,'Value');    % flag for gamma 
gamma_cst = 1.4;
Cp_cst = 1005;

CI.TP.numSection    = length(CI.CD.x_sample)-1;                 % number of sections
% 
% --------------------------inlet calculation -----------------------------
%
% The inlet gas is considered as air, since the gamma and W of the mixture
% is nearly the same as those of air, with a error of 5%.
% The burned gases are processed by the program to calculate the mass
% fraction of the mixture
CI.TP.p_mean(1:2,1)     = str2num(get(handles.edit_TP_p1,'string'));
CI.TP.T_mean(1:2,1)     = str2num(get(handles.edit_TP_T1,'string'));
% This information is given in the main program
% CI.Ru    = 8.3145;
% CI.W_air = 28.96512;
% CI.R_air = CI.Ru./CI.W_air.*1000;
% calculate  Cp and gamma of air in unburned gases
% use the program:
% [ci,co,DeltaHr,Cp_o]=Fcn_calculation_c_q_air(Ti,To)
switch CI.TP.index_gamma
    case 1
        CI.TP.c_mean(1:2,1)     = Fcn_c_gamma_cst(CI.TP.T_mean(1,1));
        CI.TP.Cp(1:2,1)         = Cp_cst;
        CI.TP.gamma(1:2,1)      = gamma_cst;
    case 2
        [CI.TP.c_mean(1:2,1),temp1,temp2,CI.TP.Cp(1:2,1)]...
        = Fcn_calculation_c_q_air(CI.TP.T_mean(1,1));
        CI.TP.gamma(1:2,1)  = CI.TP.Cp(1,1)./(CI.TP.Cp(1,1) - CI.R_air);            % gamma
end
% It is necessary to indicate that T_mean(a1:a2,b), herein, a1 denotes the 
% final state in the current interface, and a2 denotes the intermediate state
% b means the index number of the section
%
%
CI.TP.M1_u1_style   = get(handles.pop_TP_M1_u1,'Value');                    % Heat addition style
%
switch CI.TP.M1_u1_style
    case 1
        CI.TP.M_mean(1:2,1)     = str2num(get(handles.edit_TP_M1,'string')); % mean Mach number
        CI.TP.u_mean(1:2,1)     = CI.TP.M_mean(1,1).*CI.TP.c_mean(1,1);      % mean velocity
    case 2
        CI.TP.u_mean(1:2,1)     = str2num(get(handles.edit_TP_M1,'string'));
        CI.TP.M_mean(1:2,1)     = CI.TP.u_mean(1,1)./CI.TP.c_mean(1,1);
end
CI.TP.rho_mean(1:2,1)           = CI.TP.p_mean(1,1)./(CI.R_air.*CI.TP.T_mean(1,1)); % mean density
% 
% -----------calculate the properties in different sections----------------
%
for ss = 1:CI.TP.numSection-1 
    CI.TP.Theta(ss)  = (CI.CD.r_sample(ss+1)./CI.CD.r_sample(ss)).^2;       % Surface area ratio S2/S1
end
%
% --------------------------------
% Begin adding by Dong Yang
Liner_Flag = 1; % Initialize the flag of the Liner considered by the loop
HR_Flag = 1;    % Initialize the flag of the HR considered by the loop
% End adding by Dong Yang
% --------------------------------
indexHA_num = 0;                % set the initial value to zero and it will be increased by 1 after every HA interface
for ss = 1:CI.TP.numSection-1 
    % In every interface, the changes are splitted to two steps:
    % 1. cross sectional surface area change
    % 2. HA or .....
    % --------------step 1-------------------------
    %
    [   CI.TP.p_mean(1:2,ss+1),...
        CI.TP.rho_mean(1:2,ss+1),...
        CI.TP.u_mean(1:2,ss+1) ]... 
      = Fcn_calculation_TP_mean_WO_HeatAddition(CI.TP.p_mean(1,ss),...
                                                CI.TP.rho_mean(1,ss),...
                                                CI.TP.u_mean(1,ss),...
                                                CI.TP.Theta(ss),...
                                                CI.TP.gamma(1,ss));
    % ----------
    CI.TP.gamma(1:2,ss+1)   = CI.TP.gamma(1,ss);
    CI.TP.Cp(1:2,ss+1)      = CI.TP.Cp(1,ss);
    CI.TP.T_mean(1:2,ss+1)  = CI.TP.gamma(1,ss+1)/(CI.TP.gamma(1,ss+1)-1)...
                             *CI.TP.p_mean(1,ss+1)./(CI.TP.Cp(1,ss+1).*CI.TP.rho_mean(1,ss+1));
    CI.TP.c_mean(1:2,ss+1)  = ((CI.TP.gamma(1,ss+1) - 1).*CI.TP.Cp(1,ss+1).*CI.TP.T_mean(1,ss+1)).^0.5;
    CI.TP.M_mean(1:2,ss+1)  = CI.TP.u_mean(1,ss+1)./CI.TP.c_mean(1,ss+1);
    % ----------
    %
    % --------------step 2-------------------------
    %
    switch CI.CD.SectionIndex(ss+1)
        case 0
            % in case 0, no changes
        case {10,11}                                    % with HA
            indexHA_num = indexHA_num + 1;              % this number is increased by 1
            % ---------first calculate the temperature, Cp, Rg after the
            % heat addition interface -------------------------------------
            [   CI.TP.TRatio(indexHA_num),...
                CI.TP.c_mean(1,ss+1),...
                CI.TP.DeltaHr(indexHA_num),...          % heat release rate per mass flow rate
                CI.TP.gamma(1,ss+1),...
                CI.TP.Cp(1,ss+1)] =...
                Fcn_GUI_INI_TP_calculation_products_after_HA(   handles,...
                                                                indexHA_num,...
                                                                CI.TP.T_mean(2,ss+1),...
                                                                CI.TP.p_mean(2,ss+1));
            CI.TP.T_mean(1,ss+1)    = CI.TP.TRatio(indexHA_num).*CI.TP.T_mean(2,ss+1);
            switch CI.TP.index_gamma
                case 1
                    CI.TP.gamma(1,ss+1)     = gamma_cst;
                    CI.TP.Cp(1,ss+1)        = Cp_cst;
                    CI.TP.c_mean(1,ss+1)    = Fcn_c_gamma_cst(CI.TP.T_mean(1,ss+1));                
                case 2
                    % nothing happens!
            end
            Rg2 = CI.TP.Cp(1,ss+1)./(CI.TP.gamma(1,ss+1)./(CI.TP.gamma(1,ss+1)-1));
                    % ---------then, use the resolved temperature, Rg and the mean
                    % properties after the area changes to calculate the mean
                    % properties after HA ----------------------------------------
                    [   CI.TP.p_mean(1,ss+1),...               
                        CI.TP.rho_mean(1,ss+1),...
                        CI.TP.u_mean(1,ss+1)] = ...
                    Fcn_calculation_TP_mean_W_HeatAddition( CI.TP.p_mean(2,ss+1),...
                                                            CI.TP.rho_mean(2,ss+1),...
                                                            CI.TP.u_mean(2,ss+1),...
                                                            Rg2,...
                                                            CI.TP.T_mean(1,ss+1));
                    CI.TP.M_mean(1,ss+1)  = CI.TP.u_mean(1,ss+1)./CI.TP.c_mean(1,ss+1);
            switch CI.TP.index_gamma
                case 1  % constant Cp and gamma
                    CI.TP.DeltaHr(indexHA_num) = CI.TP.Cp(1,ss+1).*(CI.TP.T_mean(1,ss+1) - CI.TP.T_mean(2,ss+1))...
                        + 0.5*(CI.TP.u_mean(1,ss+1).^2 - CI.TP.u_mean(1,ss).^2);              
                case 2
                    % nothing happens!
            end
                    mass = CI.TP.rho_mean(2,ss+1).*CI.TP.u_mean(2,ss+1).*CI.CD.r_sample(ss+1).^2.*pi; % mass flow rate before HA
                    CI.TP.Q(indexHA_num)  = CI.TP.DeltaHr(indexHA_num).*mass;       % heat release rate
                    
        % Begin adding by Dong Yang
        case 2            
            % Typically, it can be assumed that mean flow will not be affect by the HR, anyway, 
            % the following procedure is an alternate way to calculate mean parameters across the crosssection where the resonator is installed.      
            
            section_Num=ss;
            [   CI.TP.T_mean(1,ss+1),...               
                CI.TP.rho_mean(1,ss+1),...
                CI.TP.u_mean(1,ss+1),...
                CI.TP.Cp(1,ss+1),...
                CI.TP.gamma(1,ss+1)] = ...
            Fcn_calculation_TP_mean_across_HR( HR_Flag,...
                                               section_Num,...
                                               CI.TP.T_mean(2,ss+1),...
                                               CI.TP.rho_mean(2,ss+1),...
                                               CI.TP.u_mean(2,ss+1));
                CI.TP.p_mean(1:2,ss+1)  = CI.TP.T_mean(1,ss+1).*CI.TP.rho_mean(1,ss+1)*(CI.TP.gamma(1,ss+1)-1)./CI.TP.gamma(1,ss+1)*CI.TP.Cp(1,ss+1);
                CI.TP.c_mean(1:2,ss+1)  = ((CI.TP.gamma(1,ss+1) - 1).*CI.TP.Cp(1,ss+1).*CI.TP.T_mean(1,ss+1)).^0.5;
                CI.TP.M_mean(1:2,ss+1)  = CI.TP.u_mean(1,ss+1)./CI.TP.c_mean(1,ss+1); 
            HR_Flag = HR_Flag+1;
        case 30
            % Nothing happens because mean flow properties are constant
            % across the inlet side cross-section of the lined duct
        case 31
            % It is assumed only mean flow velocity will increase a bit in
            % the lined duct region.
            [ CI.TP.u_mean(1,ss+1)] = Fcn_calculation_TP_mean_across_Liner(Liner_Flag,ss+1);
            CI.TP.M_mean(1,ss+1)=CI.TP.u_mean(1,ss+1)/CI.TP.c_mean(1,ss+1);
            Liner_Flag = Liner_Flag+1;
        % End adding by Dong Yang                
    end
end
%
% ---------after the calculation, we need to update the first HA-----------
%
guidata(hObject, handles);
Fcn_GUI_INI_TP_HA_Para2UI(hObject)
handles = guidata(hObject);
guidata(hObject, handles);
%
assignin('base','CI',CI)


function c = Fcn_c_gamma_cst(T)
gamma = 1.4;
Rs = 287;
c = (gamma*Rs*T).^0.5;
%
% -------------------------------end --------------------------------------