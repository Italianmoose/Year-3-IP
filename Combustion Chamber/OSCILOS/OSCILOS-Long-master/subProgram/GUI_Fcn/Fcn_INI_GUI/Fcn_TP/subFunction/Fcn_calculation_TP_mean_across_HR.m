function [T_mean2,rho_mean2,u_mean2,Cp_2,gamma_2]...
        = Fcn_calculation_TP_mean_across_HR(HR_Num, section_Num,T_mean1,rho_mean1,u_mean1)
% This function is developed to calculate mean flow parameters after the
% HR.
global CI

T_HR                      =CI.CD.HR_config(HR_Num,7);                       % Mean temperature in the cavity
[c_mean,temp1,temp2,Cp_HR]= Fcn_calculation_c_q_air(T_HR);                  % Calculate Cp in the cavity based on the given mean temperature || Is this right for burnt gas?
gamma_HR                  = Cp_HR/(Cp_HR - CI.R_air);                       % Using liner temperature to calculate Gamma outside the liner
c_HR                      =sqrt(gamma_HR*CI.R_air*T_HR);                    % Sound speed in the HR cavity
v_HR                      =CI.CD.HR_config(HR_Num,8)*c_HR;                  % Mean velocity in the neck
S_2                       =CI.CD.HR_config(HR_Num,5);                       % Neck area
P_HR =CI.TP.p_mean(1,(section_Num))/(1-1/2/CI.R_air/T_HR*(v_HR/0.64)^2);    % Calculate mean pressure in the cavity. The discharge coefficient is assumed to be 0.64 according to Falkovich (Fluid mechanics, 2011)
rho_v                     =P_HR/CI.R_air/T_HR;                              % Mean density in the cavity
S_1                       =pi*(CI.CD.r_sample(section_Num))^2;              % Cross sectional area of the duct

mass1_mean                =rho_mean1*u_mean1*S_1;                           % Mean mass flux in duct 1                 
mass_neck_mean            =rho_v*v_HR*S_2;                                  % Mean mass flux in HR neck
mass2_mean                =mass1_mean+mass_neck_mean;                       % Mean mass flux in duct 2

enthalpy_stag_mean1       =CI.TP.Cp(1,section_Num)*T_mean1+0.5*u_mean1^2;   % Mean stagnation enthalpy in duct 1
enthalpy_stag_mean_neck   =Cp_HR*T_HR+0.5*v_HR^2;                           % Mean stagnation enthalpy in HR neck
enthalpy_stag_mean2       =(enthalpy_stag_mean1*mass1_mean+enthalpy_stag_mean_neck*mass_neck_mean)/mass2_mean; % Mean stagnation enthalpy in duct 2
Momen_mean1               =CI.TP.p_mean(1,section_Num)+CI.TP.rho_mean(1,section_Num)*CI.TP.u_mean(1,section_Num)^2; % Momentum flux in duct 1
R_2                       =CI.TP.Cp(1,section_Num)*(CI.TP.gamma(1,section_Num)-1)/CI.TP.gamma(1,section_Num);

global X
X.mass2     =mass2_mean;
X.enthalpy2 =enthalpy_stag_mean2;
X.Momen2    =Momen_mean1;
X.S_1       =S_1;
X.Cp        =CI.TP.Cp(1,section_Num);
X.R_2       =R_2;
assignin('base','X',X); 

%% Calculation procedure
options = optimoptions('fsolve','Display','off');  % Turn off display
x0=[CI.TP.rho_mean(1,section_Num); CI.TP.u_mean(1,section_Num); CI.TP.T_mean(1,section_Num)];
[x,fval]=fsolve(@myfun,x0, options);

rho_mean2 =x(1);
u_mean2   =x(2);
T_mean2   =x(3);
Cp_2      =CI.TP.Cp(1,section_Num);
gamma_2   =CI.TP.gamma(1,section_Num);



function F=myfun(x)
 global X
 F = [x(1)*x(2)*X.S_1-X.mass2;
      X.Cp*x(3)+0.5*x(2)^2-X.enthalpy2;
      x(1)*X.R_2*x(3)+x(1)*x(2)^2-X.Momen2];

%--------------------------------end---------------------------------------%


      