function [Mtr_HR_33] = Fcn_calculation_Oscillations_across_HR(section_Num,HR_Num,i_omega,A_1_p,A_1_n)
% This function is used to calculate the frequency domain oscillations
% after the HR resonator cross-section.

% section number denotes the section before the HR
global CI

if CI.CD.HR_config(HR_Num,3)==0  % Linear HR
%% Begin the linear HR model
T_HR                      =CI.CD.HR_config(HR_Num,7);                       % Mean temperature in the cavity
[c_mean,temp1,temp2,Cp_HR]= Fcn_calculation_c_q_air(T_HR);                  % Calculate Cp in the cavity based on the given mean temperature || Is this right for burnt gas?
gamma                     = Cp_HR/(Cp_HR - CI.R_air);                       % Using liner temperature to calculate Gamma outside the liner
c_HR                      =sqrt(gamma*CI.R_air*T_HR);                       % Sound speed in the HR cavity
v_HR                      =CI.CD.HR_config(HR_Num,8)*c_HR;                  % Mean velocity in the neck
V_HR                      =CI.CD.HR_config(HR_Num,6);                       % Cavity Volume of the HR
P_HR =CI.TP.p_mean(1,section_Num)/(1-1/2/CI.R_air/T_HR*(v_HR/0.64)^2);  % Calculate mean pressure in the cavity. The discharge coefficient is assumed to be 0.64 according to Falkovich (Fluid mechanics, 2011)
rho_v=P_HR/CI.R_air/T_HR;                                                   % Mean density in the cavity

%% Calculate Rayleigh conductivity, then neck mouth impedance
S_2    =CI.CD.HR_config(HR_Num,5);
r_H    =sqrt(S_2/pi);                                                       % Neck radius
l_H    =CI.CD.HR_config(HR_Num,4);                                          % Neck length
St     =i_omega/1i*r_H/v_HR;                                                % Strouhal number (!!! The sign of this expression needs to be checked very carefully)
Rayl_con_value_zerolength =2*r_H*(1+(pi/2*besseli(1,St)*exp(-St)-1i*besselk(1,St)*sinh(St))/(St*(pi/2*besseli(1,St)*exp(-St)+1i*besselk(1,St)*cosh(St)))); % Rayleigh conductivity value
Rayl_con_value_zerolength =real(Rayl_con_value_zerolength)-1i*imag(Rayl_con_value_zerolength); 
Rayl_con_value_revised    =1/(1/(Rayl_con_value_zerolength)+l_H/(pi*r_H^2));% Revised Rayleigh conductivity
% Impedance ratio of the neck mouth
Z_HR   =i_omega*rho_v*S_2/Rayl_con_value_revised+c_HR^2*rho_v*S_2/(i_omega*V_HR);

%% Oscillation transfer matrix form left to right
c_1=CI.TP.c_mean(1,section_Num);                                        % Sound speed in the upstream duct
c_2=CI.TP.c_mean(1,section_Num+1);                                            % Sound speed in the downstream duct
rho_1=CI.TP.rho_mean(1,section_Num);                                    % Density in the upstream duct
rho_2=CI.TP.rho_mean(1,section_Num+1);                                        % Density in the downstream duct
u_1=CI.TP.u_mean(1,section_Num);                                        % Velocity in the upstream duct
u_2=CI.TP.u_mean(1,section_Num+1);                                            % Velocity in the downstream duct
gamma_1=CI.TP.gamma(1,section_Num);                                     % gamma in the upstream duct
gamma_2=CI.TP.gamma(1,section_Num+1);                                         % gamma in the downstream duct

p_1=CI.TP.p_mean(1,section_Num);                                        % p in the upstream duct
p_2=CI.TP.p_mean(1,section_Num+1);                                            % p in the downstream duct
S_1=pi*(CI.CD.r_sample(section_Num))^2;                                   % Cross sectional area of the duct

% A full stagnation enthalpy conservation equation is used and a transfer 
% matrix and an inverse transfer matrix are used to solve the conservation 
% equations

% Flux oscillations from the left duct
T_Wave_to_Flow_1=[1, 1, 0;
                  1/c_1^2, 1/c_1^2, -1/c_1^2;
                  1/rho_1/c_1, -1/rho_1/c_1, 0];
T_Flow_to_Flux_1=[0, u_1, rho_1;
                  1, u_1^2, 2*rho_1*u_1;
                  gamma_1/(gamma_1-1)*u_1, u_1^3/2, (gamma_1/(gamma_1-1)*p_1+1.5*rho_1*u_1^2)];
Osc_duct_factor_matrix1=S_1*T_Flow_to_Flux_1*T_Wave_to_Flow_1;
                            
% Flux oscillations from the neck              
Mass_neck_osc_factor        =-rho_v*S_2/Z_HR;
Mass_neck_osc_factor_m_mean =(rho_v*S_2*v_HR)*(c_HR^2*rho_v*S_2/Z_HR/(i_omega*V_HR)/rho_v); 
Mass_neck_osc_factor_m_osc  =-(rho_v*S_2/Z_HR)*(Cp_HR*T_HR+0.5*v_HR^2);                                
Enthalpy_neck_osc_factor    =Mass_neck_osc_factor_m_mean+Mass_neck_osc_factor_m_osc; % Enthalpy_neck_osc=Enthalpy_neck_osc_factor*(A_1_p+A_1_n)
                             
Osc_neck_factor_matrix=[Mass_neck_osc_factor, Mass_neck_osc_factor,0;
                        0, 0, 0;
                        Enthalpy_neck_osc_factor, Enthalpy_neck_osc_factor, 0];

% Total flux oscillations for the right duct                                        
Matrix_wave_1_to_Flux_2=Osc_duct_factor_matrix1+Osc_neck_factor_matrix; % Flux_2=Matrix_wave_1_to_Flux_2*[A_1_p; A_1_n; A_E_1];


T_Wave_to_Flow_2=[1, 1, 0;
                  1/c_2^2, 1/c_2^2, -1/c_2^2;
                  1/rho_2/c_2, -1/rho_2/c_2, 0];
T_Flow_to_Flux_2=[0, u_2, rho_2;
                  1, u_2^2, 2*rho_2*u_2;
                  gamma_2/(gamma_2-1)*u_2, u_2^3/2, (gamma_2/(gamma_2-1)*p_2+1.5*rho_2*u_2^2)];
Osc_duct_factor_matrix2=S_1*T_Flow_to_Flux_2*T_Wave_to_Flow_2;              

%% Final results
Mtr_HR_33=Osc_duct_factor_matrix2\Matrix_wave_1_to_Flux_2;

%% End the liner HR model


else  % Nonlinear HR
%% Begin nonlinea HR model   
   %% Pressure at the nonlinear HR neck mouth
   p_xh=A_1_p+A_1_n; 
   
   %% Calculate the nonlinear u_HR
   % Parameter settings
   S_1        =pi*(CI.CD.r_sample(section_Num))^2;                          % Cross sectional area of the duct
   S_2        =CI.CD.HR_config(HR_Num,5);                                   % Cross-sectional neck area
   V_HR       =CI.CD.HR_config(HR_Num,6);                                   % Cavity Volume of the HR
   alpha      =CI.CD.HR_config(HR_Num,9);                                   % Discharge coefficient of the this HR
   T_HR       =CI.CD.HR_config(HR_Num,7);                                   % Mean temperature in the cavity
   [c_mean,temp1,temp2,Cp_HR]= Fcn_calculation_c_q_air(T_HR);               % Calculate Cp in the cavity based on the given mean temperature || Is this right for burnt gas?
   gamma      = Cp_HR/(Cp_HR - CI.R_air);                                   % Using liner temperature to calculate Gamma outside the liner
   c_HR       =sqrt(gamma*CI.R_air*T_HR);                                   % Sound speed in the HR cavity
   P_HR       =CI.TP.p_mean(1,(section_Num));                               % Mean pressure in the cavity is the same with that in the duct
   rho_v      =P_HR/CI.R_air/T_HR;                                          % Air density in the cavity
   r_H=sqrt(S_2/pi);                                                        % Neck radius
   l_H=CI.CD.HR_config(HR_Num,4);                                           % Neck length
   l_eff      =l_H+2*0.68*r_H; %16*r_H/3/pi;                                % Effective neck length. Mass end correction is assumed to be 8/3/pi*D.
   %    Use fsolve to calculate \widetilde{u}_{HR} 
    u_HR_0=10;  % The initial guess
    options = optimset('Display','off');        % the calculation results by fsolve will not be shown in the workspace
    [x,fval,exitflag]=fsolve(@(x) p_xh-x*(alpha*rho_v*abs(x)+i_omega*rho_v*l_eff+c_HR^2*rho_v*S_2/(i_omega*V_HR)),u_HR_0, options);
    u_HR=x;
   % Impedance ratio of the neck mouth
    Z_HR=alpha*rho_v*abs(u_HR)+i_omega*rho_v*l_eff+c_HR^2*rho_v*S_2/(i_omega*V_HR);
    
   %% Oscillation transfer matrix from left to right
   % Mean flow parameters up and down streams
   c_1    =CI.TP.c_mean(1,section_Num);                                     % Sound speed in the upstream duct
   c_2    =CI.TP.c_mean(1,section_Num+1);                                   % Sound speed in the downstream duct
   rho_1  =CI.TP.rho_mean(1,section_Num);                                   % Density in the upstream duct
   rho_2  =CI.TP.rho_mean(1,section_Num+1);                                 % Density in the downstream duct
   u_1    =CI.TP.u_mean(1,section_Num);                                     % Velocity in the upstream duct
   u_2    =CI.TP.u_mean(1,section_Num+1);                                   % Velocity in the downstream duct
   gamma_1=CI.TP.gamma(1,section_Num);                                      % gamma in the upstream duct
   gamma_2=CI.TP.gamma(1,section_Num+1);                                    % gamma in the downstream duct
   p_1    =CI.TP.p_mean(1,section_Num);                                     % p in the upstream duct
   p_2    =CI.TP.p_mean(1,section_Num+1);                                   % p in the downstream duct
  %% Full stagnation enthalpy conservation equations
  % Flux oscillations from the left duct
  T_Wave_to_Flow_1=[1, 1, 0;
                    1/c_1^2, 1/c_1^2, -1/c_1^2;
                    1/rho_1/c_1, -1/rho_1/c_1, 0];
  T_Flow_to_Flux_1=[0, u_1, rho_1;
                    1, u_1^2, 2*rho_1*u_1;
                    gamma_1/(gamma_1-1)*u_1, u_1^3/2, (gamma_1/(gamma_1-1)*p_1+1.5*rho_1*u_1^2)];
  Osc_duct_factor_matrix1=S_1*T_Flow_to_Flux_1*T_Wave_to_Flow_1;
                            
  % Flux oscillations from the neck              
  Mass_neck_osc_factor        =-rho_v*S_2/Z_HR;                          
  Enthalpy_neck_osc_factor    =-(rho_v*S_2/Z_HR)*(Cp_HR*T_HR); % Enthalpy_neck_osc=Enthalpy_neck_osc_factor*(A_1_p+A_1_n)
                             
  Osc_neck_factor_matrix      =[Mass_neck_osc_factor, Mass_neck_osc_factor,0;
                                0, 0, 0;
                                Enthalpy_neck_osc_factor, Enthalpy_neck_osc_factor, 0];

  % Total flux oscillations for the right duct                                        
  Matrix_wave_1_to_Flux_2     =Osc_duct_factor_matrix1+Osc_neck_factor_matrix; % Flux_2=Matrix_wave_1_to_Flux_2*[A_1_p; A_1_n; A_E_1];

  T_Wave_to_Flow_2            =[1, 1, 0;
                                1/c_2^2, 1/c_2^2, -1/c_2^2;
                                1/rho_2/c_2, -1/rho_2/c_2, 0];
  T_Flow_to_Flux_2            =[0, u_2, rho_2;
                                1, u_2^2, 2*rho_2*u_2;
                                gamma_2/(gamma_2-1)*u_2, u_2^3/2, (gamma_2/(gamma_2-1)*p_2+1.5*rho_2*u_2^2)];
  Osc_duct_factor_matrix2     =S_1*T_Flow_to_Flux_2*T_Wave_to_Flow_2;              

  %% Final results
  Mtr_HR_33=Osc_duct_factor_matrix2\Matrix_wave_1_to_Flux_2;
  
  %% End Nonlinear HR model           
end

%----------------------------------end-------------------------------------
