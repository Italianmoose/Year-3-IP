function [Osc_trans_Mtr_across_Liner] = Fcn_calculation_Osc_trans_matrixes_across_Liner(section_Num,Liner_Num,i_omega)
% This function is used to calculate the oscillating transfer matrixes
% across the liner region. These matrixes are used for plotting waves in
% the lined region.

global CI

% section_Num is the number of the liner section

omega=i_omega/1i; % Angular velocity
NUM_x_d_sections=50;


% Duct length and radius
Rp          =CI.CD.r_sample(1,section_Num);    % Radius of duct
C1          =2*pi*Rp;                          % Circumference of the duct
Sp          =pi*Rp^2;                          % Cross-sectional area of the duct
c_duct      =CI.TP.c_mean(1,section_Num);      % Mean sound speed in the duct
M_duct_l    =CI.TP.M_mean(1,section_Num);      % Mean Mach number at the left side of the duct
M_duct_r    =CI.TP.M_mean(1,section_Num+1);    % Mean Mach number at the right side of the duct
rho_duct    =CI.TP.rho_mean(1,section_Num);    % Mean density in the duct
Liner_length=CI.CD.Liner_config(Liner_Num,3);  % Length of the perforated liners
% Discetilize the Liner
dx          =Liner_length/NUM_x_d_sections;
x_step_all  =round(Liner_length/dx+1);

% Mean flow velocity in the duct
s_steps=1:1:x_step_all;
x_locations=dx.*(s_steps-1);


%% Consider the four types of Liners separately
if CI.CD.Liner_config(Liner_Num,4)==0  % Single liner
%% Begin single Liner model
     if CI.CD.Liner_config(Liner_Num,5)==0  % With large cavity outside boundary
         %% Begin Type 1
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         % Parameters of the liner 1
         a1=CI.CD.Liner_config(Liner_Num,8);      % Radius of the perforated holes in liner 1
         d1=CI.CD.Liner_config(Liner_Num,9);      % Hole distance in liner 1
         t1=CI.CD.Liner_config(Liner_Num,11);     % Thickness of the liner 1
         sigma1=pi*a1^2/d1^2;                     % Open area ratio on liner 1
         
         T_Liner                      =CI.CD.Liner_config(Liner_Num,7);                       % Mean temperature in the liners
         [c_mean,temp1,temp2,Cp_Liner]   = Fcn_calculation_c_q_air(T_Liner);                  % Calculate Cp in the Liner based on the given mean temperature || Is this right for burnt gas?
         gamma_Liner                     = Cp_Liner/(Cp_Liner - CI.R_air);                    % Using liner temperature to calculate Gamma outside the liner
         c_Liner                      =sqrt(gamma_Liner*CI.R_air*T_Liner);                    % Sound speed in the HR cavity
         v_Liner1                      =CI.CD.Liner_config(Liner_Num,10)*c_Liner;             % Mean velocity in the Liner hole
         P_Liner =CI.TP.p_mean(1,section_Num)/(1-1/2/CI.R_air/T_Liner*(v_Liner1/0.64)^2);     % Calculate mean pressure in the cavity. The discharge coefficient is assumed to be 0.64 according to Falkovich (Fluid mechanics, 2011)
         rho_Liner=P_Liner/CI.R_air/T_Liner;                                                  % Mean density in the cavity
         u_x(1:x_step_all)=CI.TP.u_mean(1,section_Num)+C1*x_locations./Sp.*(v_Liner1*sigma1);    % Mean velocity distribution in the lined duct
         
        
         % Calculate compliance of Liner 1 
         St1=omega*a1/v_Liner1;
         Ray_gamma_num_Liner1=besseli(1,St1)^2*(1+1/St1)+4/pi^2*exp(2*St1)*cosh(St1)*besselk(1,St1)^2*(cosh(St1)-sinh(St1)/St1);
         Ray_delta_num_Liner1=2/pi/St1*besseli(1,St1)*besselk(1,St1)*exp(2*St1);
         Ray_Den_Liner1      =besseli(1,St1)^2+4/pi^2*exp(2*St1)*cosh(St1)^2*besselk(1,St1)^2;
         Rayl_con_value_zerolength_Liner1=2*a1*(Ray_gamma_num_Liner1/Ray_Den_Liner1+1i*Ray_delta_num_Liner1/Ray_Den_Liner1);
         Rayl_con_value_revised_Liner1=1/(1/(Rayl_con_value_zerolength_Liner1)+t1/(pi*a1^2));    % Revised Rayleigh conductivity for Liner 1
         eta1=Rayl_con_value_revised_Liner1/d1^2;  % Compliance
         %% Transfer matrixes from the begin to any point across the Liner  
         M_x_0_L(1).mat=1;
         for x_step=1:1:(x_step_all-1);
         
         u_x_forthis_M = u_x(x_step);
         M_T1_x_step_11=-(1i*omega/(c_duct+u_x_forthis_M)+(C1*c_duct*eta1*rho_Liner)/(2*1i*omega*Sp*rho_duct));
         M_T1_x_step_12=-(C1*c_duct*eta1*rho_Liner)/(2*1i*omega*Sp*rho_duct);
         M_T1_x_step_21=(C1*c_duct*eta1*rho_Liner)/(2*1i*omega*Sp*rho_duct);
         M_T1_x_step_22=1i*omega/(c_duct-u_x_forthis_M)+(C1*c_duct*eta1*rho_Liner)/(2*1i*omega*Sp*rho_duct);
         
         M_T1_x_step=[M_T1_x_step_11, M_T1_x_step_12; M_T1_x_step_21, M_T1_x_step_22];
         I_T1=[1,0;
               0,1];
         
         %%%% Forth order Runge-kutta stepping method
         % First step
         M1_x_step= M_T1_x_step;
         % Second step
         u_x_forthis_M = (u_x(x_step)+u_x(x_step+1))/2;
         M_T1_x_stepandhalf_11=-(1i*omega/(c_duct+u_x_forthis_M)+(C1*c_duct*eta1*rho_Liner)/(2*1i*omega*Sp*rho_duct));
         M_T1_x_stepandhalf_12=-(C1*c_duct*eta1*rho_Liner)/(2*1i*omega*Sp*rho_duct);
         M_T1_x_stepandhalf_21=(C1*c_duct*eta1*rho_Liner)/(2*1i*omega*Sp*rho_duct);
         M_T1_x_stepandhalf_22=1i*omega/(c_duct-u_x_forthis_M)+(C1*c_duct*eta1*rho_Liner)/(2*1i*omega*Sp*rho_duct);
         M_x_step_and_half=[M_T1_x_stepandhalf_11, M_T1_x_stepandhalf_12; M_T1_x_stepandhalf_21, M_T1_x_stepandhalf_22];
    
         M2_x_step=M_x_step_and_half*(I_T1+dx/2*M1_x_step);
         % Third step
         M3_x_step=M_x_step_and_half*(I_T1+dx/2*M2_x_step);
         % Forth step
         u_x_forthis_M = u_x(x_step+1);
         M_T1_x_nextstep_11=-(1i*omega/(c_duct+u_x_forthis_M)+(C1*c_duct*eta1*rho_Liner)/(2*1i*omega*Sp*rho_duct));
         M_T1_x_nextstep_12=-(C1*c_duct*eta1*rho_Liner)/(2*1i*omega*Sp*rho_duct);
         M_T1_x_nextstep_21=(C1*c_duct*eta1*rho_Liner)/(2*1i*omega*Sp*rho_duct);
         M_T1_x_nextstep_22=1i*omega/(c_duct-u_x_forthis_M)+(C1*c_duct*eta1*rho_Liner)/(2*1i*omega*Sp*rho_duct);
         M_x_step_next_step=[M_T1_x_nextstep_11, M_T1_x_nextstep_12; M_T1_x_nextstep_21, M_T1_x_nextstep_22];
         
         M4_x_step=M_x_step_next_step*(I_T1+dx*M3_x_step);
    
         % Transfer matrix between parameters at x_step and (x_step+1)
         M_x_step_transf=I_T1+dx/6*(M1_x_step+2*M2_x_step+2*M3_x_step+M4_x_step);
         if x_step==1
             M_x_0_L(x_step).mat=M_x_0_L(x_step).mat*M_x_step_transf;
         else
             M_x_0_L(x_step).mat=M_x_0_L(x_step-1).mat*M_x_step_transf;
         end
         M_x_0_L_type1(x_step).mat=M_x_0_L(x_step).mat;
         
         end
         
         %% Link oscillations at the begin and end of the Liner
         for x_step=1:1:(x_step_all-1);
         A_to_phi_l(x_step).mat=[(1+M_duct_l)/c_duct^2,0;
                                 0, (1-M_duct_l)/c_duct^2];
         A_to_phi_r(x_step).mat=[(1+u_x(x_step+1)/c_duct)/c_duct^2,0;
                                 0, (1-u_x(x_step+1)/c_duct)/c_duct^2];
         
         Mtr_Liner_22  =A_to_phi_r(x_step).mat\M_x_0_L_type1(x_step).mat*A_to_phi_l(x_step).mat;
         %tau_c(x_step) =x_locations(x_step)/CI.TP.u_mean(1,section_Num);    % Entropy wave convecting time
         Osc_trans_Mtr_across_Liner(x_step).mat=[Mtr_Liner_22(1,1), Mtr_Liner_22(1,2), 0;
                                                 Mtr_Liner_22(2,1), Mtr_Liner_22(2,2), 0;
                                                 0                , 0,                 1];
         end

         
         %% End Type 1
         %%%%%%%%%%%%%%%%%%%%%%
     else if CI.CD.Liner_config(Liner_Num,5)==1 % With rigid wall boundary
             %% Begin Type 2
         a1=CI.CD.Liner_config(Liner_Num,8);      % Radius of the perforated holes in liner 1
         d1=CI.CD.Liner_config(Liner_Num,9);      % Hole distance in liner 1
         t1=CI.CD.Liner_config(Liner_Num,11);     % Thickness of the liner 1
         Rw=CI.CD.Liner_config(Liner_Num,6);      % Radius of the outside annular wall
         Sp1=pi*(Rw^2-(Rp+t1)^2);                 % Cross-sectional area of the annular gap between liner1 and the rigid wall
         
         sigma1=pi*a1^2/d1^2;                     % Open area ratio on liner 1
         
         T_Liner                         =CI.CD.Liner_config(Liner_Num,7);                    % Mean temperature in the liners
         [c_mean,temp1,temp2,Cp_Liner]   = Fcn_calculation_c_q_air(T_Liner);                  % Calculate Cp in the Liner based on the given mean temperature || Is this right for burnt gas?
         gamma_Liner                     = Cp_Liner/(Cp_Liner - CI.R_air);                    % Using liner temperature to calculate Gamma outside the liner
         c_Liner                         =sqrt(gamma_Liner*CI.R_air*T_Liner);                 % Sound speed in the HR cavity
         v_Liner1                        =CI.CD.Liner_config(Liner_Num,10)*c_Liner;           % Mean velocity in the Liner hole
         P_Liner =CI.TP.p_mean(1,section_Num)/(1-1/2/CI.R_air/T_Liner*(v_Liner1/0.64)^2);     % Calculate mean pressure in the cavity. The discharge coefficient is assumed to be 0.64 according to Falkovich (Fluid mechanics, 2011)
         rho_Liner=P_Liner/CI.R_air/T_Liner;                                                  % Mean density in the cavity
         u_x(1:x_step_all)=CI.TP.u_mean(1,section_Num)+C1*x_locations./Sp.*(v_Liner1*sigma1);    % Mean velocity distribution in the lined duct
         

         % Calculate compliance of Liner 1 
         St1=omega*a1/v_Liner1;
         Ray_gamma_num_Liner1=besseli(1,St1)^2*(1+1/St1)+4/pi^2*exp(2*St1)*cosh(St1)*besselk(1,St1)^2*(cosh(St1)-sinh(St1)/St1);
         Ray_delta_num_Liner1=2/pi/St1*besseli(1,St1)*besselk(1,St1)*exp(2*St1);
         Ray_Den_Liner1      =besseli(1,St1)^2+4/pi^2*exp(2*St1)*cosh(St1)^2*besselk(1,St1)^2;
         Rayl_con_value_zerolength_Liner1=2*a1*(Ray_gamma_num_Liner1/Ray_Den_Liner1+1i*Ray_delta_num_Liner1/Ray_Den_Liner1);
         Rayl_con_value_revised_Liner1=1/(1/(Rayl_con_value_zerolength_Liner1)+t1/(pi*a1^2));    % Revised Rayleigh conductivity for Liner 1
         eta1=Rayl_con_value_revised_Liner1/d1^2;  % Compliance


         %% Transfer matrix from the begin to the end of the Liner  
         M_x_0_L(1).mat=1;
         for x_step=1:1:(x_step_all-1);
         
         u_x_forthis_M = u_x(x_step);
         M_T2_x_step_11=-(1i*omega/(c_duct+u_x_forthis_M)+(C1*c_duct*eta1*rho_Liner)/(2*1i*omega*Sp*rho_duct));
         M_T2_x_step_12=-(C1*c_duct*eta1*rho_Liner)/(2*1i*omega*Sp*rho_duct);
         M_T2_x_step_13=(C1*eta1*rho_Liner)/(2*1i*omega*c_duct*Sp*rho_duct);
         M_T2_x_step_14=0;
         M_T2_x_step_21=(C1*c_duct*eta1*rho_Liner)/(2*1i*omega*Sp*rho_duct);
         M_T2_x_step_22=1i*omega/(c_duct-u_x_forthis_M)+(C1*c_duct*eta1*rho_Liner)/(2*1i*omega*Sp*rho_duct);
         M_T2_x_step_23=-(C1*eta1*rho_Liner)/(2*1i*omega*c_duct*Sp*rho_duct);
         M_T2_x_step_24=0;
         M_T2_x_step_31=0;
         M_T2_x_step_32=0;
         M_T2_x_step_33=0;
         M_T2_x_step_34=-1i*omega;
         M_T2_x_step_41=(C1*eta1*c_duct^2)/(1i*omega*Sp1);
         M_T2_x_step_42=(C1*eta1*c_duct^2)/(1i*omega*Sp1);
         M_T2_x_step_43=-((C1*eta1)/(1i*omega*Sp1)+(1i*omega)/(c_Liner^2));
         M_T2_x_step_44=0;         
         
         M_T2_x_step=[M_T2_x_step_11, M_T2_x_step_12, M_T2_x_step_13, M_T2_x_step_14;
                      M_T2_x_step_21, M_T2_x_step_22, M_T2_x_step_23, M_T2_x_step_24;
                      M_T2_x_step_31, M_T2_x_step_32, M_T2_x_step_33, M_T2_x_step_34;
                      M_T2_x_step_41, M_T2_x_step_42, M_T2_x_step_43, M_T2_x_step_44;];
         I_T2=[1,0,0,0;
               0,1,0,0;
               0,0,1,0;
               0,0,0,1];
           
         %%%% Forth order Runge-kutta stepping method
         % First step
         M1_x_step= M_T2_x_step;
         % Second step
         u_x_forthis_M = (u_x(x_step)+u_x(x_step+1))/2;
         M_T2_x_step_11=-(1i*omega/(c_duct+u_x_forthis_M)+(C1*c_duct*eta1*rho_Liner)/(2*1i*omega*Sp*rho_duct));
         M_T2_x_step_12=-(C1*c_duct*eta1*rho_Liner)/(2*1i*omega*Sp*rho_duct);
         M_T2_x_step_13=(C1*eta1*rho_Liner)/(2*1i*omega*c_duct*Sp*rho_duct);
         M_T2_x_step_14=0;
         M_T2_x_step_21=(C1*c_duct*eta1*rho_Liner)/(2*1i*omega*Sp*rho_duct);
         M_T2_x_step_22=1i*omega/(c_duct-u_x_forthis_M)+(C1*c_duct*eta1*rho_Liner)/(2*1i*omega*Sp*rho_duct);
         M_T2_x_step_23=-(C1*eta1*rho_Liner)/(2*1i*omega*c_duct*Sp*rho_duct);
         M_T2_x_step_24=0;
         M_T2_x_step_31=0;
         M_T2_x_step_32=0;
         M_T2_x_step_33=0;
         M_T2_x_step_34=-1i*omega;
         M_T2_x_step_41=(C1*eta1*c_duct^2)/(1i*omega*Sp1);
         M_T2_x_step_42=(C1*eta1*c_duct^2)/(1i*omega*Sp1);
         M_T2_x_step_43=-((C1*eta1)/(1i*omega*Sp1)+(1i*omega)/(c_Liner^2));
         M_T2_x_step_44=0;         
         
         M_x_step_and_half=[M_T2_x_step_11, M_T2_x_step_12, M_T2_x_step_13, M_T2_x_step_14;
                            M_T2_x_step_21, M_T2_x_step_22, M_T2_x_step_23, M_T2_x_step_24;
                            M_T2_x_step_31, M_T2_x_step_32, M_T2_x_step_33, M_T2_x_step_34;
                            M_T2_x_step_41, M_T2_x_step_42, M_T2_x_step_43, M_T2_x_step_44;];
    
         M2_x_step=M_x_step_and_half*(I_T2+dx/2*M1_x_step);
         % Third step
         M3_x_step=M_x_step_and_half*(I_T2+dx/2*M2_x_step);
         % Forth step
         u_x_forthis_M = u_x(x_step+1);
         M_T2_x_step_11=-(1i*omega/(c_duct+u_x_forthis_M)+(C1*c_duct*eta1*rho_Liner)/(2*1i*omega*Sp*rho_duct));
         M_T2_x_step_12=-(C1*c_duct*eta1*rho_Liner)/(2*1i*omega*Sp*rho_duct);
         M_T2_x_step_13=(C1*eta1*rho_Liner)/(2*1i*omega*c_duct*Sp*rho_duct);
         M_T2_x_step_14=0;
         M_T2_x_step_21=(C1*c_duct*eta1*rho_Liner)/(2*1i*omega*Sp*rho_duct);
         M_T2_x_step_22=1i*omega/(c_duct-u_x_forthis_M)+(C1*c_duct*eta1*rho_Liner)/(2*1i*omega*Sp*rho_duct);
         M_T2_x_step_23=-(C1*eta1*rho_Liner)/(2*1i*omega*c_duct*Sp*rho_duct);
         M_T2_x_step_24=0;
         M_T2_x_step_31=0;
         M_T2_x_step_32=0;
         M_T2_x_step_33=0;
         M_T2_x_step_34=-1i*omega;
         M_T2_x_step_41=(C1*eta1*c_duct^2)/(1i*omega*Sp1);
         M_T2_x_step_42=(C1*eta1*c_duct^2)/(1i*omega*Sp1);
         M_T2_x_step_43=-((C1*eta1)/(1i*omega*Sp1)+(1i*omega)/(c_Liner^2));
         M_T2_x_step_44=0;         
         
         M_x_step_next_step=[M_T2_x_step_11, M_T2_x_step_12, M_T2_x_step_13, M_T2_x_step_14;
                             M_T2_x_step_21, M_T2_x_step_22, M_T2_x_step_23, M_T2_x_step_24;
                             M_T2_x_step_31, M_T2_x_step_32, M_T2_x_step_33, M_T2_x_step_34;
                             M_T2_x_step_41, M_T2_x_step_42, M_T2_x_step_43, M_T2_x_step_44;];
         M4_x_step=M_x_step_next_step*(I_T2+dx*M3_x_step);
    
         % Transfer matrix between parameters at x_step and (x_step+1)
         M_x_step_transf=I_T2+dx/6*(M1_x_step+2*M2_x_step+2*M3_x_step+M4_x_step);
         if x_step==1
             M_x_0_L(x_step).mat=M_x_0_L(x_step).mat*M_x_step_transf;
         else
             M_x_0_L(x_step).mat=M_x_0_L(x_step-1).mat*M_x_step_transf;
         end
         M_x_0_L_type2(x_step).mat=M_x_0_L(x_step).mat;
         
         end
         %% Link oscillations at the begin and end of the Liner
         for x_step=1:1:(x_step_all-1);
         M_x_0_L_type2_22(x_step).mat=[M_x_0_L(x_step).mat(1,1)-M_x_0_L(x_step).mat(1,3)*M_x_0_L(x_step_all-1).mat(4,1)/M_x_0_L(x_step_all-1).mat(4,3),   M_x_0_L(x_step).mat(1,2)-M_x_0_L(x_step).mat(1,3)*M_x_0_L(x_step_all-1).mat(4,2)/M_x_0_L(x_step_all-1).mat(4,3);
                                       M_x_0_L(x_step).mat(2,1)-M_x_0_L(x_step).mat(2,3)*M_x_0_L(x_step_all-1).mat(4,1)/M_x_0_L(x_step_all-1).mat(4,3),   M_x_0_L(x_step).mat(2,2)-M_x_0_L(x_step).mat(2,3)*M_x_0_L(x_step_all-1).mat(4,2)/M_x_0_L(x_step_all-1).mat(4,3)];
         A_to_phi_l(x_step).mat=[(1+M_duct_l)/c_duct^2,0;
                                 0, (1-M_duct_l)/c_duct^2];
         A_to_phi_r(x_step).mat=[(1+u_x(x_step+1)/c_duct)/c_duct^2,0;
                                 0, (1-u_x(x_step+1)/c_duct)/c_duct^2];
         
         Mtr_Liner_22=A_to_phi_r(x_step).mat\M_x_0_L_type2_22(x_step).mat*A_to_phi_l(x_step).mat;
         
         tau_c(x_step)=x_locations(x_step)/CI.TP.u_mean(1,section_Num);    % Entropy wave convecting time
         Osc_trans_Mtr_across_Liner(x_step).mat=[Mtr_Liner_22(1,1), Mtr_Liner_22(1,2), 0;
                                                 Mtr_Liner_22(2,1), Mtr_Liner_22(2,2), 0;
                                                 0                , 0,                 exp(-tau_c(x_step).*i_omega)];
         end   

             %% End Type 2
             %%%%%%%%%%%%%%%%
         end
         % Not probable to happen
     end
%% End single Liner model
else if CI.CD.Liner_config(Liner_Num,4)==1  % Double liners
     %% Begin double liners model 
     if CI.CD.Liner_config(Liner_Num,5)==0  % With large cavity outside boundary
         %% Begin Type 3
         
         a1=CI.CD.Liner_config(Liner_Num,12);      % Radius of the perforated holes in liner 1
         d1=CI.CD.Liner_config(Liner_Num,13);      % Hole distance in liner 1
         t1=CI.CD.Liner_config(Liner_Num,15);      % Thickness of the liner 1
         a2=CI.CD.Liner_config(Liner_Num,16);      % Radius of the perforated holes in liner 2
         d2=CI.CD.Liner_config(Liner_Num,17);      % Hole distance in liner 2
         R2=CI.CD.Liner_config(Liner_Num,18);      % Radius liner 2
         t2=CI.CD.Liner_config(Liner_Num,19);      % Thickness of the liner 2
         C2=2*pi*R2;                               % Circumference at the liner 2 inner face

         Sp1=pi*(R2^2-(Rp)^2);                  % Cross-sectional area of the annular gap between liner1 and the rigid wall
         
         sigma1=pi*a1^2/d1^2;                      % Open area ratio on liner 1
         sigma2=pi*a2^2/d2^2;                      % Open area ratio on liner 1
         
         T_Liner                         =CI.CD.Liner_config(Liner_Num,7);                    % Mean temperature in the liners
         [c_mean,temp1,temp2,Cp_Liner]   = Fcn_calculation_c_q_air(T_Liner);                  % Calculate Cp in the Liner based on the given mean temperature || Is this right for burnt gas?
         gamma_Liner                     =Cp_Liner/(Cp_Liner - CI.R_air);                     % Using liner temperature to calculate Gamma outside the liner
         c_Liner                         =sqrt(gamma_Liner*CI.R_air*T_Liner);                 % Sound speed in the HR cavity
         v_Liner1                        =CI.CD.Liner_config(Liner_Num,14)*c_Liner;           % Mean velocity in the Liner 1 hole
         v_Liner2                        =sigma1*v_Liner1/sigma2*C1/C2;                       % Mean velocity in the Liner 1 hole
         P_Liner =CI.TP.p_mean(1,section_Num)/(1-1/2/CI.R_air/T_Liner*(v_Liner1/0.64)^2);     % Calculate mean pressure in the cavity. The discharge coefficient is assumed to be 0.64 according to Falkovich (Fluid mechanics, 2011)
         rho_Liner=P_Liner/CI.R_air/T_Liner;                                                  % Mean density in the cavity
         u_x(1:x_step_all)=CI.TP.u_mean(1,section_Num)+C1*x_locations./Sp.*(v_Liner1*sigma1);   % Mean velocity distribution in the lined duct
         
         % Calculate compliance of Liner 1 and 2
         St1=omega*a1/v_Liner1;
         St2=omega*a2/v_Liner2;
         Ray_gamma_num_Liner1=besseli(1,St1)^2*(1+1/St1)+4/pi^2*exp(2*St1)*cosh(St1)*besselk(1,St1)^2*(cosh(St1)-sinh(St1)/St1);
         Ray_delta_num_Liner1=2/pi/St1*besseli(1,St1)*besselk(1,St1)*exp(2*St1);
         Ray_Den_Liner1      =besseli(1,St1)^2+4/pi^2*exp(2*St1)*cosh(St1)^2*besselk(1,St1)^2;
         Rayl_con_value_zerolength_Liner1=2*a1*(Ray_gamma_num_Liner1/Ray_Den_Liner1+1i*Ray_delta_num_Liner1/Ray_Den_Liner1);
         Rayl_con_value_revised_Liner1=1/(1/(Rayl_con_value_zerolength_Liner1)+t1/(pi*a1^2));    % Revised Rayleigh conductivity for Liner 1
         eta1=Rayl_con_value_revised_Liner1/d1^2;  % Compliance

         Ray_gamma_num_Liner2=besseli(1,St2)^2*(1+1/St2)+4/pi^2*exp(2*St2)*cosh(St2)*besselk(1,St2)^2*(cosh(St2)-sinh(St2)/St2);
         Ray_delta_num_Liner2=2/pi/St2*besseli(1,St2)*besselk(1,St2)*exp(2*St2);
         Ray_Den_Liner2      =besseli(1,St2)^2+4/pi^2*exp(2*St2)*cosh(St2)^2*besselk(1,St2)^2;
         Rayl_con_value_zerolength_Liner2=2*a2*(Ray_gamma_num_Liner2/Ray_Den_Liner2+1i*Ray_delta_num_Liner2/Ray_Den_Liner2);
         Rayl_con_value_revised_Liner2=1/(1/(Rayl_con_value_zerolength_Liner2)+t2/(pi*a2^2));    % Revised Rayleigh conductivity for Liner 1
         eta2=Rayl_con_value_revised_Liner2/d2^2;  % Compliance
         
         %% Transfer matrix from the begin to the end of the Liner  
         M_x_0_L(1).mat=1;
         for x_step=1:1:(x_step_all-1);
         
         u_x_forthis_M = u_x(x_step);
         M_T3_x_step_11=-(1i*omega/(c_duct+u_x_forthis_M)+(C1*c_duct*eta1*rho_Liner)/(2*1i*omega*Sp*rho_duct));
         M_T3_x_step_12=-(C1*c_duct*eta1*rho_Liner)/(2*1i*omega*Sp*rho_duct);
         M_T3_x_step_13=(C1*eta1*rho_Liner)/(2*1i*omega*c_duct*Sp*rho_duct);
         M_T3_x_step_14=0;
         M_T3_x_step_21=(C1*c_duct*eta1*rho_Liner)/(2*1i*omega*Sp*rho_duct);
         M_T3_x_step_22=1i*omega/(c_duct-u_x_forthis_M)+(C1*c_duct*eta1*rho_Liner)/(2*1i*omega*Sp*rho_duct);
         M_T3_x_step_23=-(C1*eta1*rho_Liner)/(2*1i*omega*c_duct*Sp*rho_duct);
         M_T3_x_step_24=0;
         M_T3_x_step_31=0;
         M_T3_x_step_32=0;
         M_T3_x_step_33=0;
         M_T3_x_step_34=-1i*omega;
         M_T3_x_step_41=(C1*eta1*c_duct^2)/(1i*omega*Sp1);
         M_T3_x_step_42=(C1*eta1*c_duct^2)/(1i*omega*Sp1);
         M_T3_x_step_43=-((C1*eta1)/(1i*omega*Sp1)+(C2*eta2)/(1i*omega*Sp1)+(1i*omega)/(c_Liner^2));
         M_T3_x_step_44=0;         
         
         M_T3_x_step=[M_T3_x_step_11, M_T3_x_step_12, M_T3_x_step_13, M_T3_x_step_14;
                      M_T3_x_step_21, M_T3_x_step_22, M_T3_x_step_23, M_T3_x_step_24;
                      M_T3_x_step_31, M_T3_x_step_32, M_T3_x_step_33, M_T3_x_step_34;
                      M_T3_x_step_41, M_T3_x_step_42, M_T3_x_step_43, M_T3_x_step_44;];
         I_T3=[1,0,0,0;
               0,1,0,0;
               0,0,1,0;
               0,0,0,1];
           
         %%%% Forth order Runge-kutta stepping method
         % First step
         M1_x_step= M_T3_x_step;
         % Second step
         u_x_forthis_M = (u_x(x_step)+u_x(x_step+1))/2;
         M_T3_x_step_11=-(1i*omega/(c_duct+u_x_forthis_M)+(C1*c_duct*eta1*rho_Liner)/(2*1i*omega*Sp*rho_duct));
         M_T3_x_step_12=-(C1*c_duct*eta1*rho_Liner)/(2*1i*omega*Sp*rho_duct);
         M_T3_x_step_13=(C1*eta1*rho_Liner)/(2*1i*omega*c_duct*Sp*rho_duct);
         M_T3_x_step_14=0;
         M_T3_x_step_21=(C1*c_duct*eta1*rho_Liner)/(2*1i*omega*Sp*rho_duct);
         M_T3_x_step_22=1i*omega/(c_duct-u_x_forthis_M)+(C1*c_duct*eta1*rho_Liner)/(2*1i*omega*Sp*rho_duct);
         M_T3_x_step_23=-(C1*eta1*rho_Liner)/(2*1i*omega*c_duct*Sp*rho_duct);
         M_T3_x_step_24=0;
         M_T3_x_step_31=0;
         M_T3_x_step_32=0;
         M_T3_x_step_33=0;
         M_T3_x_step_34=-1i*omega;
         M_T3_x_step_41=(C1*eta1*c_duct^2)/(1i*omega*Sp1);
         M_T3_x_step_42=(C1*eta1*c_duct^2)/(1i*omega*Sp1);
         M_T3_x_step_43=-((C1*eta1)/(1i*omega*Sp1)+(C2*eta2)/(1i*omega*Sp1)+(1i*omega)/(c_Liner^2));
         M_T3_x_step_44=0;         
         
         M_x_step_and_half=[M_T3_x_step_11, M_T3_x_step_12, M_T3_x_step_13, M_T3_x_step_14;
                            M_T3_x_step_21, M_T3_x_step_22, M_T3_x_step_23, M_T3_x_step_24;
                            M_T3_x_step_31, M_T3_x_step_32, M_T3_x_step_33, M_T3_x_step_34;
                            M_T3_x_step_41, M_T3_x_step_42, M_T3_x_step_43, M_T3_x_step_44;];
    
         M2_x_step=M_x_step_and_half*(I_T3+dx/2*M1_x_step);
         % Third step
         M3_x_step=M_x_step_and_half*(I_T3+dx/2*M2_x_step);
         % Forth step
         u_x_forthis_M = u_x(x_step+1);
         M_T3_x_step_11=-(1i*omega/(c_duct+u_x_forthis_M)+(C1*c_duct*eta1*rho_Liner)/(2*1i*omega*Sp*rho_duct));
         M_T3_x_step_12=-(C1*c_duct*eta1*rho_Liner)/(2*1i*omega*Sp*rho_duct);
         M_T3_x_step_13=(C1*eta1*rho_Liner)/(2*1i*omega*c_duct*Sp*rho_duct);
         M_T3_x_step_14=0;
         M_T3_x_step_21=(C1*c_duct*eta1*rho_Liner)/(2*1i*omega*Sp*rho_duct);
         M_T3_x_step_22=1i*omega/(c_duct-u_x_forthis_M)+(C1*c_duct*eta1*rho_Liner)/(2*1i*omega*Sp*rho_duct);
         M_T3_x_step_23=-(C1*eta1*rho_Liner)/(2*1i*omega*c_duct*Sp*rho_duct);
         M_T3_x_step_24=0;
         M_T3_x_step_31=0;
         M_T3_x_step_32=0;
         M_T3_x_step_33=0;
         M_T3_x_step_34=-1i*omega;
         M_T3_x_step_41=(C1*eta1*c_duct^2)/(1i*omega*Sp1);
         M_T3_x_step_42=(C1*eta1*c_duct^2)/(1i*omega*Sp1);
         M_T3_x_step_43=-((C1*eta1)/(1i*omega*Sp1)+(C2*eta2)/(1i*omega*Sp1)+(1i*omega)/(c_Liner^2));
         M_T3_x_step_44=0;         
         
         M_x_step_next_step=[M_T3_x_step_11, M_T3_x_step_12, M_T3_x_step_13, M_T3_x_step_14;
                             M_T3_x_step_21, M_T3_x_step_22, M_T3_x_step_23, M_T3_x_step_24;
                             M_T3_x_step_31, M_T3_x_step_32, M_T3_x_step_33, M_T3_x_step_34;
                             M_T3_x_step_41, M_T3_x_step_42, M_T3_x_step_43, M_T3_x_step_44;];
         M4_x_step=M_x_step_next_step*(I_T3+dx*M3_x_step);
    
         % Transfer matrix between parameters at x_step and (x_step+1)
         M_x_step_transf=I_T3+dx/6*(M1_x_step+2*M2_x_step+2*M3_x_step+M4_x_step);
         if x_step==1
             M_x_0_L(x_step).mat=M_x_0_L(x_step).mat*M_x_step_transf;
         else
             M_x_0_L(x_step).mat=M_x_0_L(x_step-1).mat*M_x_step_transf;
         end
         M_x_0_L_type3(x_step).mat=M_x_0_L(x_step).mat;
         
         end
         
         %% Link oscillations at the begin and end of the Liner
         for x_step=1:1:(x_step_all-1);
         M_x_0_L_type3_22(x_step).mat=[M_x_0_L(x_step).mat(1,1)-M_x_0_L(x_step).mat(1,3)*M_x_0_L(x_step_all-1).mat(4,1)/M_x_0_L(x_step_all-1).mat(4,3),   M_x_0_L(x_step).mat(1,2)-M_x_0_L(x_step).mat(1,3)*M_x_0_L(x_step_all-1).mat(4,2)/M_x_0_L(x_step_all-1).mat(4,3);
                                       M_x_0_L(x_step).mat(2,1)-M_x_0_L(x_step).mat(2,3)*M_x_0_L(x_step_all-1).mat(4,1)/M_x_0_L(x_step_all-1).mat(4,3),   M_x_0_L(x_step).mat(2,2)-M_x_0_L(x_step).mat(2,3)*M_x_0_L(x_step_all-1).mat(4,2)/M_x_0_L(x_step_all-1).mat(4,3)];
         A_to_phi_l(x_step).mat=[(1+M_duct_l)/c_duct^2,0;
                                 0, (1-M_duct_l)/c_duct^2];
         A_to_phi_r(x_step).mat=[(1+u_x(x_step+1)/c_duct)/c_duct^2,0;
                                 0, (1-u_x(x_step+1)/c_duct)/c_duct^2];
                             
         Mtr_Liner_22=A_to_phi_r(x_step).mat\M_x_0_L_type3_22(x_step).mat*A_to_phi_l(x_step).mat;
         
         
         tau_c(x_step)=x_locations(x_step)/CI.TP.u_mean(1,section_Num);    % Entropy wave convecting time
         Osc_trans_Mtr_across_Liner(x_step).mat=[Mtr_Liner_22(1,1), Mtr_Liner_22(1,2), 0;
                                                 Mtr_Liner_22(2,1), Mtr_Liner_22(2,2), 0;
                                                 0                , 0,                 exp(-tau_c(x_step).*i_omega)];
         end 

         %% End Type 3
         %%%%%%%%%%%%%%%%%%
         else if CI.CD.Liner_config(Liner_Num,5)==1 % With rigid wall boundary
                 %% Begin Type 4
         a1=CI.CD.Liner_config(Liner_Num,12);      % Radius of the perforated holes in liner 1
         d1=CI.CD.Liner_config(Liner_Num,13);      % Hole distance in liner 1
         t1=CI.CD.Liner_config(Liner_Num,15);      % Thickness of the liner 1
         a2=CI.CD.Liner_config(Liner_Num,16);      % Radius of the perforated holes in liner 2
         d2=CI.CD.Liner_config(Liner_Num,17);      % Hole distance in liner 2
         R2=CI.CD.Liner_config(Liner_Num,18);      % Radius liner 2
         t2=CI.CD.Liner_config(Liner_Num,19);      % Thickness of the liner 2
         C2=2*pi*R2;                               % Circumference at the liner 2 inner face
         Rw=CI.CD.Liner_config(Liner_Num,6);       % Radius of the outside annular wall
         Sp2=pi*(Rw^2-(R2+t2)^2);                  % Cross-sectional area of the annular gap between liner2 and the rigid wall
         Sp1=pi*(R2^2-(Rp+t1)^2);                  % Cross-sectional area of the annular gap between liner1 and the rigid wall
         
         sigma1=pi*a1^2/d1^2;                      % Open area ratio on liner 1
         sigma2=pi*a2^2/d2^2;                      % Open area ratio on liner 1
         
         T_Liner                         =CI.CD.Liner_config(Liner_Num,7);                    % Mean temperature in the liners
         [c_mean,temp1,temp2,Cp_Liner]   = Fcn_calculation_c_q_air(T_Liner);                  % Calculate Cp in the Liner based on the given mean temperature || Is this right for burnt gas?
         gamma_Liner                     = Cp_Liner/(Cp_Liner - CI.R_air);                    % Using liner temperature to calculate Gamma outside the liner
         c_Liner                         =sqrt(gamma_Liner*CI.R_air*T_Liner);                 % Sound speed in the HR cavity
         v_Liner1                        =CI.CD.Liner_config(Liner_Num,14)*c_Liner;           % Mean velocity in the Liner 1 hole
         v_Liner2                        =sigma1*v_Liner1/sigma2*C1/C2;                       % Mean velocity in the Liner 2 hole
         P_Liner =CI.TP.p_mean(1,section_Num)/(1-1/2/CI.R_air/T_Liner*(v_Liner1/0.64)^2);     % Calculate mean pressure in the cavity. The discharge coefficient is assumed to be 0.64 according to Falkovich (Fluid mechanics, 2011)
         rho_Liner=P_Liner/CI.R_air/T_Liner;                                                  % Mean density in the cavity
         u_x(1:x_step_all)=CI.TP.u_mean(1,section_Num)+C1*x_locations./Sp.*(v_Liner1*sigma1);   % Mean velocity distribution in the lined duct
         
         % Calculate compliance of Liner 1 and 2
         St1=omega*a1/v_Liner1;
         St2=omega*a2/v_Liner2;
         Ray_gamma_num_Liner1=besseli(1,St1)^2*(1+1/St1)+4/pi^2*exp(2*St1)*cosh(St1)*besselk(1,St1)^2*(cosh(St1)-sinh(St1)/St1);
         Ray_delta_num_Liner1=2/pi/St1*besseli(1,St1)*besselk(1,St1)*exp(2*St1);
         Ray_Den_Liner1      =besseli(1,St1)^2+4/pi^2*exp(2*St1)*cosh(St1)^2*besselk(1,St1)^2;
         Rayl_con_value_zerolength_Liner1=2*a1*(Ray_gamma_num_Liner1/Ray_Den_Liner1+1i*Ray_delta_num_Liner1/Ray_Den_Liner1);
         Rayl_con_value_revised_Liner1=1/(1/(Rayl_con_value_zerolength_Liner1)+t1/(pi*a1^2));    % Revised Rayleigh conductivity for Liner 1
         eta1=Rayl_con_value_revised_Liner1/d1^2;  % Compliance

         Ray_gamma_num_Liner2=besseli(1,St2)^2*(1+1/St2)+4/pi^2*exp(2*St2)*cosh(St2)*besselk(1,St2)^2*(cosh(St2)-sinh(St2)/St2);
         Ray_delta_num_Liner2=2/pi/St2*besseli(1,St2)*besselk(1,St2)*exp(2*St2);
         Ray_Den_Liner2      =besseli(1,St2)^2+4/pi^2*exp(2*St2)*cosh(St2)^2*besselk(1,St2)^2;
         Rayl_con_value_zerolength_Liner2=2*a2*(Ray_gamma_num_Liner2/Ray_Den_Liner2+1i*Ray_delta_num_Liner2/Ray_Den_Liner2);
         Rayl_con_value_revised_Liner2=1/(1/(Rayl_con_value_zerolength_Liner2)+t2/(pi*a2^2));    % Revised Rayleigh conductivity for Liner 1
         eta2=Rayl_con_value_revised_Liner2/d2^2;  % Compliance
         
         %% Transfer matrix from the begin to the end of the Liner  
         M_x_0_L(1).mat=1;
         for x_step=1:1:(x_step_all-1);
         
         u_x_forthis_M = u_x(x_step);
         M_T4_x_step_11=-(1i*omega/(c_duct+u_x_forthis_M)+(C1*c_duct*eta1*rho_Liner)/(2*1i*omega*Sp*rho_duct));
         M_T4_x_step_12=-(C1*c_duct*eta1*rho_Liner)/(2*1i*omega*Sp*rho_duct);
         M_T4_x_step_13=(C1*eta1*rho_Liner)/(2*1i*omega*c_duct*Sp*rho_duct);
         M_T4_x_step_14=0;
         M_T4_x_step_15=0;
         M_T4_x_step_16=0;
         
         M_T4_x_step_21=(C1*c_duct*eta1*rho_Liner)/(2*1i*omega*Sp*rho_duct);
         M_T4_x_step_22=1i*omega/(c_duct-u_x_forthis_M)+(C1*c_duct*eta1*rho_Liner)/(2*1i*omega*Sp*rho_duct);
         M_T4_x_step_23=-(C1*eta1*rho_Liner)/(2*1i*omega*c_duct*Sp*rho_duct);
         M_T4_x_step_24=0;
         M_T4_x_step_25=0;
         M_T4_x_step_26=0;
         
         M_T4_x_step_31=0;
         M_T4_x_step_32=0;
         M_T4_x_step_33=0;
         M_T4_x_step_34=-1i*omega;
         M_T4_x_step_35=0;
         M_T4_x_step_36=0;
         
         M_T4_x_step_41=(C1*eta1*c_duct^2)/(1i*omega*Sp1);
         M_T4_x_step_42=(C1*eta1*c_duct^2)/(1i*omega*Sp1);
         M_T4_x_step_43=-((C1*eta1)/(1i*omega*Sp1)+(C2*eta2)/(1i*omega*Sp1)+(1i*omega)/(c_Liner^2));
         M_T4_x_step_44=0;  
         M_T4_x_step_45=C2*eta2/(1i*omega*Sp1);
         M_T4_x_step_46=0;
         
         M_T4_x_step_51=0;
         M_T4_x_step_52=0;
         M_T4_x_step_53=0;
         M_T4_x_step_54=0;
         M_T4_x_step_55=0;
         M_T4_x_step_56=-1i*omega;
         
         M_T4_x_step_61=0;
         M_T4_x_step_62=0;
         M_T4_x_step_63=C2*eta2/(1i*omega*Sp2);
         M_T4_x_step_64=0;  
         M_T4_x_step_65=-(C2*eta2/(1i*omega*Sp2)+(1i*omega)/(c_Liner^2));
         M_T4_x_step_66=0;
         
         M_T4_x_step=[M_T4_x_step_11, M_T4_x_step_12, M_T4_x_step_13, M_T4_x_step_14, M_T4_x_step_15, M_T4_x_step_16;
                      M_T4_x_step_21, M_T4_x_step_22, M_T4_x_step_23, M_T4_x_step_24, M_T4_x_step_25, M_T4_x_step_26;
                      M_T4_x_step_31, M_T4_x_step_32, M_T4_x_step_33, M_T4_x_step_34, M_T4_x_step_35, M_T4_x_step_36;
                      M_T4_x_step_41, M_T4_x_step_42, M_T4_x_step_43, M_T4_x_step_44, M_T4_x_step_45, M_T4_x_step_46;
                      M_T4_x_step_51, M_T4_x_step_52, M_T4_x_step_53, M_T4_x_step_54, M_T4_x_step_55, M_T4_x_step_56;
                      M_T4_x_step_61, M_T4_x_step_62, M_T4_x_step_63, M_T4_x_step_64, M_T4_x_step_65, M_T4_x_step_66];
         I_T4=[1,0,0,0,0,0;
               0,1,0,0,0,0;
               0,0,1,0,0,0;
               0,0,0,1,0,0;
               0,0,0,0,1,0;
               0,0,0,0,0,1];
           
         %%%% Forth order Runge-kutta stepping method
         % First step
         M1_x_step= M_T4_x_step;
         % Second step
         u_x_forthis_M = (u_x(x_step)+u_x(x_step+1))/2;
         M_T4_x_step_11=-(1i*omega/(c_duct+u_x_forthis_M)+(C1*c_duct*eta1*rho_Liner)/(2*1i*omega*Sp*rho_duct));
         M_T4_x_step_12=-(C1*c_duct*eta1*rho_Liner)/(2*1i*omega*Sp*rho_duct);
         M_T4_x_step_13=(C1*eta1*rho_Liner)/(2*1i*omega*c_duct*Sp*rho_duct);
         M_T4_x_step_14=0;
         M_T4_x_step_15=0;
         M_T4_x_step_16=0;
         
         M_T4_x_step_21=(C1*c_duct*eta1*rho_Liner)/(2*1i*omega*Sp*rho_duct);
         M_T4_x_step_22=1i*omega/(c_duct-u_x_forthis_M)+(C1*c_duct*eta1*rho_Liner)/(2*1i*omega*Sp*rho_duct);
         M_T4_x_step_23=-(C1*eta1*rho_Liner)/(2*1i*omega*c_duct*Sp*rho_duct);
         M_T4_x_step_24=0;
         M_T4_x_step_25=0;
         M_T4_x_step_26=0;
         
         M_T4_x_step_31=0;
         M_T4_x_step_32=0;
         M_T4_x_step_33=0;
         M_T4_x_step_34=-1i*omega;
         M_T4_x_step_35=0;
         M_T4_x_step_36=0;
         
         M_T4_x_step_41=(C1*eta1*c_duct^2)/(1i*omega*Sp1);
         M_T4_x_step_42=(C1*eta1*c_duct^2)/(1i*omega*Sp1);
         M_T4_x_step_43=-((C1*eta1)/(1i*omega*Sp1)+(C2*eta2)/(1i*omega*Sp1)+(1i*omega)/(c_Liner^2));
         M_T4_x_step_44=0;  
         M_T4_x_step_45=C2*eta2/(1i*omega*Sp1);
         M_T4_x_step_46=0;
         
         M_T4_x_step_51=0;
         M_T4_x_step_52=0;
         M_T4_x_step_53=0;
         M_T4_x_step_54=0;
         M_T4_x_step_55=0;
         M_T4_x_step_56=-1i*omega;
         
         M_T4_x_step_61=0;
         M_T4_x_step_62=0;
         M_T4_x_step_63=C2*eta2/(1i*omega*Sp2);
         M_T4_x_step_64=0;  
         M_T4_x_step_65=-(C2*eta2/(1i*omega*Sp2)+(1i*omega)/(c_Liner^2));
         M_T4_x_step_66=0;
         
         M_x_step_and_half=[M_T4_x_step_11, M_T4_x_step_12, M_T4_x_step_13, M_T4_x_step_14, M_T4_x_step_15, M_T4_x_step_16;
                            M_T4_x_step_21, M_T4_x_step_22, M_T4_x_step_23, M_T4_x_step_24, M_T4_x_step_25, M_T4_x_step_26;
                            M_T4_x_step_31, M_T4_x_step_32, M_T4_x_step_33, M_T4_x_step_34, M_T4_x_step_35, M_T4_x_step_36;
                            M_T4_x_step_41, M_T4_x_step_42, M_T4_x_step_43, M_T4_x_step_44, M_T4_x_step_45, M_T4_x_step_46;
                            M_T4_x_step_51, M_T4_x_step_52, M_T4_x_step_53, M_T4_x_step_54, M_T4_x_step_55, M_T4_x_step_56;
                            M_T4_x_step_61, M_T4_x_step_62, M_T4_x_step_63, M_T4_x_step_64, M_T4_x_step_65, M_T4_x_step_66];
    
         M2_x_step=M_x_step_and_half*(I_T4+dx/2*M1_x_step);
         % Third step
         M3_x_step=M_x_step_and_half*(I_T4+dx/2*M2_x_step);
         % Forth step
         u_x_forthis_M = u_x(x_step+1);
         M_T4_x_step_11=-(1i*omega/(c_duct+u_x_forthis_M)+(C1*c_duct*eta1*rho_Liner)/(2*1i*omega*Sp*rho_duct));
         M_T4_x_step_12=-(C1*c_duct*eta1*rho_Liner)/(2*1i*omega*Sp*rho_duct);
         M_T4_x_step_13=(C1*eta1*rho_Liner)/(2*1i*omega*c_duct*Sp*rho_duct);
         M_T4_x_step_14=0;
         M_T4_x_step_15=0;
         M_T4_x_step_16=0;
         
         M_T4_x_step_21=(C1*c_duct*eta1*rho_Liner)/(2*1i*omega*Sp*rho_duct);
         M_T4_x_step_22=1i*omega/(c_duct-u_x_forthis_M)+(C1*c_duct*eta1*rho_Liner)/(2*1i*omega*Sp*rho_duct);
         M_T4_x_step_23=-(C1*eta1*rho_Liner)/(2*1i*omega*c_duct*Sp*rho_duct);
         M_T4_x_step_24=0;
         M_T4_x_step_25=0;
         M_T4_x_step_26=0;
         
         M_T4_x_step_31=0;
         M_T4_x_step_32=0;
         M_T4_x_step_33=0;
         M_T4_x_step_34=-1i*omega;
         M_T4_x_step_35=0;
         M_T4_x_step_36=0;
         
         M_T4_x_step_41=(C1*eta1*c_duct^2)/(1i*omega*Sp1);
         M_T4_x_step_42=(C1*eta1*c_duct^2)/(1i*omega*Sp1);
         M_T4_x_step_43=-((C1*eta1)/(1i*omega*Sp1)+(C2*eta2)/(1i*omega*Sp1)+(1i*omega)/(c_Liner^2));
         M_T4_x_step_44=0;  
         M_T4_x_step_45=C2*eta2/(1i*omega*Sp1);
         M_T4_x_step_46=0;
         
         M_T4_x_step_51=0;
         M_T4_x_step_52=0;
         M_T4_x_step_53=0;
         M_T4_x_step_54=0;
         M_T4_x_step_55=0;
         M_T4_x_step_56=-1i*omega;
         
         M_T4_x_step_61=0;
         M_T4_x_step_62=0;
         M_T4_x_step_63=C2*eta2/(1i*omega*Sp2);
         M_T4_x_step_64=0;  
         M_T4_x_step_65=-(C2*eta2/(1i*omega*Sp2)+(1i*omega)/(c_Liner^2));
         M_T4_x_step_66=0;
         
         M_x_step_next_step=[M_T4_x_step_11, M_T4_x_step_12, M_T4_x_step_13, M_T4_x_step_14, M_T4_x_step_15, M_T4_x_step_16;
                            M_T4_x_step_21, M_T4_x_step_22, M_T4_x_step_23, M_T4_x_step_24, M_T4_x_step_25, M_T4_x_step_26;
                            M_T4_x_step_31, M_T4_x_step_32, M_T4_x_step_33, M_T4_x_step_34, M_T4_x_step_35, M_T4_x_step_36;
                            M_T4_x_step_41, M_T4_x_step_42, M_T4_x_step_43, M_T4_x_step_44, M_T4_x_step_45, M_T4_x_step_46;
                            M_T4_x_step_51, M_T4_x_step_52, M_T4_x_step_53, M_T4_x_step_54, M_T4_x_step_55, M_T4_x_step_56;
                            M_T4_x_step_61, M_T4_x_step_62, M_T4_x_step_63, M_T4_x_step_64, M_T4_x_step_65, M_T4_x_step_66];
         M4_x_step=M_x_step_next_step*(I_T4+dx*M3_x_step);
    
         % Transfer matrix between parameters at x_step and (x_step+1)
         M_x_step_transf=I_T4+dx/6*(M1_x_step+2*M2_x_step+2*M3_x_step+M4_x_step);
         
         if x_step==1
             M_x_0_L(x_step).mat=M_x_0_L(x_step).mat*M_x_step_transf;
         else
             M_x_0_L(x_step).mat=M_x_0_L(x_step-1).mat*M_x_step_transf;
         end
         M_x_0_L_type4(x_step).mat=M_x_0_L(x_step).mat;
    
         end
         
         %% Link oscillations at the begin and any point across the Liner
         % Calculate B1(0) and B2(0)
         M_L=[M_x_0_L(x_step_all-1).mat(4,3), M_x_0_L(x_step_all-1).mat(4,5);
              M_x_0_L(x_step_all-1).mat(6,3), M_x_0_L(x_step_all-1).mat(6,5)];
         M_R=[-M_x_0_L(x_step_all-1).mat(4,1), -M_x_0_L(x_step_all-1).mat(4,2);
              -M_x_0_L(x_step_all-1).mat(6,1), -M_x_0_L(x_step_all-1).mat(6,2)];
         M_phi0_to_B0=M_L\M_R;
        
         for x_step=1:1:(x_step_all-1);
         M_x_0_L_type4_22(x_step).mat=[M_x_0_L(x_step).mat(1,1)+M_x_0_L(x_step).mat(1,3)*M_phi0_to_B0(1,1)+M_x_0_L(x_step).mat(1,5)*M_phi0_to_B0(2,1), M_x_0_L(x_step).mat(1,2)+M_x_0_L(x_step).mat(1,3)*M_phi0_to_B0(1,2)+M_x_0_L(x_step).mat(1,5)*M_phi0_to_B0(2,2);
                                       M_x_0_L(x_step).mat(2,1)+M_x_0_L(x_step).mat(2,3)*M_phi0_to_B0(1,1)+M_x_0_L(x_step).mat(2,5)*M_phi0_to_B0(2,1), M_x_0_L(x_step).mat(2,2)+M_x_0_L(x_step).mat(2,3)*M_phi0_to_B0(1,2)+M_x_0_L(x_step).mat(2,5)*M_phi0_to_B0(2,2)];
         A_to_phi_l(x_step).mat=[(1+M_duct_l)/c_duct^2,0;
                                 0, (1-M_duct_l)/c_duct^2];
         A_to_phi_r(x_step).mat=[(1+u_x(x_step+1)/c_duct)/c_duct^2,0;
                                 0, (1-u_x(x_step+1)/c_duct)/c_duct^2];
         Mtr_Liner_22=A_to_phi_r(x_step).mat\M_x_0_L_type4_22(x_step).mat*A_to_phi_l(x_step).mat;
         
         tau_c(x_step)=x_locations(x_step)/CI.TP.u_mean(1,section_Num);    % Entropy wave convecting time
         Osc_trans_Mtr_across_Liner(x_step).mat=[Mtr_Liner_22(1,1), Mtr_Liner_22(1,2), 0;
                                                 Mtr_Liner_22(2,1), Mtr_Liner_22(2,2), 0;
                                                 0                , 0,                 exp(-tau_c(x_step).*i_omega)];
         end
        
              %% End Type 4
              %%%%%%%%%%%%%%%%%%%
             end
             % Not probable to happen
    end
     %% End double liners model
    end
% Not probable to happen
end

%----------------------------------end-------------------------------------