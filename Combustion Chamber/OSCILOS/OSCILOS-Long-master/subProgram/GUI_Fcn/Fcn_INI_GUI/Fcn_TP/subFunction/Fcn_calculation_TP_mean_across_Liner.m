function [u_mean2]...
        = Fcn_calculation_TP_mean_across_Liner(Liner_Num, section_Num)
% This function is developed to calculate mean flow velocity immediately
% after the lined duct.

global CI

if CI.CD.Liner_config(Liner_Num,4)==0 % Single liner is used
                T_l       =CI.CD.Liner_config(Liner_Num,7);                       % Mean temperature in the liner system
                [c_mean,temp1,temp2,Cp]   = Fcn_calculation_c_q_air(T_l);         % Calculate Cp in the liner system based on the given mean temperature
                gamma     = Cp/(Cp - CI.R_air);                                   % Using liner temperature to calculate Gamma outside the liner
                c_l       =sqrt(gamma*CI.R_air*T_l);                              % Sound speed in the liner system
                v_1_l     =CI.CD.Liner_config(Liner_Num,10)*c_l;                  % Mean hole velocity in the line holes
                P_l       =CI.TP.p_mean(1,section_Num)/(1-1/2/CI.R_air/T_l*(v_1_l/0.64)^2);% Calculate mean pressure in the liner system. The discharge coefficient is assumed to be 0.64 according to Falkovich (Fluid mechanics, 2011)
                rho_l     =P_l/CI.R_air/T_l;                                      % Mean density in the liner system
                sigma_1_l =pi*CI.CD.Liner_config(Liner_Num,8)^2/CI.CD.Liner_config(Liner_Num,9)^2; % Open area ratio of the liner
                
                u_mean2   =CI.TP.u_mean(1,section_Num)+rho_l*(2*pi*CI.CD.r_sample(section_Num))*(sigma_1_l*v_1_l)*CI.CD.Liner_config(Liner_Num,3)/CI.TP.rho_mean(1,section_Num)/pi/CI.CD.r_sample(section_Num)^2; % Mean velocity after the liner
else if CI.CD.Liner_config(Liner_Num,4)==1 % Double liners is used 
 
                    T_l       =CI.CD.Liner_config(Liner_Num,7);                       % Mean temperature in the liner system
                    [c_mean,temp1,temp2,Cp]   = Fcn_calculation_c_q_air(T_l);         % Calculate Cp in the liner system based on the given mean temperature
                    gamma     = Cp/(Cp - CI.R_air);                                   % Using liner temperature to calculate Gamma outside the liner
                    c_l       =sqrt(gamma*CI.R_air*T_l);                              % Sound speed in the liner system
                    v_1_l     =CI.CD.Liner_config(Liner_Num,14)*c_l;                  % Mean hole velocity in the first line holes
                    P_l       =CI.TP.p_mean(1,section_Num)/(1-1/2/CI.R_air/T_l*(v_1_l/0.64)^2);% Calculate mean pressure in the liner system. The discharge coefficient is assumed to be 0.64 according to Falkovich (Fluid mechanics, 2011)
                    rho_l     =P_l/CI.R_air/T_l;                                      % Mean density in the liner system
                    sigma_1_l =pi*CI.CD.Liner_config(Liner_Num,12)^2/CI.CD.Liner_config(Liner_Num,13)^2; % Open area ratio of the first liner
                
                    u_mean2   =CI.TP.u_mean(1,section_Num)+rho_l*2*pi*CI.CD.r_sample(section_Num)*sigma_1_l*v_1_l*CI.CD.Liner_config(Liner_Num,3)/CI.TP.rho_mean(1,section_Num)/pi/CI.CD.r_sample(section_Num)^2; % Mean velocity after the liner 
     end
                 % Unlikely to happen
end
%--------------------------------end---------------------------------------%


      
