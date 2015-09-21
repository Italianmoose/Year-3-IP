function [ xi ] = Fcn_TD_Gequ_RK(SU,xi,y_vec,dt,current_time,U1,area_ratio,uratio,u_gutter )
%This function uses the Runge kutta method to compute a new value of the
%flame
      
% Stage 1
[K1] = Fcn_TD_Gequ_calc_dxidt ( u_gutter,SU,xi,y_vec);
XI1 = Fcn_TD_Gequ_forward_euler( xi,0.5 * dt,K1 );
current_time = current_time + 0.5*dt;

% Stage 2
u_gutter = Fcn_TD_Gequ_calc_ugutter( U1,area_ratio,uratio,current_time);
[K2] = Fcn_TD_Gequ_calc_dxidt (  u_gutter,SU,XI1,y_vec);
XI2 = Fcn_TD_Gequ_forward_euler( xi,0.5 * dt,K2);

% Stage 3
u_gutter = Fcn_TD_Gequ_calc_ugutter( U1,area_ratio,uratio,current_time);
[K3] = Fcn_TD_Gequ_calc_dxidt (  u_gutter,SU,XI2,y_vec);
XI3 = Fcn_TD_Gequ_forward_euler( xi,dt,K3);
current_time = current_time + 0.5 * dt;

% Stage 4
u_gutter = Fcn_TD_Gequ_calc_ugutter( U1,area_ratio,uratio,current_time);
[K4] = Fcn_TD_Gequ_calc_dxidt (  u_gutter,SU,XI3,y_vec);

% Run runge kutta iteration
xi = xi + 1/6 * dt * (K1+ 2*K2 +  2*K3 + K4);

end

