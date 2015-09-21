function [ q_ratio,xi,bashforth_data ] = Fcn_TD_Gequ_interface( SU, xi, y_vec, dt, current_time, U1, area_ratio, uratio,Q_mean,delta_h,rho1,bashforth_data,IT,time_integration )
%This is the G-equation interface with OSCILOS. It asks in all of the
%required variables, and gives out the heat release ratio to be used in
%the acoustics calculations

% The G_equation is a three step operation, First find the velocity right before the flame
u_gutter = Fcn_TD_Gequ_calc_ugutter( U1,area_ratio,uratio,current_time);

% Then compute the new flame shape
% For this we have a choice of different numerical methods
switch time_integration 
    case(1)
        [ dxidt ] = Fcn_TD_Gequ_calc_dxidt( u_gutter,SU,xi,y_vec );
        [ xi ] = Fcn_TD_Gequ_forward_euler( xi,dt,dxidt );
    case(2)
        [ xi ] = Fcn_TD_Gequ_RK(SU,xi,y_vec,dt,current_time,U1,area_ratio,uratio,u_gutter );
    case (3) %multistep method, adam bashforth
        [ dxidt ] = Fcn_TD_Gequ_calc_dxidt( u_gutter,SU,xi,y_vec );
        % bashforth data is a three dimentional matrix (lines are different flames, columns are different r, last dimention are different iterations 
        bashforth_data(:,:,2:3) = bashforth_data(:,:,1:2);
        bashforth_data(:,:,1) =  dxidt;
        xi = Fcn_TD_Gequ_adam_bashforth( xi, bashforth_data, dt, IT);
end

% Then compute the heat release from the area
[ Q_fluct ] = Fcn_TD_Gequ_calc_Q_fluct( y_vec,SU,xi,u_gutter, Q_mean, delta_h,rho1,y_vec);

% and finally compute the Q_ratio
q_ratio = Q_fluct./Q_mean; % Q_mean can be a vector

end

