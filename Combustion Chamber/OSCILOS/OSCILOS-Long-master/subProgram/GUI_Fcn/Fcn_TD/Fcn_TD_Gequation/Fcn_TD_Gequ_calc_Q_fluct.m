function [ Q_fluct ] = Fcn_TD_Gequ_calc_Q_fluct( y_vect,SU,xi,u_gutter, Q_mean, delta_h,rho1,y_vec)
%Compute the fluctuating heat release from the flame shape

grad_xi = Fcn_TD_Gequ_calc_grad_xi( u_gutter,SU,xi,y_vec );

flame_area = zeros(size(grad_xi,1),1);
diff_y_vec = diff(y_vec);
% The formula for the area of the total flame is defined as int_{ra}^{rb}
% (2 pi r sqrt(1+dxidr^2))
for runner = 1:length(y_vect)-1 % Integral using trapezoidal method
    flame_area = flame_area +  diff_y_vec(runner) .* ...
        pi *(y_vect(runner).*sqrt(1.0+grad_xi(:,runner).^2) + y_vect(runner+1)*sqrt(1.0+grad_xi(:,runner+1).^2))/2;
end

flame_area = flame_area * 2; % multiply by 2 so that the flame area for top and bottom half of the duct are accounted for

% The fomula for the heat release is defined as Q = flame_area *
% laminar_burning_vel * rho1 * delta_h  (where delta_h is the heat
% release per unit mass of fuel mix)
Q_total = flame_area .*SU.*rho1.*delta_h;
Q_fluct = Q_total - Q_mean; % Q_mean and delta_h can be a vector

end

