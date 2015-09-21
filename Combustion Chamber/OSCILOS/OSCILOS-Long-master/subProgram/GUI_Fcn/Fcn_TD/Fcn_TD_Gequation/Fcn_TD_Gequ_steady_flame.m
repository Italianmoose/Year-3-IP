function [ xi_steady ] = Fcn_TD_Gequ_steady_flame( U_gutter_steady,SU,y_vec )
%This function computes the one sided steady flame shape of the flame using the vale
%of the laminar burning velocity SU and the mean flow. U1

slope = sqrt((U_gutter_steady./SU).^2 - 1); % If there are multple values of U_gutter_steady and SU, this is a vector of size(1,n)
xi_steady = zeros(size(y_vec)); % Matrix of size (num flames, num_points)
for runner = 1:length(slope)
    xi_steady(runner,:) = slope(runner) * y_vec(runner,:) - y_vec(runner,1) * slope(runner);
end


end

