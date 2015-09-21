function [ u_gutter ] = Fcn_TD_Gequ_calc_ugutter( U1,area_ratio,uratio,current_time)
%Computes the value of u_gutter from uratio, U_upstream and aratio

% Find the values of uratio_now
%uratio_now = uratio(1); % should depend on current time
uratio_now = uratio'; % Needs to be in format (n,1)

%Comupute u_gutter
u_gutter = U1 .*(1+uratio_now) ./ area_ratio;

end

