function [ new_field ] = Fcn_TD_Gequ_forward_euler( old_field,dt,f )
%Forward Euler Routine goes here.

new_field = old_field + dt * f;

end

