function [ dxidt ] = Fcn_TD_Gequ_calc_dxidt( u_gutter,SU,xi,y_vec )
%This computes the rate of change dxidt, and applies boundary conditions

% Compute the gradient of the flame
grad_xi = Fcn_TD_Gequ_calc_grad_xi(u_gutter,SU,xi,y_vec);

prod_vec = ones(1,size(grad_xi,2));
dxidt = u_gutter * prod_vec - SU * sqrt(1 + grad_xi.^2);

end

