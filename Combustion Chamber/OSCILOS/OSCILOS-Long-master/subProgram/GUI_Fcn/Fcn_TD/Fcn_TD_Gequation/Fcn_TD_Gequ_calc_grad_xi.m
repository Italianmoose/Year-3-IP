function [ grad_xi ] = Fcn_TD_Gequ_calc_grad_xi( u_gutter,SU,xi,y_vec )
%Simple backward euler first oder graidnt method. The boundary condition is
%applied such that the flame remains anchored at y(1) if SU < u_gutter. If
%this is not the case, then the gradinet at y(1) is set to zero

% Compute BC
BC_grad = zeros(size(u_gutter));
for runner = 1:length(u_gutter)
    if u_gutter(runner) >= SU(runner)
        BC_grad(runner) = (xi(runner,1) >= 0 ) .* sqrt(( u_gutter(runner)./SU(runner)).^2 - 1); % Implicit test is to determine whether or not we have reached the anchoring point
    else
        BC_grad(runner) = 0;
    end
end

prod_vec = ones(length(u_gutter),1);
grad_xi = cat(2,BC_grad, diff(xi,1,2)./(prod_vec * diff(y_vec)));

end

