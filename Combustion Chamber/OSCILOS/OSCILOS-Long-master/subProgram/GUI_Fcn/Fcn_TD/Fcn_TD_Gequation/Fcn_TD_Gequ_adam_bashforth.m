function [ new_vector ] = Fcn_TD_Gequ_adam_bashforth(old_vector, bash_data, dt, IT)
% Compute y of the next iteration using data from the previous iteration
% set the coefficients for 3 step bashforth

if IT >= 3
    coefs = [23/12 ; -4/3 ; 5/12];
else
    if IT == 2
        coefs = [3/2 ; -1/2 ; 0];
    else
        coefs = [1 ; 0 ; 0];
    end
end

bash_sum = 0;
for runner = 1:size(bash_data,3)
    bash_sum = bash_sum + bash_data(:,:,runner) .* coefs(runner);
end

new_vector = Fcn_TD_Gequ_forward_euler( old_vector, dt , bash_sum);

end

