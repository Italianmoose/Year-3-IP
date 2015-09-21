function x = Fcn_saturation(x,a)
% This function is an abrupt saturation filter, x is the orginal signal,
% a is the bound and is a non-negative real value
x(abs(x)>a) = sign(x(abs(x)>a)).*a;
% --------------------------end--------------------------------------------
