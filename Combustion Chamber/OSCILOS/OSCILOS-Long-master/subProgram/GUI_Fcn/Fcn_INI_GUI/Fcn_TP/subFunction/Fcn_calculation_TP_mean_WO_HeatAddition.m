function [p_mean2,rho_mean2,u_mean2] =...
    Fcn_calculation_TP_mean_WO_HeatAddition(p_mean1,rho_mean1,u_mean1,Theta,gamma)
% This function is used to calculate the mean thermoproperties changing
% across the interface between two sections with different section surface
% area. 
% The heat addition is not accounted for in this program
% Theta=S2/S1;
% gamma is considered constant before and after the flame separately.
global MF
MF.x1 = [p_mean1,rho_mean1,u_mean1];
MF.Theta = Theta;
MF.gamma = gamma;
assignin('base','MF',MF);                   % save the current information to the workspace
%--------------------------------------------------------------------------
options = optimset('Display','off');        % the calculation results by fsolve will not be shown in the workspace
x0 = MF.x1;                                 % Make a starting guess at the solution
[x2,fval] = fsolve(@myfun,x0,options);      % Call solver
%
p_mean2     = x2(1);               
rho_mean2   = x2(2);
u_mean2     = x2(3);
clear MF

function F = myfun(x2)
global MF
% p2   = x2(1);
% rho2 = x2(2);
% u2   = x2(3);
x1      = MF.x1;
Theta   = MF.Theta;
gamma   = MF.gamma;

F(1) = Theta*x2(2)*x2(3) - x1(2)*x1(3);
F(3) = gamma/(gamma-1)*x2(1)/x2(2)+0.5*x2(3)^2 - (gamma/(gamma-1)*x1(1)/x1(2)+0.5*x1(3)^2);

if Theta>=1              % Area increasing
    F(2) = x2(1)+x2(2)*x2(3)^2 - (x1(1)+x1(2)*x1(3)^2/Theta);
elseif Theta<1       % Area decreasing
    F(2) = x2(1)*x1(2)^gamma - x1(1)*x2(2)^gamma;
end
%----------------------------------end-------------------------------------