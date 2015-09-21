function [uRatio,Lf] = Fcn_GUI_INI_FM_Nonlinear_model_J_Li_A_Morgas(alpha,beta,DuRatio,uRatioMax)
%
% nonlinear model: from J.Li and A.Morgans
uRatio  = DuRatio:DuRatio:uRatioMax; 
slope   = 1./(1+(uRatio + alpha).^beta);
qRatio  = zeros(1,length(uRatio));
Lf      = zeros(1,length(uRatio));
for ss=1:length(uRatio)
    qRatio(ss) = sum(slope(1:ss)).*DuRatio;                     % Intergrated q_ratio
    Lf(ss)=qRatio(ss)./uRatio(ss);                              % q_ratio/y_ratio
end
%
% ----------------------------end------------------------------------------

