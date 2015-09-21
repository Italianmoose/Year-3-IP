function F = Fcn_DetEqn(s)
global CI
% boundary condition
[R1,R2]     = Fcn_boundary_condition(s);
Rs          = -0.5*CI.TP.M_mean(end)./(1 + 0.5*(CI.TP.gamma(end) - 1 ).*CI.TP.M_mean(end)).*CI.TP.rho_mean(1,end).*CI.TP.c_mean(1,end).^2./CI.TP.p_mean(1,end);;
%
%
Te = Fcn_TF_entropy_convection(s);
%--------------------------------
tau_plus    = CI.TP.tau_plus;
tau_minus   = CI.TP.tau_minus;
tau_c       = CI.TP.tau_c;
%--------------------------------
G = eye(3);
for ss = 1:CI.TP.numSection-1 
    D1 = diag([ exp(-s*tau_plus(ss)),...
                exp( s*tau_minus(ss)),...
                exp(-s*tau_c(ss))]);
    switch CI.CD.SectionIndex(ss+1)
        case 0
            Z       = CI.TPM.BC{ss}*D1;
        case 10
            B2b     = zeros(3);
            B2b(3,2)= 0;
            Bsum    = CI.TPM.B1{2,ss}*(CI.TPM.B2{1,ss}\CI.TPM.B1{1,ss}) + B2b;
            BC1     = Bsum*CI.TPM.C1;
            BC2     = CI.TPM.B2{2,ss}*CI.TPM.C2;
            Z       = (BC2\BC1)*D1;
        case 11
            % Flame model
            FDF     = Fcn_flame_model(s);
            B2b     = zeros(3);
            B2b(3,2)= CI.TP.DeltaHr./CI.TP.c_mean(2,ss+1)./CI.TP.c_mean(1,ss)./CI.TP.Theta(ss).*FDF;
            Bsum    = CI.TPM.B1{2,ss}*(CI.TPM.B2{1,ss}\CI.TPM.B1{1,ss}) + B2b;
            BC1     = Bsum*CI.TPM.C1;
            BC2     = CI.TPM.B2{2,ss}*CI.TPM.C2;
            Z       = (BC2\BC1)*D1;
    end
    G = Z*G;
end
%
A1_minus        = 1;
A1_plus         = R1.*A1_minus;
E1              = 0;
Array_LeftBD    = [A1_plus, A1_minus, E1]';
%
D1End           = diag([    exp(-s*tau_plus(end)),...
                            exp( s*tau_minus(end)),...
                            Te.*exp(-s*tau_c(end))]);
%
Array_RightBD   = D1End*G*Array_LeftBD;
AN_plus         = Array_RightBD(1);
AN_minus        = Array_RightBD(2);
EN_plus         = Array_RightBD(3);

F = (R2.*AN_plus + Rs.*EN_plus) - AN_minus;  

%
%----------------------Pressure Reflection coefficients -------------------
%
function [R1,R2] = Fcn_boundary_condition(s)
global CI
R1      = polyval(CI.BC.num1,s)./polyval(CI.BC.den1,s).*exp(-CI.BC.tau_d1.*s);
R2      = polyval(CI.BC.num2,s)./polyval(CI.BC.den2,s).*exp(-CI.BC.tau_d2.*s);
%
%----------------------Entropy convection transfer function ---------------
%
function Te = Fcn_TF_entropy_convection(s)
global CI
tau     = CI.BC.ET.Dispersion.Delta_tauCs;
k       = CI.BC.ET.Dissipation.k;
switch CI.BC.ET.pop_type_model
    case 1
        Te = k;
    case 2
        Te = k.*exp((tau.*s).^2./4);
    case 3
        if tau == 0
            tau = eps;
        end
%         Te = k.*(exp(tau*s) - exp(-tau*s))./(2*tau);
        Te = k.*sinc(tau*s./pi);
end                        
%
% ----------------------Flame transfer function ---------------------------
%
function F = Fcn_flame_model(s)
global FDF
global CI
F = polyval(FDF.num,s)./polyval(FDF.den,s).*exp(-s.*FDF.tauf);
if CI.indexFM == 0
    switch CI.FM.NL.style
        case 2
            uRatio = FDF.uRatio;
            if FDF.uRatio == 0
                uRatio = eps;
            end
        qRatioLinear = abs(F).*uRatio; 
        Lf = interp1(  CI.FM.NL.Model2.qRatioLinear,...
                       CI.FM.NL.Model2.Lf,...
                       qRatioLinear,'linear','extrap');                                    
        F = F.*Lf; 
    %     otherwise
    end
end
clear FDF
%
% -----------------------------end-----------------------------------------
