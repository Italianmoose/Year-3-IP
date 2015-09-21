function F = Fcn_DetEqn_nonlinear_dampers(s)
global CI
[R1,R2]     = Fcn_boundary_condition(s);
Rs          = -0.5*CI.TP.M_mean(end)./(1 + 0.5*(CI.TP.gamma(end) - 1 ).*CI.TP.M_mean(end));
%
Te = Fcn_TF_entropy_convection(s);
%--------------------------------
tau_plus    = CI.TP.tau_plus;
tau_minus   = CI.TP.tau_minus;
tau_c       = CI.TP.tau_c;
% -------------------------------------------------------------------------
A_minus(1)=CI.EIG.PR.A1minus_this_step;
   if abs(A_minus(1)) < 0.1
      A_minus(1) = 0.1;
   end
A_plus(1)         = R1.*A_minus(1);
E(1)              = 0;  
Array(:,1)  = [A_plus(1); A_minus(1); E(1)];
%
indexHA = 0;            % index of heat addition
indexHP = 0;            % index of heat perturbation
HR_Flag=1;
Liner_Flag=1;
% -------------------------------------------------------------------------
G = eye(3);
for ss = 1:CI.TP.numSection-1 
    D1 = diag([ exp(-s*tau_plus(ss)),...
                exp( s*tau_minus(ss)),...
                exp(-s*tau_c(ss))]);
    switch CI.CD.SectionIndex(ss+1)
        case 0
            Z       = CI.TPM.BC{ss}*D1;
        case 10
            indexHA = indexHA + 1;
            B2b     = zeros(3);
            B2b(3,2)= 0;
            Bsum    = CI.TPM.B1{2,ss}*(CI.TPM.B2{1,ss}\CI.TPM.B1{1,ss}) + B2b;
            BC1     = Bsum*CI.TPM.C1;
            BC2     = CI.TPM.B2{2,ss}*CI.TPM.C2;
            Z       = (BC2\BC1)*D1;
        case 11
            indexHA = indexHA + 1;
            indexHP = indexHP + 1;
            B2b     = zeros(3);
            Array_LeftInt    = D1*G*Array(:,1);
            uRatio  = abs((Array_LeftInt(1) - Array_LeftInt(2))./(CI.TP.c_mean(1,ss).*CI.TP.rho_mean(1,ss).*CI.TP.u_mean(1,ss)));         % velocity ratio before the flame
            if CI.FM.indexFM(indexHP)==1
                 FTF= Fcn_linear_flame_model(s,indexHP);
             else
                 FTF= Fcn_flame_model_nonlinear_Damper(s,uRatio,indexHP);
            end
            B2b(3,2)= CI.TP.DeltaHr(indexHA)./CI.TP.c_mean(2,ss+1)./CI.TP.c_mean(1,ss)./CI.TP.Theta(ss).*FTF;
            Bsum    = CI.TPM.B1{2,ss}*(CI.TPM.B2{1,ss}\CI.TPM.B1{1,ss}) + B2b;
            BC1     = Bsum*CI.TPM.C1;
            BC2     = CI.TPM.B2{2,ss}*CI.TPM.C2;
            Z       = (BC2\BC1)*D1;
        case 2
            HR_Num         =HR_Flag;
            i_omega        =s;                               
            Array_LeftInt  =D1*G*Array(:,1); % These values are useless in the linear HR model, but need to be defined for the HR function
            Mtr_HR_33=Fcn_calculation_Oscillations_across_HR(ss,HR_Num,i_omega, Array_LeftInt(1),Array_LeftInt(2));
            CI.TPM.BC{ss}  =Mtr_HR_33; 
            Z   =CI.TPM.BC{ss}*D1;
            HR_Flag  =HR_Flag+1;
        case 30
            % Refer to the transfer matrix calculated by an outer liner
            % wave calculation function
            CI.TPM.BC{ss}  =diag([1,...
                                  1,...
                                  1]);
            Z   = CI.TPM.BC{ss}*D1;
        case 31
            Liner_Num      =Liner_Flag;
            i_omega        =s;                           % s is used here directory as the real part of it is assumed to be small and thus will not affect the liner performance much.
            Mtr_Liner_33   = Fcn_calculation_Oscillations_across_Liner(ss,Liner_Num,i_omega);
            CI.TPM.BC{ss}  =Mtr_Liner_33;
            D_Liner        =diag([ 1    1    exp(-tau_c(ss).*s)]);  % It is assumed that only acoustic waves are considered in the liner damping model. Entropy wave convects from left side to the right side.
            Z   =CI.TPM.BC{ss}*D_Liner;
            Liner_Flag     =Liner_Flag+1;
            
    end
    G = Z*G;
end
%
%
D1End           = diag([    exp(-s*tau_plus(end)),...
                            exp( s*tau_minus(end)),...
                            Te.*exp(-s*tau_c(end))]);
%
Array_RightBD   = D1End*G*Array(:,1);
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
        Te = 0;
    case 2
        Te = k.*exp((tau.*s).^2./4);
    case 3
        if tau == 0
            tau = eps;
        end
        Te = k.*(exp(tau*s) - exp(-tau*s))./(2*tau);
end                        
%

% ----------------------linear Flame transfer function --------------------
%         
function F = Fcn_linear_flame_model(s,indexHP)
global CI
HP      = CI.FM.HP{indexHP};
num     = HP.FTF.num;
den     = HP.FTF.den;
tauf    = HP.FTF.tauf;
F       = polyval(num,s)./polyval(den,s).*exp(-s.*tauf);          
% %
% % ---------------------- Nonlinear flame transfer function ----------------
% %
% function F = Fcn_nonlinear_flame_model(s)
% global FDF
% global CI
% global HP
% F = polyval(FDF.num,s)./polyval(FDF.den,s).*exp(-s.*FDF.tauf);
% if  CI.EIG.APP_style == 21    % from model
%     switch HP.NL.style
%         case 2
%             uRatio = FDF.uRatio;
%             if FDF.uRatio == 0
%                 uRatio = eps;
%             end
%         qRatioLinear = abs(F).*uRatio; 
%         Lf = interp1(  HP.NL.Model2.qRatioLinear,...
%                        HP.NL.Model2.Lf,...
%                        qRatioLinear,'linear','extrap');                                    
%         F = F.*Lf; 
%     %     otherwise
%     end
% end
% clear FDF
%----------------------- Nonlinear flame transfer function ----------------
function F = Fcn_flame_model_nonlinear_Damper(s,uRatio,indexHP)
global CI
HP      = CI.FM.HP{indexHP}; % This can be used to distinguish different Nonlinear 
                             % flame models in the system 
num     = HP.FTF.num;
den     = HP.FTF.den;
tauf    = HP.FTF.tauf;
F       = polyval(num,s)./polyval(den,s).*exp(-s.*tauf); 

% Nonlinear flame models can be added below
switch HP.NL.style
    case 1 
    case 2
        if uRatio > 1                                
            uRatio = 1;
        end
        qRatioLinear = abs(F).*uRatio; 
        Lf = interp1(  HP.NL.Model2.qRatioLinear,...
                       HP.NL.Model2.Lf,...
                       qRatioLinear,'linear','extrap');                                    
        F = F.*Lf; 
    case 3 % Should be checked carefully!
        if uRatio > 1                                
            uRatio = 1;
        end
        Lf              = interp1(  HP.NL.Model3.uRatio,...
                                    HP.NL.Model3.Lf,...
                                    uRatio,'linear','extrap');
        taufNSp         = HP.NL.Model3.taufN.*(1-Lf);
        F = F.*Lf.*exp(-s.*taufNSp);
end
% -----------------------------end-----------------------------------------
