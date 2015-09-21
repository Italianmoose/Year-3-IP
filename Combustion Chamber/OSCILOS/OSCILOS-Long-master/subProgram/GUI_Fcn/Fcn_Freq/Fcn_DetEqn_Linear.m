function F = Fcn_DetEqn_Linear(s)
% This function only accounts for linear cases
% 
global CI
% boundary condition
[R1,R2]     = Fcn_boundary_condition(s);
Rs          = -0.5*CI.TP.M_mean(end)./(1 + 0.5*(CI.TP.gamma(end) - 1 ).*CI.TP.M_mean(end)).*CI.TP.rho_mean(1,end).*CI.TP.c_mean(1,end).^2./CI.TP.p_mean(1,end);
%
%
Te = Fcn_TF_entropy_convection(s);
%--------------------------------
tau_plus    = CI.TP.tau_plus;
tau_minus   = CI.TP.tau_minus;
tau_c       = CI.TP.tau_c;
%--------------------------------
indexHA = 0;            % index of heat addition
indexHP = 0;            % index of heat perturbation
HR_Flag=1;
Liner_Flag=1;
% Initial wave amplitudes
A_0_n=1;
A_0_p=A_0_n*R1;
A_E_0=0;
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
            FTF     = Fcn_flame_model(s,indexHP);
            B2b     = zeros(3);
            % in case there are two heat addition, but the first one is a
            % steady one and the second one is unsteady  ????
            % this problem is solved in GUI_FREQ_EigCal_pannel_appearance
            B2b(3,2)= CI.TP.DeltaHr(indexHA)./CI.TP.c_mean(2,ss+1)./CI.TP.c_mean(1,ss)./CI.TP.Theta(ss).*FTF;
            Bsum    = CI.TPM.B1{2,ss}*(CI.TPM.B2{1,ss}\CI.TPM.B1{1,ss}) + B2b;
            BC1     = Bsum*CI.TPM.C1;
            BC2     = CI.TPM.B2{2,ss}*CI.TPM.C2;
            Z       = (BC2\BC1)*D1;
        case 2
            HR_Num   =HR_Flag;
            i_omega  =s;                               %  s is used here directory as the real part of it is assumed to be small and thus will not affect the HR performance much.
            % Calculate pressure wave values at the inlet of the HR (not useful for linear model but needed for nonlinear model)
            w_1=D1*G*[A_0_p; A_0_n; A_E_0];
            A_1_p=w_1(1);
            A_1_n=w_1(2);
            Mtr_HR_33=Fcn_calculation_Oscillations_across_HR(ss,HR_Num,i_omega, A_1_p,A_1_n);
            E=[1,0,0;0,1,0;0,0,0];
            Z        =E*Mtr_HR_33*D1; 
            HR_Flag  =HR_Flag+1;
        case 30
            Z=eye(3)*D1;
            
        case 31
            Liner_Num    =Liner_Flag;
            i_omega      =s;                           % s is used here directory as the real part of it is assumed to be small and thus will not affect the liner performance much.
            Mtr_Liner_33 = Fcn_calculation_Oscillations_across_Liner(ss,Liner_Num,i_omega);
            D_Liner      =diag([ 1,...
                                 1,...
                                 exp(-s*tau_c(ss))]);  % It is assumed that only acoustic waves are considered in the liner damping model. Entropy wave convects from left side to the right side.
            Z            =Mtr_Liner_33*D_Liner;
            Liner_Flag   =Liner_Flag+1;
    end
    G = Z*G;
end
%
A1_minus        = 1000;
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
%         if tau == 0
%             tau = eps;
%         end
%         Te = k.*(exp(tau*s) - exp(-tau*s))./(2*tau);
        Te = k.*sinc(tau*s./pi);
end                        
%
% ----------------------linear Flame transfer function --------------------
%         
function F = Fcn_flame_model(s,indexHP)
global CI
HP      = CI.FM.HP{indexHP};
num     = HP.FTF.num;
den     = HP.FTF.den;
tauf    = HP.FTF.tauf;
F       = polyval(num,s)./polyval(den,s).*exp(-s.*tauf);          
%
% -----------------------------end-----------------------------------------
