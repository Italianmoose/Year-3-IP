function [x_resample,p,u]=Fcn_calculation_eigenmode_Linear(s_star)
% This function is used to plot the modeshape of slected mode, 
% Linear
global CI
[R1,R2]     = Fcn_boundary_condition(s_star);
Rs          = -0.5*CI.TP.M_mean(end)./(1 + 0.5*(CI.TP.gamma(end) - 1 ).*CI.TP.M_mean(end)).*CI.TP.rho_mean(1,end).*CI.TP.c_mean(1,end).^2./CI.TP.p_mean(1,end);
%
%
Te          = Fcn_TF_entropy_convection(s_star);
%--------------------------------
tau_plus    = CI.TP.tau_plus;
tau_minus   = CI.TP.tau_minus;
tau_c       = CI.TP.tau_c;
%
% -------------------------------------------------------------------------
A_minus(1)  = 1;
E(1)        = 0;
A_plus(1)   = R1.*A_minus(1);
Array(:,1)  = [A_plus(1),A_minus(1),E(1)]';
%
indexHA = 0;            % index of heat addition
indexHP = 0;            % index of heat perturbation
HR_Flag=1;
Liner_Flag=1;
% -------------------------------------------------------------------------
for ss = 1:CI.TP.numSection-1 
    D1 = diag([ exp(-s_star*tau_plus(ss)),...
                exp( s_star*tau_minus(ss)),...
                exp(-s_star*tau_c(ss))]);
    switch CI.CD.SectionIndex(ss+1)
        case 0
            CI.TPM.Z{ss}       = CI.TPM.BC{ss}*D1;
        case 10
            indexHA = indexHA + 1;
            B2b     = zeros(3);
            B2b(3,2)= 0;
            Bsum    = CI.TPM.B1{2,ss}*(CI.TPM.B2{1,ss}\CI.TPM.B1{1,ss}) + B2b;
            BC1     = Bsum*CI.TPM.C1;
            BC2     = CI.TPM.B2{2,ss}*CI.TPM.C2;
            CI.TPM.Z{ss}       = (BC2\BC1)*D1;
        case 11
            % linear flame transfer function
            indexHA = indexHA + 1;
            indexHP = indexHP + 1;
            FTF     = Fcn_flame_model(s_star,indexHP);
            B2b     = zeros(3);
            temp    = D1*Array(:,ss);
            uf      = abs((temp(1) - temp(2))./(CI.TP.c_mean(1,ss).*CI.TP.rho_mean(1,ss)));         % velocity before the flame
            B2b(3,2)= CI.TP.DeltaHr(indexHA)./CI.TP.c_mean(2,ss+1)./CI.TP.c_mean(1,ss)./CI.TP.Theta(ss).*FTF;
            Bsum    = CI.TPM.B1{2,ss}*(CI.TPM.B2{1,ss}\CI.TPM.B1{1,ss}) + B2b;
            BC1     = Bsum*CI.TPM.C1;
            BC2     = CI.TPM.B2{2,ss}*CI.TPM.C2;
            CI.TPM.Z{ss}       = (BC2\BC1)*D1;
        case 2
            HR_Num         =HR_Flag;
            i_omega        =s_star;                               
            w_0            =[Array(1,ss); Array(2,ss); Array(3,ss)]; % These values are useless in the linear HR model, but need to be defined for the HR function
            Array_LeftInt  =D1*w_0;
            Mtr_HR_33=Fcn_calculation_Oscillations_across_HR(ss,HR_Num,i_omega, Array_LeftInt(1),Array_LeftInt(2));
            CI.TPM.BC{ss}  =Mtr_HR_33; 
            CI.TPM.Z{ss}   =CI.TPM.BC{ss}*D1;
            HR_Flag  =HR_Flag+1;
        case 30
            % Refer to the transfer matrix calculated by an outer liner
            % wave calculation function
            CI.TPM.BC{ss}  =diag([1,...
                                  1,...
                                  1]);
            CI.TPM.Z{ss}   = CI.TPM.BC{ss}*D1;
        case 31
            Liner_Num      =Liner_Flag;
            i_omega        =s_star;                           % s is used here directory as the real part of it is assumed to be small and thus will not affect the liner performance much.
            Mtr_Liner_33   = Fcn_calculation_Oscillations_across_Liner(ss,Liner_Num,i_omega);
            CI.TPM.BC{ss}  =Mtr_Liner_33;
            D_Liner        =diag([ 1    1    exp(-tau_c(ss).*s_star)]);  % It is assumed that only acoustic waves are considered in the liner damping model. Entropy wave convects from left side to the right side.
            CI.TPM.Z{ss}   =CI.TPM.BC{ss}*D_Liner;
            Liner_Flag     =Liner_Flag+1;
    end
    Array(:,ss+1)    = CI.TPM.Z{ss}*Array(:,ss);
end
%
for k = 1:CI.TP.numSection-1
    A_plus(k+1)     = Array(1,k+1);
    A_minus(k+1)    = Array(2,k+1);
    E(k+1)          = Array(3,k+1);
end
%
for k=1:length(CI.CD.x_sample)-1
    x_resample(k,:)=linspace(CI.CD.x_sample(k),CI.CD.x_sample(k+1),51);    % resample of x-coordinate
end
% speed of sound and mean velocity
c_mean      = CI.TP.c_mean(1,:);
u_mean      = CI.TP.u_mean(1,:);
rho_mean    = CI.TP.rho_mean(1,:);
%     
for k = 1:length(CI.CD.x_sample)-1
    Liner_Flag=1;
    if CI.CD.SectionIndex(k)==30
        Liner_Num    =Liner_Flag;
        [Osc_trans_Mtr_across_Liner] = Fcn_calculation_Osc_trans_matrixes_across_Liner(k,Liner_Num,s_star);
        A_plus_k(1)  =A_plus(k);
        A_minus_k(1) =A_minus(k);
        E_k(1)       =0;
        for step_k=2:1:51
            wave_k=Osc_trans_Mtr_across_Liner(step_k-1).mat*[A_plus_k(1);A_minus_k(1);E_k(1)];
            A_plus_k(step_k)   =wave_k(1);
            A_minus_k(step_k)  =wave_k(2);
        end   
        p(k,:)      =A_plus_k+A_minus_k;
        u(k,:)      =(A_plus_k-A_minus_k)./rho_mean(k)./c_mean(k);
        Liner_Flag  =Liner_Flag+1;
    else
    kw1_plus(k)  = s_star./(c_mean(k)+u_mean(k));
    kw1_minus(k) = s_star./(c_mean(k)-u_mean(k));
    p(k,:) =    A_plus(k).*exp(-kw1_plus(k).*(x_resample(k,:)-x_resample(k,1)))+...
                A_minus(k).*exp( kw1_minus(k).*(x_resample(k,:)-x_resample(k,1)));
    u(k,:) =    (A_plus(k).*exp(-kw1_plus(k).*(x_resample(k,:)-x_resample(k,1)))-...
                A_minus(k).*exp( kw1_minus(k).*(x_resample(k,:)-x_resample(k,1))))./rho_mean(k)./c_mean(k);
    end
end


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