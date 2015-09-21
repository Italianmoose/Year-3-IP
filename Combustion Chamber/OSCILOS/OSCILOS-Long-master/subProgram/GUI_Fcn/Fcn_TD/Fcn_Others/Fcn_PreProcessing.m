function Fcn_PreProcessing
% This function is used to calculate some mean matrix coefficients
global CI
%
Fcn_calculation_Matrix_elements
%
x_diff = diff(CI.CD.x_sample);
CI.TP.tau_plus      = x_diff./(CI.TP.c_mean(1,:)+CI.TP.u_mean(1,:));
CI.TP.tau_minus     = x_diff./(CI.TP.c_mean(1,:)-CI.TP.u_mean(1,:));
CI.TP.tau_c         = x_diff./ CI.TP.u_mean(1,:);       % convection time delay
assignin('base','CI',CI)
%
% -------------------------------------------------------------------------
%
function Fcn_calculation_Matrix_elements
%
global CI
CI.TPM.C1 =  [      1   1   0;...
                    1  -1   0;...
                    1   1   0];
CI.TPM.C2 =  [      1   1   0;...
                    1  -1   0;...
                    1   1  -1];
%
% -------------------------------------
%
indexHA = 0;
for ss = 1:CI.TP.numSection-1 
    switch CI.CD.SectionIndex(ss+1)
        case 0
            [B1,B2] = Fcn_Matrix_calculation_WO_Addition_effect(ss,1,1);
            CI.TPM.B1{1,ss}         = B1;               % first step
            CI.TPM.B2{1,ss}         = B2;
            CI.TPM.B1{2,ss}         = eye(3);           % second step
            CI.TPM.B2{2,ss}         = eye(3);
            CI.TPM.BC{ss}           = (B2*CI.TPM.C2)\(B1*CI.TPM.C1);
        case {10,11} 
            indexHA = indexHA + 1;
            [B1,B2] = Fcn_Matrix_calculation_WO_Addition_effect(ss,1,2);
            CI.TPM.B1{1,ss}         = B1;               % first step
            CI.TPM.B2{1,ss}         = B2;
            [B1a,B2]= Fcn_Matrix_calculation_W_HA(ss,indexHA);
            CI.TPM.B1{2,ss}         = B1a;              % second step
            CI.TPM.B2{2,ss}         = B2;
            CI.TPM.BC{ss}           = eye(3);           % assign an identy matrix, cannot be directly used
    end
end
assignin('base','CI',CI)
%
% -------------------------------------------------------------------------
%
function [B1a,B2] = Fcn_Matrix_calculation_W_HA(ss,indexHA)
global CI
cRatio  = CI.TP.c_mean(1,ss+1)./CI.TP.c_mean(2,ss+1);
M1      = CI.TP.M_mean(2,ss+1);
M2      = CI.TP.M_mean(1,ss+1);
gamma1  = CI.TP.gamma(2,ss+1);
gamma2  = CI.TP.gamma(1,ss+1);
HA      = CI.TP.DeltaHr(indexHA)./CI.TP.c_mean(2,ss+1).^2;
% ----------------------------------       
B1a(1,1) =   0;
B1a(1,2) =   cRatio;
B1a(1,3) =   M1.*cRatio;
B1a(2,1) =   1;
B1a(2,2) = 2*M1;
B1a(2,3) =   M1.^2;
B1a(3,1) =   M1.*gamma1./(gamma1-1);
B1a(3,2) =   M1.^2 - HA;
B1a(3,3) =  -M1.*(1/(gamma1-1) + HA);
%
B2(1,1) =   0;
B2(1,2) =   1;
B2(1,3) =   M2;
B2(2,1) =   1;
B2(2,2) = 2*M2;
B2(2,3) =   M2.^2;
B2(3,1) =   cRatio.*M2.*gamma2./(gamma2-1);
B2(3,2) =   cRatio.*M2.^2;
B2(3,3) =  -cRatio.*M2./(gamma2-1); 
%
% -------------------------------------------------------------------------
%
function [B1,B2] = Fcn_Matrix_calculation_WO_Addition_effect(ss,index1,index2)
global CI
cRatio  = CI.TP.c_mean(index2,ss+1)./CI.TP.c_mean(index1,ss);
M1      = CI.TP.M_mean(index1,ss);
M2      = CI.TP.M_mean(index2,ss+1);
gamma1  = CI.TP.gamma(index1,ss);
gamma2  = CI.TP.gamma(index2,ss+1);
Theta   = CI.TP.Theta(ss);
rho1    = CI.TP.rho_mean(index1,ss);
rho2    = CI.TP.rho_mean(index2,ss+1);
% ----------------------------------       
B1(1,1) =   0;
B1(1,2) =   cRatio;
B1(1,3) =   M1.*cRatio;
B1(3,1) =   M1.*gamma1./(gamma1-1);
B1(3,2) =   M1.^2;
B1(3,3) =  -M1./(gamma1-1);
%
B2(1,1) =   0;
B2(1,2) =   Theta;
B2(1,3) =   Theta.*M2;
B2(3,1) =   cRatio.* Theta.*M2.*gamma2./(gamma2-1);
B2(3,2) =   cRatio.* Theta.*M2.^2;
B2(3,3) =  -cRatio.* Theta.*M2./(gamma2-1); 
% ---------------------------------- 
 if Theta>=1
    B1(2,1) =   Theta;
    B1(2,2) = 2*M1;
    B1(2,3) =   M1.^2;
    B2(2,1) =   Theta;
    B2(2,2) = 2*Theta.*M2;
    B2(2,3) =   Theta.*M2.^2;
 elseif Theta<1
    B1(2,1) = 1/rho1.^gamma1;
    B1(2,2) = 0;
    B1(2,3) =-1/rho1.^gamma1;
    B2(2,1) = 1/rho2.^gamma2;
    B2(2,2) = 0;
    B2(2,3) =-1/rho2.^gamma2;
 end
%
%----------------------------------end-------------------------------------