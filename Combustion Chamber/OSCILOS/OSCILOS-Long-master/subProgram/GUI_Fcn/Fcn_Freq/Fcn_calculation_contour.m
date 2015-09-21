% ----------------------------Calculate contour values---------------------
%
function F = Fcn_calculation_contour(CalStyle)
global CI
%
n   = length(CI.EIG.Cont.GRSp);
m   = length(CI.EIG.Cont.FreqSp);
F   = zeros(n,m);
for ss = 1:n
    GR = CI.EIG.Cont.GRSp(ss);
    for kk  = 1:m
        omega       = 2*pi*CI.EIG.Cont.FreqSp(kk);
        s           = GR+1i*omega;
        switch CalStyle
            case 1                                                                
                F(ss,kk) = Fcn_DetEqn_Linear(s);                             
            case 2   
                F(ss,kk) = Fcn_DetEqn_frozen_nonlinear(s); 
            case 3
                F(ss,kk) = Fcn_DetEqn_nonlinear_dampers(s);
        end
    end
end
% -------------------------------end--------------------------------------