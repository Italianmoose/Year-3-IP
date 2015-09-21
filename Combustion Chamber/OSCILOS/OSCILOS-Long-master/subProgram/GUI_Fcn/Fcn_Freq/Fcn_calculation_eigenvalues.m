% ----------------------------Calculate eigenvalue-------------------------
%
function E = Fcn_calculation_eigenvalues(CalStyle)
% This function is used to calculate the eigenvalue 
global CI
%---------------------------------
% set the scan domain
FreqMin     = CI.EIG.Scan.FreqMin;   
FreqMax     = CI.EIG.Scan.FreqMax;
FreqNum     = CI.EIG.Scan.FreqNum;
GRMin       = CI.EIG.Scan.GRMin;   
GRMax       = CI.EIG.Scan.GRMax;
GRNum       = CI.EIG.Scan.GRNum;
FreqSp      = linspace(FreqMin, FreqMax,    FreqNum);       % sample number of scan frequency
GRSp        = linspace(GRMin,   GRMax,      GRNum);         % sample number of growth rate
%---------------------------------
options = optimset('Display','off');        % the calculation results by fsolve will not be shown in the workspace
eigen_num=1;
for ss = 1:length(GRSp)
    GR = GRSp(ss);
    for kk = 1:length(FreqSp)
        omega = 2*pi*FreqSp(kk);
        if GR == 0
            GR = 1;
        end
        if omega == 0
            omega =10;
        end
        s0 = GR+1i*omega;                                                           % initial value
        switch CalStyle
            case 1                                                                  % linear flame model
                [x,fval,exitflag] = fsolve(@Fcn_DetEqn_Linear,s0,options);                             
            case 2                                                                  % nonlinear flame model
                [x,fval,exitflag] = fsolve(@Fcn_DetEqn_frozen_nonlinear,s0,options); 
            case 3                                                                  % inlcude two situations: there are nonlinear dampers; there are more than one NL explicit FM (and no EXP FM)
                [x,fval,exitflag] = fsolve(@Fcn_DetEqn_nonlinear_dampers,s0,options); 
        end
        if exitflag==1;
            EIG.eigenvalue(eigen_num) = x;
            EIG.eigenvalue_prec(eigen_num) = floor(EIG.eigenvalue(eigen_num)./10).*10;      % this is used to set the precision
            eigen_num = eigen_num+1;
        end
    end
end
[b,m,n] = unique(EIG.eigenvalue_prec);                      % unique function is used to get rid of the same value
EIG.eigenvalue_unique = EIG.eigenvalue(m);                  % this is the eigenvalue
EIG.eigenvalue_unique_prec = EIG.eigenvalue_prec(m);
%---------------------------------
% this is used to get rid of the zero frequency value
s_null = [];
cal = 0;
for ss = 1:length(EIG.eigenvalue_unique)
    if(abs(imag(EIG.eigenvalue_unique(ss)))<1)
        cal = cal+1;
        s_null(cal) = ss;
    end
end
EIG.eigenvalue_unique(s_null) = [];
%---------------------------------
% the previous processing is still not enough to get rid of the same
% value,such as 999.9 and 1000.1 corresponding to the function `floor'
% such as 995.1 and 994.9 corresponding to the function `round'
cal = 0;
EIG.eigenvalue_unique = sort(EIG.eigenvalue_unique);
EIG.eigenvalue_unique_diff = diff(EIG.eigenvalue_unique);
EIG.index_NULL = [];
for kk = 1:length(EIG.eigenvalue_unique_diff)
    if(abs(EIG.eigenvalue_unique_diff(kk))<10)
        cal = cal+1;
        EIG.index_NULL(cal) = kk;
    end
end
EIG.eigenvalue_unique(EIG.index_NULL) = []; % this is the eigenvalue we want

EIG.growthrate_limit    = [GRMin    GRMax];
EIG.Omega_limit         = [FreqMin  FreqMax].*2.*pi;
%-------------------------------------
s_null=[];
cal=0;
for ss=1:length(EIG.eigenvalue_unique)
    if(real(EIG.eigenvalue_unique(ss))<EIG.growthrate_limit(1)||real(EIG.eigenvalue_unique(ss))>EIG.growthrate_limit(2)||abs(imag(EIG.eigenvalue_unique(ss)))<EIG.Omega_limit(1)||abs(imag(EIG.eigenvalue_unique(ss)))>EIG.Omega_limit(2) )
        cal=cal+1;
        s_null(cal)=ss;
    end
end
EIG.eigenvalue_unique(s_null)=[];        % Do not use semicolon since I want to show the final value
E = EIG.eigenvalue_unique; 
clear EIG