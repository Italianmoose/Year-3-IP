function y = Fcn_calculation_envelope(y0,nPadding)
% This function is used to calculate the envelope of a singal
% y0 is the original signal
% nPadding is the length of zero padding signals before and after the
% original signal y0
%
% author: Jingxuan LI (jingxuan.li@imperial.ac.uk)
% first created:    2014-11-18
% last edited:      2014-11-19
%
IsRowArray  = 1;
SizeY       = size(y0);
if SizeY(1)>SizeY(2)
    y0 = y0';               % The array should be 1xn style
    IsRowArray = 0;
end
%
yP      = [zeros(1,nPadding), y0, zeros(1,nPadding)];                                       
yEnv    = abs(hilbert(detrend(yP)));
y       = yEnv(nPadding+1:end-nPadding);
%
if IsRowArray == 0
    y = y';
end
%
% ----------------------------end------------------------------------------



