function yf = Fcn_filtfilt_butterworth(y,fSample,fPass,fStop,Rp,Rs,filterType)
% This function is a Zero-phase Butterworth filter.  
% fSample represents the sampling rate
% fPass represents the cut-off frequency
% fStop represents the stop frequency
% Rp denotes Passband ripple in decibels
% Rs: denotes Stopband attenuation in decibels. This value is the number of decibels
% the stopband is down from the passband. 
% filterType is the type of filter
% FilterGroup = {'low','bandpass','high','stop' };
%
% author: Jingxuan LI (jingxuan.li@imperial.ac.uk)
% first created:    2014-11-18
% last edited:      2014-11-19
%
FilterGroup = {'low','bandpass','high','stop' };
%
switch filterType
    case 1
        Wp = fPass(1)/(fSample/2);
        Ws = fStop(1)/(fSample/2);
    case 2
        Wp = fPass(1)/(fSample/2);
        Ws = fStop(1)/(fSample/2);
    case 3
        Wp = fPass(1:2)./(fSample/2);
        Ws = fStop(1:2)./(fSample/2);
    case 4
        Wp = fPass(1:2)./(fSample/2);
        Ws = fStop(1:2)./(fSample/2);
end
[n,Wn]  = buttord(Wp,Ws,Rp,Rs);
[b,a]   = butter(n,Wn,FilterGroup{filterType}); 
yf      = filtfilt(b,a,y); 
%
% -----------------------------end-----------------------------------------
