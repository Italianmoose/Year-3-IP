function z = Fcn_lsim_varied_td(y,Sys,tSys,Var,td,dt)
% This function is used to simulate time response of dynamic system to
% arbitary inputs with varied time delay.
% y denotes the arbitary input, which is generally a signal including the
% previous "history" of time period (Var(1):Var(2)).*dt.
% Var(1):Var(2) denotes the index of a period of y.
% dt is the time step.
% Sys is the transfer function of the dynamic system.
% tSys is a suggested length of the Green's function of Sys, which is
% suggested to calculate prior of this function to save the calculation
% time.
% td denotes the time delay, which can be an array of different values.
%
% ** In case the dynamic system is the function of s, lsim is a perfect
% function to calculate the response to arbitary inputs.
% ** In case the dynamic system is the function of f, lsim cannot be used
% and it is better to first to calculate the Green's function, then do
% convolution between the arbitary input and the Green's function. In order
% to have a good precesion, interpolation is necessary to process the
% arbitary input and the Green's function should have a higher resolution.
% After convolution, the resolution is then changed to the required one.
%
% author: Jingxuan LI (jingxuan.li@imperial.ac.uk)
% first created:    2014-11-18
% last edited:      2014-11-19
%
nfloor      = floor(td./dt);
nMax        = max(nfloor);                          % maximun discreted samples of td
nMin        = min(nfloor);                          % minimum discreted samples of td
nSys        = ceil(1.2*tSys./dt);                   % increase the length of Green's function to have a better precesion
%
index1      = Var(1) - (nMax + 1) - (nSys - 1);     % start index of padded input
index2      = Var(2) - nMin;                        % end index of padded input

if index1 == index2 % This occurs if we have one operation per time step, as in the G-equation
    index2 = index2 + 1 ; % in this case simulate two time steps. The interpolation will take care of getting the correct value
    index1 = index1 - 1;
end

% if index2 >= min(Var(1):Var(2)) % not correct
%     error('You are trying to use data which has not yet been calculated in an interpolation. Try reducing your time step')
% end

t0          = (0:index2-index1).*dt;                % corresponding time samples

tIntp       = ((Var(1):Var(2)) - index1).*dt - td;  % time samples corresponded to the arbitary input

% lsim
z0          = lsim(Sys,y(index1:index2),t0);        % lsim
%
z           = interp1(t0,z0,tIntp,'linear','extrap'); % do interpolation to have required response




%
% ------------------------ end --------------------------------------------
