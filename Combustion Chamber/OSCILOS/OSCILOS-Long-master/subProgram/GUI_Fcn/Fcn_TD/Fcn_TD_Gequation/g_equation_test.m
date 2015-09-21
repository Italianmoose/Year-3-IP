clear all
close all

% G_equation test function
% This is a test function for the g-equation SOCILOS implementation

% generate required variables for test 
ra = 0.035;
rb = 0.07;
U1 = 12;
rho1 = 1.12;
area_ratio = 1.0 -(ra/rb)^2;
Ugs = Fcn_TD_Gequ_calc_ugutter( U1,area_ratio,0,0);
SU = 0.088 * Ugs;
nb_points = 60;
y_vec = linspace(ra,rb,nb_points);
deltah = 1E4;
Q_mean = pi * rb^2 * rho1 * Ugs * deltah * area_ratio; % Q_mean needs to be calculated based on the actual area of the combustor, where the flame holder is (hence the multiplication by aratio)
dt = 1E-4;
nb_iter = 1E3;

% Compute the initial steady flame shape
[ xi_steady ] = Fcn_TD_Gequ_steady_flame( Ugs,SU,y_vec );

figure
title('Steady flame shape')
plot(xi_steady,y_vec);
xlabel('xi');
ylabel('r [m]');

% Now define the time iteration, and teh value of uratio
end_time = dt*(nb_iter - 1);
t_vec = linspace(0,end_time,nb_iter);
amp = 0.2;
freq = 333;
alpha = 2;
beta = 5.4;
uratio = amp * (alpha * cos(freq* t_vec) + beta * sin(freq * t_vec));

% and compute the value of xi based on these values of uratio
% First with a loop
figure
xi = xi_steady;
q_ratio = zeros(1,nb_iter);
time_integration = 1;
bash_data = zeros(1,nb_points,3);
for runner = 2:nb_iter % first time step is steady state
    [q_ratio(runner),xi,bash_data] = Fcn_TD_Gequ_interface( SU, xi, y_vec, dt, t_vec(runner), U1, area_ratio, uratio(runner),Q_mean, deltah,rho1,bash_data,runner-1,time_integration );
    
    title('real time flame')
    plot(xi,y_vec);
    xlabel('xi');
    ylabel('r [m]');
    drawnow 
end

% Sinon on essaye avec un calcul vectoriel
% xi =ones(nb_iter,1) * xi_steady; % The xi is initialised to xi_steady for all iterations (values of uratio) as xi is not updated between these values of uratio
% bash_data = zeros(nb_iter,nb_points,3);
% IT = 1;
% time_integration = 1;
% [q_ratio,xi, bash_data] = Fcn_TD_Gequ_interface( SU, xi, y_vec, dt, 0, U1, area_ratio, uratio,Q_mean, deltah,rho1,bash_data,IT,time_integration );

figure
plot(t_vec,q_ratio)
title('q ratio as a function of time')
xlabel('time [s]')
ylabel('q ratio')