clear; clc; close all;
warning off backtrace
%% Example continuation problem of a 2 DoF "infinite" broaching operation
% an infinitely long tool with a uniform pitch is considered to enable 
% asymptotic stability analysis of periodic orbits in the "plateau" phase 
% of the broaching process, where cutting forces are close to periodic
% state dependent delays are expected due to axial vibrations of the tool
% the tool is infinitely rigid (the cutting edges move in unison)
% a similar model is studied in: https://doi.org/10.1007/s11012-024-01831-0
% THIS FORM OF THE DELAY DEFINITION IS SUSCEPTIBLE TO DRIFTING!

% Dependencies
addpath(genpath('_toolbox'))
addpath(genpath('pwsdde cont'))
addpath(genpath('plot tools'))

% Default plot options
set(0, 'DefaultLineLineWidth', 1);
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultTextInterpreter','latex');
set(0,'DefaultAxesTickLabelInterpreter','latex');
set(0,'DefaultLegendInterpreter','latex');


%% Governing PWS-DDE with state dependent delays
% Equation of motion
% m*y''(t) + cy*y'(t) + ky*y(t) = chi*F_cut     (normal direction)
% m*z''(t) + cz*z'(t) + kz*z(t) = -F_cut        (axial direction)
% Linear cutting force (Kc: cutting coefficient, w: bore widht, h0: feed)
% F_cut(t) = n(t)*Kc*w*(h0 + y(t-tau(t)) - y(t))
% Implicitely defined state dependent delay:
% v*tau(t) + z(t) - z(t-tau) = d
% Number of teeth in contact (L: workpiece length, d: gap between teeth)
% n(t) = floor(L/d) + 1 if v*t + z(t) < mod(L,d)
%      = floor(L/d) if v*t + z(t) > mod(L,d)
% z(t) + v*t has to be reset to zero at every new tooth entry event!
%
% Dimensionless, autonomous, 1st order form (yu=h0, zu=d, tu=d/v) 
% x1'(t) = x2(t)
% x2'(t) = -c1*x2(t) - k1*x1(t) + xi*f(t)
% x3'(t) = x4(t)
% x4'(t) = -c2*x4(t) - k2*x3(t) - f(t)
% x5'(t) = x4(t-x5(t)) - x4(t)
% x6'(t) = 1
% f(t) = w0*(1 + x1(t-x5(t)) - x1(t)) * (nt or nt + 1), 
% c1 = cy*d/(m*v), k1 = ky*d^2/(m*v^2), c2 = cz*d/(m*v), k2 = kz*d^2/(m*v^2)
% w0 = Kc*w*d*h0/(m*v^2), xi = chi*d/h0, nt = floor(L/d)
% Tooth entry event: n_min -> n_max
% h1(t) = x6(t) + x3(t) - 1, g1(x) = [x1, x2, x3, x4, x5, -x3]^T
% Tooth exit event: n_max -> n_min
% h2(t) = x6(t) + x3(t) - delta, g2(x) = x, d0 = mod(L,d)/d


%% System definition
addpath('_system def/infinite broaching sddde');
sys.f = 'inf_broach_f';     % vector field modes
sys.e = 'inf_broach_e';     % event functions
sys.tau = 'inf_broach_tau'; % time delays
sys.tau_no = 1;             % number of time delays
sys.mode_no = 2;            % number of modes
sys.event_no = 2;           % number of possible events
sys.sd_delay = true;        % state depentent delays are present
sys.aut = false; % (just to push the eigenvalue at 1 into the unit circle)

% Dimensionless system paramters (non-physical example dataset)
c1 = 0.5;   % 1) damping in normal direction 
k1 = 10;    % 2) stiffness in normal direction 
c2 = 1;     % 3) damping in axial direction 
k2 = 100;   % 4) stiffness in normal direction
w0 = 0.01;  % 5) cutting constant for the axial direction
xi = 30;    % 6) cutting force ratio
nt = 2;     % 7) minimum number or teeth in contact (n_min)
d0 = 0.7;   % 8) distance between tooth entries and exits
par = [c1 k1 c2 k2 w0 xi nt d0];


%% Finding a periodic orbit

% Periodic orbit initialisation via simulation
y0 = @(t) [0 0 0 0 1 t];        % initial history function
o_init.N = 2;                   % number of events expected in the periodic orbit
o_init.M = 20;                  % Chebyshev mesh dimension
sim_opts.m0 = 1;                % starting mode of the simulation
sim_opts.t_end = 20;            % simulation time
[orb0,res]= sim_ns_sd_dde(y0,par,sys,o_init,sim_opts);
figure(); plot_res(res,1); title('Transient simulation result');
% figure(); plot_orb(orb0,sys,'inf_broach_p'); title('Simulation guess');

% Correct solution guess with Newton iteration
[orb1,err] = orb_corr(orb0,sys);
figure(); plot_orb(orb1,sys,'inf_broach_p'); title('Corrected orbit');

% Visualize the stability of the found orbit
 % [mu1,~,~,v1] = orb_stab(orb1,sys);
 % figure(); plot_spectrum(mu1); % slightly inaccurate trivial eigenvalue at 1
 % figure(); plot(v1); % critical eigenvector points in the direction of x5(t)


%% Follow the periodic orbit in bore width

% Continuation run 1) in w0
opts1 = br12_opts(5); % position of w0 in par
opts1.stop.p_lim = [0 0.1]; % [w0_min, w0_max]
opts1.psa.ds0 = 1e-3; % starting stepsize
opts1.psa.ds_lim = [1e-3 0.1]; % [ds_min ds_max]
branch1 = br12_cont_adapt(orb1,sys,opts1); 
figure(); plot_br1_norm(branch1,opts1.pi);
% anim_br_spectrum(branch1)
% anim_br_orb(branch1,sys,'inf_broach_p')

% test the stability at the critical point
% [~,~,~,~,ibif] = get_br_data(branch1,5,1);
% orbp = branch1(ibif.i_sc+1); % unstable orbit right after the Hopf point
% sim_opts.m0 = 2; sim_opts.t_end = 50; 
% res = sim_stab_test(orbp,sys,sim_opts,1e-3); % perturbed simulation
% figure(); plot_res(res,1); title('Stability test');
% figure(); plot_spectrum(orbp.mu);

%% Follow the periodic orbit in bore length

% Continuation run 2) in d0
opts2 = br12_opts(8); % position of d0 in par
opts2.stop.p_lim = [0 1]; % [d0_min, d0_max]
opts2.psa.ds0 = 1e-3; % starting stepsize
opts2.psa.ds_lim = [1e-3 0.1]; % [ds_min ds_max]
branch2 = br12_cont_adapt(orb1,sys,opts2); 
figure(); plot_br1_ampl(branch2,opts2.pi,1);
% anim_br_spectrum(branch2)
% anim_br_orb(branch2,sys,'inf_broach_p')

