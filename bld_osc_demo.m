clear; clc; close all;
warning off backtrace
%% Example continuation problem of 1 DoF bilinear delayed oscillator
% considering harmonic excitation and delayed velocity feedback control
% using nondimensionalized form of governing equation of motion
% further details on the model at: https://doi.org/10.1016/j.cnsns.2019.105095

% Dependencies
addpath('cheb colloc')
addpath('cheb bvp')
addpath('cont routines')
addpath('plot tools')

% Default plot options
set(0, 'DefaultLineLineWidth', 1);
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultTextInterpreter','latex');
set(0,'DefaultAxesTickLabelInterpreter','latex');
set(0,'DefaultLegendInterpreter','latex');


%% Governing PWS-DDE
% Equation of motion
% 1) free mode (x<g)
% m*x''(t) + c*x'(t) + k1*x(t) = F*cos(Omega*t) + K*(x'(t-tau)-x'(t))
% 2) contact mode (x>=g)
%  m*x''(t) + c*x'(t) + k1*x(t) + k2*(x(t)-g) = F*cos(Omega*t) + K*(x'(t-tau)-x'(t))
% Dimensionless (xu=g,tu=1/om_n), autonomous, 1st order form
% 1) free mode (x<1)
% x1'(t) = x2(t)
% x2'(t) = -2*zeta*x2(t) - x1(t) + f*x3(t) + k*(x2(t-tau)-x2(t))
% x3'(t) = -om*x4(t) + (1-x3(t)^2-x4(t)^2)*x3(t)
% x4'(t) = om*x3(t) + (1-x3(t)^2-x4(t)^2)*x4(t)
% 2) contact mode (x>=1)
% x1'(t) = x2(t)
% x2'(t) = -2*zeta*x2(t) - x1(t) -beta*(x1(t)-1) + f*x3(t) + k*(x2(t-tau)-x2(t))
% x3'(t) = -om*x4(t) + (1-x3(t)^2-x4(t)^2)*x3(t)
% x4'(t) = om*x3(t) + (1-x3(t)^2-x4(t)^2)*x4(t)
% Contact surface and event map:
% h(x(t),x(t-tau)) = x1(t)-1
% g(x(t),x(t-tau)) = x(t) 
% System paramters (om_n=sqrt(k1/m):
% 1) tau = om_n*Tau         dimensionless time delay
% 2) k = K/(m*om_n)         dimensionless control parameter
% 3) zeta = c/(2*m*om_n)    damping ratio
% 4) beta = k2/k1           stiffness ratio
% 5) f = F/k1               dimensionless excitation amplitude
% 6) om = Omega/om_n        dimensionless excitation frequency


%% System definition
addpath('#demo');
sys.f = 'f_bldosc';     % vector field modes
sys.e = 'e_bldosc';     % event functions
sys.tau = 'tau_bldosc'; % time delays
sys.mode_no = 2;        % number of modes
sys.event_no = 3;       % number of possible events


%% Finding a periodic orbit

% Periodic orbit initialisation via simulation
p0 = [3.8 0.12 0.01 29 0.9*0.8^2/1.26 0.8]; % parameter vector [tau, k, zeta, beta, f, om]
y0 = @(t) [0 0 cos(p0(6)*t) sin(p0(6)*t)];  % initial history function
o_init.N = 2;                   % number of events
o_init.M = 20;                  % Chebyshev mesh dimension
sim_opts.h_act = [1 1 0];       % active events during simulation
sim_opts.m0 = 1;                % starting mode of the simulation
sim_opts.t_end = 50*2*pi/p0(6); % simulation time
[orb0,res]= sim_ns_dde(y0,p0,sys,o_init,sim_opts);
figure(); plot_res(res,1); title('Transient simulation result');
% figure(); plot_orb(orb0,'p_bldosc'); title('Simulation guess');

% Correct solution guess with Newton iteration
corr_opts.nr.logs = true; % log iteration progress
% corr_opts.nr.plots = true;  % plot iteration progress
[orb1,err] = orb_corr(orb0,sys,corr_opts);
figure(); plot_orb(orb1,'p_bldosc'); title('Corrected orbit');
% figure(); plot(err);


%% Follow a periodic orbit in a system parameter

% Continuation options
cont_opts.psa.ds = 0.1;  % arclenght step
cont_opts.psa.pi = 1;    % continuation parameter index
cont_opts.psa.np = 100;  % number of contination steps

% Continuation run 1)
branch1 = br12_cont(orb1,sys,cont_opts);
% figure(); plot_br1_norm(branch1,1,2);
% figure(); plot_br1_ampl(branch1,1,1);
% anim_br_orb(branch1,'p_bldosc','test.avi')

% Continuation run 2)
cont_opts.psa.pi = 2; % now in k
branch2 = br12_cont(orb1,sys,cont_opts);
% figure(); plot_br1_norm(branch2,2,2);
% figure(); plot_orb(branch2(end),'p_bldosc');


%% Follow a grazing event in 2 parameters

% Initial Grazing point (converted to new solution signature)
orb_gr = orb_convert(branch1(end),sys,0); % only 1 event at the end/beginning
figure(); plot_orb(orb_gr,'p_bldosc'); title('Grazing orbit');

% Correct grazing orbit
orb_gr.sig = 3; % new solution signature with no mode transition
bifs.type = 1;  % grazing bifurcation
bifs.ind = 1;   % at the first (and only) event
bifs.pi = 1;    % free paramter needed to locate the bifurcation point
orb2 = orb_corr(orb_gr,sys,corr_opts,bifs);
% figure(); plot_orb(orb2,'p_bldosc')

% Follow grazing orbit in 2 parameters
cont_opts.psa.pi = [1 2]; % two continuation parameters
branch_gr1 = br12_cont(orb2,sys,cont_opts,bifs);
% figure(); plot_br2_par(branch_gr1,[1 2],2)


%% Plot results
pind = [1,2];

% combined bifurcation diagram in 2D
figure();
plot_br2_par(branch1,pind); hold on
plot_br2_par(branch2,pind);
plot_br2_par(branch_gr1,pind); hold off
xlabel('$\tau$'); ylabel('$k$'); title('Continuation results')

% combined bifurcation diagram in 3D
figure();
plot_br2_3D(branch1,pind,1,2); hold on
plot_br2_3D(branch2,pind,1,2);
plot_br2_3D(branch_gr1,pind,1,2); hold off
xlabel('$\tau$'); ylabel('$k$'); title('Continuation results')


%% Check orbit stability

% normal orbit with crossing
i0 = 10; br = branch1; 
sim_opts.m0 = 1; sim_opts.h_act = [1 1 0];
figure(); 
subplot(1,2,1); plot_spectrum(br(i0).mu); % Floquet multipliers
res = sim_stab_test(br,sys,i0,sim_opts,1e-3); % Simulation started from periodic orbit
subplot(1,2,2); plot_res(res,'p_bldosc',br(i0).p); hold on
plot_orb(br(i0),'p_bldosc'); hold off

% grazing orbit
i0 = 50; br = branch_gr1; 
sim_opts.m0 = 1; sim_opts.h_act = [0 0 0];
figure(); 
subplot(1,2,1); plot_spectrum(br(i0).mu);
res = sim_stab_test(br,sys,i0,sim_opts,1e-1);
subplot(1,2,2); plot_res(res,'p_bldosc',br(i0).p); hold on
plot_orb(br(i0),'p_bldosc'); hold off
