clear; clc; close all;
warning off backtrace
%% Example continuation problem of 1 DoF dry friction oscillator
% using nondimensionalized form of governing equation of motion
% artificial DDE formalism, phase shift cosidered as time delay

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
% 1) positive sliding (x'>0)
% m*x''(t) + c*x'(t) + k*x(t) = F_0*cos(Omega*t-phi) - mu*m*g
% 2) negative sliding (x'<0)
% m*x''(t) + c*x'(t) + k*x(t) = F_0*cos(Omega*t-phi) + mu*m*g
% 3) sticking (x'=0)
% x(t) = 0
% Dimensionless (xu=mu*m*g/k, tu=1/om_n), autonomous, 1st order form
% 1) positive sliding (x'>0)
% x1'(t) = x2(t)
% x2'(t) = -2*zeta*x2(t) - x1(t) + f_0*cos(x3(t-tau)) - 1
% x3'(t) = om
% 2) negative sliding (x'<0)
% x1'(t) = x2(t)
% x2'(t) = -2*zeta*x2(t) - x1(t) + f_0*cos(x3(t-tau)) + 1
% x3'(t) = om
% 3) sticking (x'=0)
% x1'(t) = x2(t)
% x2'(t) = 0
% x3'(t) = om
% Contact surfaces and event maps:
% - direction reversal, transition to sticking
%       h1(x(t),x(t-tau) = x2(t) 
% - separation from sticking
%       h2(x(t),x(t-tau)) = -x1(t)-2*zeta*x2(t)+f_0*cos(x3(t-tau)) +- eta
% - phase reset 
%       h3(x(t),x(t-tau) = x3(t) 
%       g3(x(t),x(t-tau) = [x1(t); x2(t); 0];
% System paramters (om_n=sqrt(k1/m):
% 1) zeta = c/(2*m*om_n)    damping ratio
% 2) om = Omega/om_n        dimensionless excitation frequency
% 3) tau = phi/om           dimensionless time delay (phase shift)
% 4) f_0 = F_0/(mu*m*g)     dimensionless excitation amplitude
% 5) eta = mu_0/mu          friction coefficient ratio


%% System definition
addpath('#demo');
sys.f = 'f_fricosc';        % vector field modes
sys.e = 'e_fricosc';        % event functions
sys.tau = 'tau_fricosc';    % time delays
sys.mode_no = 3;            % number of modes
sys.event_no = 9;           % number of possible events


%% Find a simple periodic orbit

% Periodic orbit initialisation via simulation
p0 = [0.05, 0.9, 1.0, 1.5, 1.1];   % parameter vector [zeta, om, tau, f_0, eta]
y0 = @(t) [0 5.0 p0(2)*t];  % initial history function
o_init.N = 3;   % number of events
o_init.M = 20;  % Chebyshev mesh dimension
sim_opts.h_act = [1 1 0 0 0 0 1 1 0];   % active events during simulation (no sliding motion)
sim_opts.m0 = 1;    % initial vector field mode
sim_opts.t_end = 50*2*pi/p0(2); % simulation time
[orb1,res]= sim_ns_dde(y0,p0,sys,o_init,sim_opts);
figure(); plot_res(res,2); title('Transient simulation result');
figure(); plot_orb(orb1); title('Simulation guess');

% Correct solution guess with Newton iteration
corr_opts.nr.logs = true; % log iteration progress
% corr_opts.nr.plots = true;  % plot iteration progress
[orb1c,err] = orb_corr(orb1,sys,corr_opts);
figure(); plot_orb(orb1c,'p_fricosc'); title('Corrected orbit');
% figure(); plot(err);


%% Follow simple periodic orbits in excitation frequency

% Continuation options
cont_opts.psa.ds = -0.2;    % arclenght step (negative direction)
cont_opts.psa.pi = 2;       % continuation parameter index (om)
cont_opts.psa.np = 100;     % number of contination steps
cont_opts.psa.igr_stop = false; % dont stop at false interior grazing points

% Continuation run 1)
branch1 = br12_cont(orb1c,sys,cont_opts);
% figure(); plot_br1_ampl(branch1,cont_opts.psa.pi,1,false);
% figure(); plot_br1_norm(branch1,cont_opts.psa.pi,2);
% anim_br_orb(branch1,'p_fricosc','test.avi')

% Continuation run 2)
cont_opts.psa.ds = 0.5; % (positive direction)
branch2 = br12_cont(orb1c,sys,cont_opts);
% figure(); plot_br1_ampl(branch2,cont_opts.psa.pi,1,false);
% figure(); plot_br1_norm(branch1,cont_opts.psa.pi,2);
% figure(); plot_orb(branch2(end));


%% Follow orbits with sticking segments

% project sliding orbit to a new mesh to follow sticking orbits
t_guess = [0.0, 0.1996, 2.2077, 4.3246, 4.5423, 8.7036]; % guess for event locations 
orb2 = orb_convert(branch1(end),sys,t_guess); % orbit conversion
% orb2 = orb_convert(branch1(end),sys,6); % with manual event selection
orb2.sig = [4 7 3 6 5]; % new solution signature
% figure(); plot_orb(orb2)

% correct sticking orbit
orb2c = orb_corr(orb2,sys,corr_opts);
% figure(); plot_orb(orb2c)

% Continuation run 1)
cont_opts.psa.ds = -0.5;
branch_s1 = br12_cont(orb2c,sys,cont_opts);
% figure(); plot_br1_ampl(branch_s1,cont_opts.psa.pi,1,false);
% figure(); plot_orb(branch_s1(end),'p_fricosc')

% change solution signature to allow going further
orb3 = branch_s1(end);
orb3.T(orb3.T<0)=-orb3.T(orb3.T<0);
orb3.sig = [4 3 9 6 5];
orb3c = orb_corr(orb3,sys,corr_opts);
% figure(); plot_orb(orb3c)

% Continuation run 2)
cont_opts.psa.ds = -1;
branch_s2 = br12_cont(orb3c,sys,cont_opts);
% figure(); plot_br1_ampl(branch_s2,cont_opts.psa.pi,1,false);
% figure(); plot_orb(branch_s2(end),'p_fricosc')


%% Follow a sliding bifurcations

% find and correct initial sliding point
orb_sl = branch1(end);
bifs.type = 2;  % sliding bifurcation
bifs.ind = 2;   % at the second event (also third due to simmetry, but one is enough)
bifs.pi = 2;    % free paramter needed to locate the bifurcation point
orb_slc = orb_corr(orb_sl,sys,corr_opts,bifs);
figure(); plot_orb(orb_slc,'p_fricosc'); title('Sliding orbit')

% Follow sliding orbit in 2 parameters
cont_opts.psa.sl_stop = false; % don't stop at trivial sliding events (due to simmetry)
cont_opts.psa.ds = 0.5;
cont_opts.psa.pi = [4,2]; % two continuation parameters 
% (switch parameter order for continuation in the opposite direction)
branch_sl1 = br12_cont(orb_slc,sys,cont_opts,bifs);
% figure(); plot_br2_par(branch_sl1,cont_opts.psa.pi,2)
% figure(); plot_orb(branch_sl1(end),'p_fricosc')


%% Plot continuation results

% combined plot in 3D
pind = [2,4];
figure();
plot_br2_3D(branch1,pind,1,2); hold on
plot_br2_3D(branch2,pind,1,2);
plot_br2_3D(branch_s1,pind,1,2);
plot_br2_3D(branch_s2,pind,1,2);
plot_br2_3D(branch_sl1,pind,1,2); hold off
xlabel('$\omega$'); ylabel('$f_0$'); title('Continuation results')
