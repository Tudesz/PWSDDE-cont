clear; clc; close all;
warning off backtrace
%% Example continuation problem of a 1 DoF dry friction oscillator
% using a nondimensionalized form of the governing equation of motion
% artificial DDE formalism, phase shift cosidered as time delay

% Dependencies
addpath(genpath('_toolbox'))
addpath('pwsdde cont','plot tools')

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
% x(t) = const
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
addpath('_system def/dry friction oscillator');
sys.f = 'fricosc_f';        % vector field modes
sys.e = 'fricosc_e';        % event functions
sys.tau = 'fricosc_tau';    % time delays
sys.mode_no = 3;            % number of modes
sys.event_no = 9;           % number of possible events
sys.tau_no = 1;             % number of time delays


%% Find a simple periodic orbit

% Periodic orbit initialisation via simulation
p0 = [0.05 0.9 1.0 1.5 1.1];    % parameter vector [zeta, om, tau, f_0, eta]
y0 = @(t) [0 5.0 p0(2)*t];      % initial history function
o_init.N = 3;                   % number of events expected in the periodic orbit
o_init.M = 20;                  % Chebyshev mesh dimension
sim_opts.m0 = 1;                % initial vector field mode
sim_opts.t_end = 50*2*pi/p0(2); % simulation time
sim_opts.h_act = [1 1 0 0 0 0 1 1 0]; % active events during simulation (no sliding motion)
[orb1,res]= sim_ns_dde(y0,p0,sys,o_init,sim_opts);
figure(); plot_res(res,2); title('Transient simulation result');
% figure(); plot_orb(orb1,sys); title('Simulation guess');

% Correct solution guess with Newton iteration
[orb1c,err] = orb_corr(orb1,sys);
figure(); plot_orb(orb1c,sys,'fricosc_p'); title('Corrected orbit');
% figure(); plot(err);


%% Follow simple periodic orbits in excitation frequency

% Continuation run 1) in om
opts1.pi = 2; % position of om in p0
opts1.stop.int_gr = false; % do not stop at false grazing events
branch1 = br12_cont_adapt(orb1,sys,opts1);
% figure(); plot_br1_ampl(branch1,opts1.pi,1,false);
% figure(); plot_br1_norm(branch1,opts1.pi,2);
% anim_br_orb(branch1,sys,'fricosc_p','test.avi')
% figure(); plot_orb_events(branch1(end),sys)

% Using a fixed stepsize (run 1)
% cont_opts = opts1;  % extend the previous options structure
% cont_opts.np = 100; % number of contination steps
% cont_opts.ds = -0.2;% arclenght step (negative direction)
% branch1 = br12_cont_fix(orb1c,sys,cont_opts);
% cont_opts.ds = 0.5; % flip the continuation direction
% branch2 = br12_cont_fix(orb1c,sys,cont_opts);
% branch1 = cat(1,flip(branch1),branch2); % concatenate the two branches


%% Follow orbits with sticking segments

% project sliding orbit to a new mesh to follow sticking orbits
t_guess = [0.0 0.163 2.196 4.349 4.446 8.667]; % guess for adaptive branch1
% t_guess = [0.0 0.2 2.207 4.324 4.542 8.704]; % guess for fixed step branch1
orb2 = orb_convert(branch1(1),sys,t_guess); % orbit conversion
% orb2 = orb_convert(orb20,sys,6); % with manual event selection
orb2.sig = [4 7 3 6 5]; % new solution signature (with sticking events)
% figure(); plot_orb(orb2,sys)

% correct sticking orbit
orb2c = orb_corr(orb2,sys);
% figure(); plot_orb(orb2c,sys)

% Continuation run (2) in om
opts_s1.pi = 2; % position of om in p0
opts_s1.stop.int_gr = false; % do not stop at false grazing events
opts_s1.stop.n_step = [100 0]; % only do continuation in -om
branch_s1 = br12_cont_adapt(orb2,sys,opts_s1);
% figure(); plot_br1_ampl(branch_s1,opts_s1.pi,1,false);
% figure(); plot_orb(branch_s1(end),sys,'fricosc_p')

% Using a fixed stepsize (run 2)
% cont_opts.ds = -0.5;
% branch_s1 = br12_cont_fix(orb2c,sys,cont_opts);

% Change solution signature to allow going further
orb3 = branch_s1(1); % initial orbit from adaptive branch_s1
% orb3 = branch_s1(end); % initial orbit from fixed step branch_s1
orb3.T(orb3.T<0)=-orb3.T(orb3.T<0); % remove negative segment lengths
orb3.sig = [4 3 9 6 5]; % update the solution signature
orb3c = orb_corr(orb3,sys);
% figure(); plot_orb(orb3c,sys)

% Continuation run (3) in om
opts_s2.pi = 2; % position of om in p0
opts_s2.stop.int_gr = false; % do not stop at false grazing events
opts_s2.stop.n_step = [100 0]; % only do continuation in -om
branch_s2 = br12_cont_adapt(orb3c,sys,opts_s2);
% figure(); plot_br1_ampl(branch_s2,opts_s2.pi,1,false);
% figure(); plot_orb(branch_s2(end),sys,'fricosc_p')

% Using a fixed stepsize (run 3)
% cont_opts.ds = -1;
% branch_s2 = br12_cont_fix(orb3c,sys,cont_opts);


%% Follow a branch of sliding bifurcations

% find and correct initial sliding point
orb_sl = branch1(1); % initial orbit from adaptive branch1
% orb_sl = branch1(end); % initial orbit from fixed step branch1
bifs.type = 2;  % sliding bifurcation
bifs.ind = 2;   % at the second event (also third due to simmetry, but one is enough)
bifs.pi = 2;    % free paramter needed to locate the bifurcation point
orb_slc = orb_corr(orb_sl,sys,bifs);
% figure(); plot_orb(orb_slc,sys,'fricosc_p'); title('Sliding orbit')

% Continuation run (4) in f_0 and om
opts_sl1.pi = [4 2]; % position of f_0 and om in p0
opts_sl1.stop.int_gr = false; % do not stop at false grazing events
opts_sl1.stop.slide = false; % do not stop at false sliding events
opts_sl1.stop.p_lim = [-inf 1.5; -inf inf]; % [f_0_min f_0_max; om_min om_max]
opts_sl1.psa.ds_lim = [1e-3 0.5]; % [ds_min ds_max]
branch_sl1 = br12_cont_adapt(orb_slc,sys,opts_sl1,bifs);
% figure(); plot_br1_ampl(branch_sl1,opts_sl1.pi,1,false);
% figure(); plot_orb(branch_sl1(end),sys,'fricosc_p')

% Using a fixed stepsize (run 4)
% cont_opts = opts_sl1;   % extend the previous options structure
% cont_opts.np = 100;     % number of contination steps
% cont_opts.ds = 0.5;     % arclenght step
% branch_sl1 = br12_cont_fix(orb_slc,sys,cont_opts,bifs);


%% Plot continuation results

% combined bifurcation diagram in 3D
pind = [2,4];
figure();
plot_br2_3D(branch1,pind,1,2); hold on
plot_br2_3D(branch_s1,pind,1,2);
plot_br2_3D(branch_s2,pind,1,2);
plot_br2_3D(branch_sl1,pind,1,2); hold off
xlabel('$\omega$'); ylabel('$f_0$'); title('Continuation results')


%% Recreate the continuation results using COCO
% 
% % Initialize COCO
% addpath('<COCO_dir>'); % location of COCO installation
% startup;
% 
% % define a COCO compatible problem
% p_names = {'$\zeta$', '$\omega$', '$\tau$', '$f_0$', '$\eta$'};
% prob1 = pwsdde_coco_prob(sys,orb1,p_names);
% 
% % Define events for vanishing segments and sliding bifurcations
% prob1 = pwsdde_coco_ev_vanish(prob1,1e-3);
% prob1 = pwsdde_coco_ev_slide(prob1,sys);
% 
% % COCO numeric options
% prob1 = coco_set(prob1,'cont','h_min',1e-3,'h_max',0.5);
% prob1 = coco_set(prob1,'cont','ItMX', 100);
% 
% % create a branch of solutions with COCO (in om)
% bd1 = coco(prob1,'fro_r1',[],1,{'$\omega$','sl_2'},[0.5, 1.5]);
% 
% % convert to a sliding solution
% lab = coco_bd_labs('fro_r1','SL');
% orb_sl = pwsdde_coco_2orb('fro_r1',prob1,lab(1));
% % plot_orb(orb_sl,sys,'fricosc_p')
% 
% % define a new COCO compatible problem
% prob2 = pwsdde_coco_prob(sys,orb_sl,p_names);
% 
% % Include an extra zero condition to ensure sliding
% prob2 = pwsdde_coco_h_slide(prob2,sys,2);
% 
% % create a branch of sliding solutions with COCO (in om and f0)
% bd_sl = coco(prob2,'fro_sl',[],1,{'$\omega$','$f_0$'},[0.5, 1.5]);
% 
% % Plot continuation results with COCO
% thm = struct();
% thm.special = {'EP','VA','SL'};
% thm.VA = {'dk','MarkerSize',7,'MarkerFaceColor','y'};
% thm.SL = {'ok','MarkerSize',7,'MarkerFaceColor','g'};
% figure(); box on; hold on
% coco_plot_bd(thm,'fro_r1','$\omega$','$f_0$','Ampl')
% coco_plot_bd(thm,'fro_sl','$\omega$','$f_0$','Ampl')
% hold off
% 
% % Plot an individual orbit
% % lab = coco_bd_labs('fro_r1','SL');
% % orbp = pwsdde_coco_2orb('fro_r1',prob1,lab(1));
% % figure();
% % plot_orb(orbp,sys,'fricosc_p');
