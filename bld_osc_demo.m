clear; clc; close all;
warning off backtrace
%% Example continuation problem of 1 a DoF bilinear delayed oscillator
% considering harmonic excitation and delayed velocity feedback control
% using a nondimensionalized form of the governing equation of motion
% further details on the model at: https://doi.org/10.1016/j.cnsns.2019.105095

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


%% Governing PWS-DDE
% Equation of motion
% 1) free mode (x<g)
% m*x''(t) + c*x'(t) + k1*x(t) = F*cos(Omega*t) + K*(x'(t-tau)-x'(t))
% 2) contact mode (x>=g)
%  m*x''(t) + c*x'(t) + k1*x(t) + k2*(x(t)-g) = F*cos(Omega*t) + K*(x'(t-tau)-x'(t))
% Dimensionless, autonomous, 1st order form (xu=g,tu=1/om_n)
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
addpath('_system def/delayed bilinear oscillator');
sys.f = 'bldosc_f';     % vector field modes
sys.e = 'bldosc_e';     % event functions
sys.tau = 'bldosc_tau'; % time delays
sys.tau_no = 1;         % number of time delays
sys.mode_no = 2;        % number of modes
sys.event_no = 3;       % number of considered events

% add a monitor function for tracking the mean kinetic energy
% sys.q = 'bldosc_mon_Ek_mean'; % q = int(1/2*x2(t)^2)/T


%% Finding a periodic orbit

% Periodic orbit initialisation via simulation
p0 = [3.8 0.12 0.01 29 0.9*0.8^2/1.26 0.8]; % parameter vector [tau, k, zeta, beta, f, om]
y0 = @(t) [0 0 cos(p0(6)*t) sin(p0(6)*t)];  % initial history function
o_init.N = 2;                   % number of events expected in the periodic orbit
o_init.M = 20;                  % Chebyshev mesh dimension
sim_opts.h_act = [1 1 0];       % active events during simulation (grazing event turned off)
sim_opts.m0 = 1;                % starting mode of the simulation
sim_opts.t_end = 50*2*pi/p0(6); % simulation time
[orb0,res]= sim_ns_dde(y0,p0,sys,o_init,sim_opts);
figure(); plot_res(res,1); title('Transient simulation result');
% figure(); plot_orb(orb0,sys,'bldosc_p'); title('Simulation guess');

% Correct solution guess with Newton iteration
[orb1,err] = orb_corr(orb0,sys);
figure(); plot_orb(orb1,sys,'bldosc_p'); title('Corrected orbit');
% figure(); plot(err);


%% Follow a periodic orbit in a system parameter

% Continuation run 1) in tau
opts1 = br12_opts(1); % position of tau in p0
opts1.stop.p_lim = [3.5 5]; % [tau_min, tau_max]
branch1 = br12_cont_adapt(orb1,sys,opts1); 
% figure(); plot_br1_norm(branch1,opts1.pi,2);
% figure(); plot_br1_ampl(branch1,opts1.pi,1);
% anim_br_orb(branch1,sys,'bldosc_p',[],'test.avi')
% figure(); plot_orb_events(branch1(end),sys)

% Continuation run 2) in k
opts2 = br12_opts(2); % position of k in p0
opts2.stop.p_lim = [0.1 0.2]; % [k_min, k_max]
branch2 = br12_cont_adapt(orb1,sys,opts2); 
% figure(); plot_br1_norm(branch2,opts2.pi,2);
% figure(); plot_orb(branch2(end),sys,'bldosc_p');
% figure(); plot_orb_events(branch2(end),sys)

% Using a fixed stepsize
% opts1 = br12_opts(1,0.1,100); % index of tau, ds, n_step
% branch1 = br12_cont_fix(orb1,sys,opts1);
% opts2 = br12_opts(2,0.1,100); % index of k, ds, n_step
% branch2 = br12_cont_fix(orb1,sys,opts2);


%% Follow a grazing event in 2 parameters

% Initial Grazing point (converted to new solution signature) 
orb_gr = orb_convert(branch1(end),sys,0,[],20); % only 1 event at the end/beginning
% figure(); plot_orb(orb_gr,sys,'bldosc_p'); title('Grazing orbit');

% Correct grazing orbit
orb_gr.sig = 3; % new solution signature with no mode transition
bifs.type = 1;  % grazing bifurcation
bifs.ind = 1;   % at the first (and only) event
bifs.pi = 1;    % free paramter needed to locate the bifurcation point
orb2 = orb_corr(orb_gr,sys,[],bifs);
% figure(); plot_orb(orb2,sys,'bldosc_p')

% Adaptive stepsize continuation run
optsgr = br12_opts([1,2]); % position of tau and k in p0
optsgr.stop.p_lim = [1 5; -inf inf]; % [tau_min tau_max; k_min k_max];
optsgr.psa.ds_lim = [1e-3 0.1]; % stepsize limits [ds_min ds_max]
optsgr.psa.ds0 = 0.05; % starting stepsize
branch_gr1 = br12_cont_adapt(orb2,sys,optsgr,bifs);
% figure(); plot_br2_par(branch_gr1,optsgr.pi,2)

% Fixed stepsize continuation run
% optsgr = br12_opts([1,2],0.1,100); % position of tau and k in p0, ds, n_step
% optsgr.stop.p_lim = [1 5; -inf inf]; % [tau_min tau_max; k_min k_max];
% branch_gr1 = br12_cont_fix(orb2,sys,optsgr,bifs);


%% Follow a grazing orbits considering double events
% 
% % Initial Grazing point (converted to new solution signature)
% sys2 = sys; sys2.event_no = 4; % append a dummy event for tracking smooth orbits
% gr_cond = @(x,xd,p) 1 - max(x(1,:)); % add a monitor function to validate grazing
% sys2.q = @(y,orb,sys,pind) po_f_ev(y,orb,sys,pind,gr_cond); % append monitor function
% t_guess = [0.0 3.940 7.868]; % guess for event times [0.0 3.940 7.868]
% orb_gr = orb_convert(branch1(end),sys2,t_guess); % two velocity zero crossings
% % figure(); plot_orb(orb_gr,sys2,'bldosc_p'); title('Grazing orbit');
% 
% % Correct grazing orbit (grazing defined as a double event)
% orb_gr.sig = [4 4]; % new solution signature with no mode transition
% bifs.type = 3;  % user defined bifurcation
% bifs.ind = 2;   % at the second (grazing) event
% bifs.pi = 1;    % free paramter needed to locate the bifurcation point
% bifs.f = @(U,T,p,orb,sys,del,bif_ind,type) ...
%     add_event_func(U,T,p,orb,sys,del,bif_ind,type,3);
% orb2 = orb_corr(orb_gr,sys2,[],bifs); % append the third event condition to MP-BVP
% % figure(); plot_orb(orb2,sys,'bldosc_p')
% 
% % Adaptive stepsize continuation run
% optsgr = br12_opts([1,2]); % position of tau and k in p0
% optsgr.stop.p_lim = [1 5; -inf inf]; % [tau_min tau_max; k_min k_max];
% optsgr.psa.ds_lim = [1e-3 0.2]; % stepsize limits [ds_min ds_max]
% branch_gr1 = br12_cont_adapt(orb2,sys2,optsgr,bifs);
% % figure(); plot_br2_par(branch_gr1,optsgr.pi,2)
% % anim_br_spectrum(branch_gr1,'test.avi')


%% Plot results
pind = [1,2];

% combined bifurcation diagram in 2D
figure();
plot_br2_par(branch1,pind); hold on
plot_br2_par(branch2,pind);
plot_br2_par(branch_gr1,pind); hold off
xlabel('$\tau$'); ylabel('$k$'); title('Continuation results')

% combined bifurcation diagram in 3D
% uind = 1;
% figure();
% plot_br2_3D(branch1,pind,uind,2); hold on
% plot_br2_3D(branch2,pind,uind,2);
% plot_br2_3D(branch_gr1,pind,uind,2); hold off
% xlabel('$\tau$'); ylabel('$k$'); zlabel('$|x|$'); title('Continuation results')

% mean kinetic energy (if sys.q is defined)
% qind = 1;
% figure();
% plot_br2_qmon(branch1,pind,qind,2); hold on
% plot_br2_qmon(branch2,pind,qind,2);
% plot_br2_qmon(branch_gr1,pind,qind,0); hold off
% xlabel('$\tau$'); ylabel('$k$'); zlabel('$||E_{\mathrm{k}}||$'); 
% title('Monitor function values')


%% Validate orbit stability via simulation
% 
% %normal orbit with crossing
% orbp = branch1(10); sim_opts.h_act = [1 1 0]; % turn off grazing event
% res = sim_stab_test(orbp,sys,sim_opts,1e-3); % Simulation started from perturbed periodic orbit
% [mu,~,~,~] = orb_stab(orbp,sys); % find Floquet multipliers
% figure(); subplot(1,2,1); plot_spectrum(mu);
% subplot(1,2,2); plot_res(res,'bldosc_p',orbp.p); hold on
% plot_orb(orbp,sys,'bldosc_p'); hold off
% 
% % grazing orbit
% orbp = branch_gr1(10); sim_opts.h_act = [0 0 0]; % turn off all events (smooth system)
% res = sim_stab_test(orbp,sys,sim_opts,1e-1); % Simulation started from perturbed periodic orbit
% [mu,~,~,~] = orb_stab(orbp,sys); % Floquet multipliers
% figure(); subplot(1,2,1); plot_spectrum(mu);
% subplot(1,2,2); plot_res(res,'bldosc_p',orbp.p); hold on
% plot_orb(orbp,sys,'bldosc_p'); hold off


%% Advanced plotting options

% Visualization of solution branches with a unique linestyle
% ls = {'b-.','LineWidth',2};
% figure(); plot_br1_ampl(branch1,1,1,0,ls);
% figure(); plot_br1_norm(branch1,1,0,ls);
% figure(); plot_br2_3D(branch1,[1 2],1,2,ls);
% figure(); plot_br2_par(branch_gr1,[1 2],0,ls);
% figure(); plot_orb(orb1,sys,'bldosc_p',100,ls)
% ms = {'rd','filled'};
% figure(); plot_br2_3D(branch_gr1,[1 2],1,0,ls); hold on
% plot_bif_points(branch_gr1,[1 2],0,4,0,ms); % mark only boundary points


%% Recreate the continuation results using COCO
% 
% % Initialize COCO
% addpath('<COCO_dir>'); % location of COCO installation
% startup;
% 
% % define a COCO compatible problem
% p_names = {'$\tau$', '$k$', '$\zeta$', '$\beta$', '$f$', '$\omega$'};
% prob1 = pwsdde_coco_prob(sys,orb0,p_names);
% 
% % Define an event for detecting vanishing segments
% prob1 = pwsdde_coco_ev_vanish(prob1,1e-2);
% 
% % COCO numeric options
% prob1 = coco_set(prob1,'cont','h_min',1e-3,'h_max',0.5);
% prob1 = coco_set(prob1,'cont','ItMX',[0 100]);
% 
% % create a branch of solutions with COCO (in tau)
% bd1 = coco(prob1,'bld_r1',[],1,{'$\tau$','mu_crit','T_min'},[3.5, 5]);
% 
% % create a second branch of solutions with COCO (in k)
% bd2 = coco(prob1,'bld_r2',[],1,{'$k$','mu_crit','T_min'},[0.1, 0.2]);
% 
% % convert to a grazing solution
% lab = coco_bd_labs('bld_r1','VA');
% orb1 = pwsdde_coco_2orb('bld_r1',prob1,lab);
% orb_gr = orb_convert(orb1,sys,0);
% % plot_orb(orb_gr,sys,'bldosc_p')
% 
% % define a new COCO compatible problem
% prob2 = pwsdde_coco_prob(sys,orb_gr,p_names);
% 
% % Include an extra zero condition to ensure grazing
% prob2 = pwsdde_coco_h_graze(prob2,sys,1); 
% 
% % create a branch of grazing solutions with COCO (in tau and k)
% bd_gr = coco(prob2,'bld_gr',[],1,{'$\tau$','$k$','mu_crit','T_min'},[1, 5]);
% 
% % Plot continuation results obtained with COCO
% thm = struct();
% thm.special = {'EP','VA'};
% thm.VA = {'dk','MarkerSize',7,'MarkerFaceColor','y'};
% figure(); box on; hold on
% coco_plot_bd(thm,'bld_r1','$\tau$','$k$','Ampl')
% coco_plot_bd(thm,'bld_r2','$\tau$','$k$','Ampl')
% coco_plot_bd(thm,'bld_gr','$\tau$','$k$','Ampl')
% hold off
% 
% % figure(); coco_plot_bd(thm,'bld_r1','$\tau$','Ampl')
% % figure(); coco_plot_bd(thm,'bld_r2','$k$','Ampl')
% 
% % Plot an individual orbit
% % lab = coco_bd_labs('bld_r1','VA');
% % orbp = pwsdde_coco_2orb('bld_r1',prob1,lab(1));
% % figure();
% % plot_orb(orbp,sys,'bldosc_p');