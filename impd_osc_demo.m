clear; clc; close all;
warning off backtrace
%% COCO compatible continuation problem of 2 a DoF impact Duffing oscillator
% considering harmonic excitation and asymetric air gaps
% equation of motion shown in "impduff_def.m"
% using a nondimensionalized form and the old system definition structure

% Dependencies
addpath(genpath('_toolbox'))
addpath('pwsdde cont','plot tools')

% Default plot options
set(0, 'DefaultLineLineWidth', 1);
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultTextInterpreter','latex');
set(0,'DefaultAxesTickLabelInterpreter','latex');
set(0,'DefaultLegendInterpreter','latex');

%% System definition
addpath('_system def/two dof impact oscillator');
impduff_def;   % system definition file (returns 'struct')
sys = struct2func(struct); % conversion from old system definition structure


%% Finding a periodic orbit

% Periodic orbit initialisation based on simulation
p0 = [1.5 1.0 1.0 0.05 0.9 0.02 0.127 0.952 0.1 0.1]; % parameter set
y0 = @(t) [0; 0; 0; 0; 1; 0];   % initial conditions
o_init.N = 2;                   % number of events expected in the periodic orbit
o_init.M = 20;                  % Chebyshev mesh dimension
sim_opts.t_end = 150*pi/p0(1);  % simulation time
sim_opts.h_act = [1 1 0 0];     % only normal events are active
[orb0,res]= sim_ns_dde(y0,p0,sys,o_init,sim_opts);
figure(); plot_res(res,1); title('Transient simulation result');
% figure(); plot_orb(orb0,sys,'impduff_p'); title('Simulation guess');

% initialize a second orbit
p0(1) = 0.7; sim_opts.t_end = 150*pi/p0(1);
[orb1,res]= sim_ns_dde(y0,p0,sys,o_init,sim_opts);
% figure(); plot_res(res,1); title('Transient simulation result');
% figure(); plot_orb(orb0,sys,'impduff_p'); title('Simulation guess');


%% Follow periodic orbits in a system parameter

% Continuation run 1)
opts1.pi = 1; % index of om in p0
branch1 = br12_cont_adapt(orb0,sys,opts1);
% figure(); plot_br1_norm(branch1,opts1.pi,2);
% figure(); plot_orb(branch1(1),sys,'impduff_p');
% anim_br_orb(branch1,sys,'impduff_p')

% Continuation run 2)
opts2.pi = 1; % index of om in p0
branch2 = br12_cont_adapt(orb1,sys,opts2);
% figure(); plot_br1_norm(branch2,opts2.pi,2);
% figure(); plot_orb(branch2(end),sys,'impduff_p');
% anim_br_orb(branch2,sys,'impduff_p')


%% Switch to new branches at a symmetry breaking pitchfork bifurcations

% Perturb a detected saddle node from branch 1
[~,~,~,~,ibif] = get_br_data(branch1,opts1.pi); % indices of bifurcation points
pert_opts.type = 2; % perturbation type (saddle node)
orbpd = bif_orb_perturb(branch1(ibif.i_sc(1)),sys,pert_opts);
corr_opts.psa.ds = 1e-5; corr_opts.psa.pi = 1;
orbpdc = orb_corr_psa(orbpd,sys,corr_opts); % correct via a psa step

% Follow the new branch of periodic orbits in om
opts3 = opts1; % samme options as previously
branch3 = br12_cont_adapt(orbpdc,sys,opts3);
% figure(); plot_br1_norm(branch3,opts3.pi,2);
% figure(); plot_orb(branch3(end),sys,'impduff_p');
% anim_br_orb(branch3,sys,'impduff_p')

% Perturb a detected saddle node from branch 2
[~,~,~,~,ibif] = get_br_data(branch2,opts2.pi); % indices of bifurcation points
pert_opts.type = 2; % perturbation type (saddle node)
orbpd = bif_orb_perturb(branch2(ibif.i_sc(2)),sys,pert_opts);
corr_opts.psa.ds = 1e-5; corr_opts.psa.pi = 1;
orbpdc = orb_corr_psa(orbpd,sys,corr_opts); % correct via a psa step

% Follow the new branch of periodic orbits in om
opts4 = opts1; % samme options as previously
branch4 = br12_cont_adapt(orbpdc,sys,opts4);
% figure(); plot_br1_norm(branch4,opts4.pi,2);
% figure(); plot_orb(branch4(end),sys,'impduff_p');
% anim_br_orb(branch4,sys,'impduff_p')


%% Resume continuation with new solution signatures at grazing points

% Restart continuation from the end of branch 3
t_guess = [0.0 2.100 3.451 4.620]; % guess for the event times
orbgr = orb_convert(branch3(end),sys,t_guess); % convert the orbit
opts5 = opts3; opts5.stop.ext_gr = false; % don't stop at exterior grazing events
branch5 = br12_cont_adapt(orbgr,sys,opts5);
% figure(); plot_orb(branch5(end),sys,'impduff_p');
% anim_br_orb(branch5,sys,'impduff_p')

% Restart continuation from the end of branch 4
t_guess = [0.0 2.859 4.617 5.980]; % guess for the event times
orbgr = orb_convert(branch4(end),sys,t_guess); % convert the orbit
opts6 = opts4; opts6.stop.ext_gr = false;% don't stop at exterior grazing events
branch6 = br12_cont_adapt(orbgr,sys,opts6);
% figure(); plot_orb(branch6(end),sys,'impduff_p');
% anim_br_orb(branch6,sys,'impduff_p')


%% Plot continuation results

% combined bifurcation diagram
pind = 1;
figure();
plot_br1_ampl(branch1,pind,1); hold on
plot_br1_ampl(branch2,pind,1);
plot_br1_ampl(branch3,pind,1);
plot_br1_ampl(branch4,pind,1);
plot_br1_ampl(branch5,pind,1);
plot_br1_ampl(branch6,pind,1);
hold off




%% Recreate the first two branches with COCO
% 
% % Initialize COCO
% addpath('C:\Users\Tudo\Documents\MATLAB\COCO');
% startup;
% 
% % define a COCO compatible problem
% p_names = {'$\omega$','$\delta_1$','$\delta_2$', '$\mu$', '$r$', ...
%     '$\zeta$', '$\xi$', '$\phi$', '$\eta_1$', '$\eta_2$'};
% prob1 = pwsdde_coco_prob(sys,orb0,p_names);
% 
% % Add events for grazing detection
% prob1 = pwsdde_coco_ev_graze_int(prob1,sys);
% prob1 = pwsdde_coco_ev_graze_bd(prob1,sys);
% 
% % COCO numeric options
% prob1 = coco_set(prob1,'cont','h_min',1e-3,'h_max',5);
% prob1 = coco_set(prob1,'cont','ItMX',100);
% 
% % create a branch of solutions with COCO
% bd1 = coco(prob1,'nsd_r1',[],1,{'$\omega$','mu_crit'},[1, 5]);
% 
% % define a second COCO compatible problem
% prob2 = pwsdde_coco_prob(sys,orb1,p_names);
% prob2 = pwsdde_coco_ev_graze_int(prob2,sys);
% prob2 = pwsdde_coco_ev_graze_bd(prob2,sys);
% prob2 = coco_set(prob2,'cont','h_min',1e-3,'h_max',5);
% prob2 = coco_set(prob2,'cont','ItMX',100);
% 
% % create a second branch of solutions with COCO
% bd2 = coco(prob2,'nsd_r2',[],1,{'$\omega$','mu_crit'},[0.3, 3]);
% 
% % Plot continuation results with COCO
% 
% thm = struct();
% thm.special = {'EP','BGR','IGR'};
% thm.BGR = {'dk','MarkerSize',7,'MarkerFaceColor','r'};
% thm.IGR = {'ok','MarkerSize',7,'MarkerFaceColor','r'};
% thm.ustab = 'mu_crit';
% thm.ustabfun =  @(x) (abs(x)<1) + 1;
% thm.lspec = {{'b--', 'LineWidth', 1}, {'b-', 'LineWidth', 1}};
% figure(); box on; hold on
% coco_plot_bd(thm,'nsd_r1','$\omega$','Ampl')
% coco_plot_bd(thm,'nsd_r2','$\omega$','Ampl')
% hold off
% 
% % Plot an individual orbit
% % lab = coco_bd_labs('nsd_r1','IGR');
% % orbp = pwsdde_coco_2orb('nsd_r1',prob1,lab(1));
% % figure();
% % plot_orb(orbp,sys,'impduff_p');