clear; clc; close all;
warning off backtrace
%% Example continuation problem of a nonliear robotic arm
% actively controlled, considering feedback delay and compliant coupling
% governed by a 2 DoF smooth neutral delay differential equation
% more information on the model is available in:
% https://doi.org/10.1016/j.ijnonlinmec.2022.104239

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


%% Governing NDDE
% Dimensionless equation of motion
% x1''(t) + 2*khi*r*gamma*(x1'(t)-x2'(t)) + gamma^2*r*(x1(t)-x2(t)) + 
%   x1(t) - mu*r*(x1(t)-x2(t))^2 - (k*x1''(t-tau) + knl*x1''(t-tau)^3)
% x2''(t) + 2*khi*gamma*(x2'(t)-x1'(t)) + gamma^2*(x2(t)-x1(t)) + 
%   mu*r*(x2(t)-x1(t))^2
% Dimensionless system paramters:
% 1) khi:   standalone damping ratio
% 2) r:     mass ratio
% 3) gamma: frequency ratio
% 4) mu:    quadratic nonlinear term
% 5) k:     linear feedback gain
% 6) knl:   nonlinear feedback gain
% 7) tau:   feedback delay


%% System definition
addpath('_system def/nonlinear robot arm');
sys.f = 'nlne_rob_f';       % vector field modes
sys.tau = 'nlne_rob_tau';   % time delays
sys.tau_no = 1;             % number of time delays
sys.mode_no = 1;            % number of modes

% dummy event for compatibility (zero crossing of x1')
sys.event_no = 1;   % number of possible events (one dummy)
sys.e = @(x,xd,p,id,type,l) sys_dummy_e(x,xd,p,id,type,l,2);


%% Finding a periodic orbit

% Periodic orbit initialisation via simulation
p0 = [0.05 1 1 -1 0.48 0 0.1];  % parameter vector close to a Hopf bifurcation [khi, r, gamma, mu, k, knl, tau] 
y0 = @(t) [p0(7); 0; 0; -0.1];  % initial history function
o_init.N = 2;               % number of events expected in the periodic orbit
o_init.M = 20;              % Chebyshev mesh dimension
sim_opts.t_end = 100;       % simulation time
[orb0,res]= sim_ns_dde(y0,p0,sys,o_init,sim_opts);
figure(); plot_res(res,2); title('Transient simulation result');
% figure(); plot_orb(orb0,sys,'nlne_rob_p'); title('Simulation guess');

% Correct solution guess with Newton iteration
[orb1,err] = orb_corr(orb0,sys);
figure(); plot_orb(orb1,sys,'nlne_rob_p'); title('Corrected orbit');
% figure(); plot(err);


%% Follow the periodic orbit in a system parameter

% Follow periodic orbits in k
opts1 = br12_opts(5); % index of k in p0
opts1.stop.p_lim = [0.435 inf]; % [k_min k_max]
opts1.psa.ds_lim = [1e-3 0.1]; % [ds_min ds_max]
branch1 = br12_cont_adapt(orb1,sys,opts1);
% figure(); plot_br1_norm(branch1,opts1.pi,2);
% figure(); plot_orb(branch1(end),sys,'nlne_rob_p');


%% Switch to a new branch at a symmetry breaking pitchfork bifurcation

% Initialize the starting point of a different assymetric branch
[~,~,~,~,ibif] = get_br_data(branch1,opts1.pi); % find the index of the saddle node
pert_opts.type = 2; % perturbation type (saddle node, in this case it is a pitchfork scenario)
orb2 = bif_orb_perturb(branch1(ibif.i_sc),sys,pert_opts); % perturb the orbit
corr_opts = br12_opts(5,1e-5);  % correction in k with ds = 1e-5
orb2c = orb_corr_psa(orb2,sys,corr_opts); % correct via a psa step

% Follow the new branch of periodic orbits in k
opts2 = opts1; % samme options as previously
branch2 = br12_cont_adapt(orb2c,sys,opts2);
% figure(); plot_br1_norm(branch2,opts2.pi,2);
% figure(); plot_orb(branch2(end),sys,'nlne_rob_p');
% anim_br_orb(branch2,sys,'nlne_rob_p','test.avi')


%% Plot continuation results

% Combined bifurcation diagram
pind = 5; uind = 1;
figure();
plot_br1_ampl(branch1,pind,uind,1); hold on
plot_br1_ampl(branch2,pind,uind,1); hold off
xlabel('$w$'); ylabel('$|x_1|$'); title('Continuation results')


%% Recreate the first branch using COCO
% 
% % Initialize COCO
% addpath('<COCO_dir>'); % location of COCO installation
% startup;
% 
% % define a COCO compatible problem
% p_names = {'$\xi$', '$r$', '$\gamma$', '$\mu$', '$k$', '$k_{nl}$', '$\tau$'};
% prob1 = pwsdde_coco_prob(sys,orb0,p_names);
% 
% % COCO numeric options
% prob1 = coco_set(prob1,'cont','h_min',1e-3,'h_max',0.5);
% prob1 = coco_set(prob1,'cont','ItMX',[100 100]);
% 
% % create a branch of solutions with COCO
% bd1 = coco(prob1,'nlro_r1',[],1,{'$k$','mu_crit','T_min'},[0.4, 0.5]);
% 
% % plot continuation results
% thm = struct();
% thm.special = {'EP'};
% thm.ustab = 'mu_crit';
% thm.ustabfun =  @(x) (abs(x)<1) + 1;
% thm.lspec = {{'b--', 'LineWidth', 1}, {'b-', 'LineWidth', 1}};
% figure(); box on; hold on
% coco_plot_bd(thm,'nlro_r1','$k$','Ampl')
% hold off
% 
% % Plot an individual orbit
% % lab = coco_bd_labs('nlro_r1','BP');
% % orbp = pwsdde_coco_2orb('nlro_r1',prob1,lab(2));
% % figure();
% % plot_orb(orbp,sys,'nlne_rob_p');