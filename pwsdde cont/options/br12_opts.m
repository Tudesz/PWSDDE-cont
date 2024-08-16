function opts = br12_opts(pind,ds,np)
%BR12_OPTS Set default options for continuation runs conducted with 
% "br12_cont_adapt.m", "br12_cont_fix.m", or "br12_cont_nat.m", 
% Input:
%   pind: index of the continuation parameters (1 or 2)
%   ds: stepsize for non-adaptive continuation runs (default 0.1)
%   np: number of steps for non-adaptive continuation runs (default 100)
% Output:
%   opts: numerical method parameters for continuation runs
%    -> pi: indicies of continuation parameters (length of 1 or 2)
%    -> ds: arclength stepsize (non-adaptive only)
%    -> np: number of continuation steps (non-adaptive only)
%    -> psa: pseudo-arclength method parameters
%      -> ds0: default arclength stepsize (default 0.1)
%      -> ds_lim: minimum and maximum allowed stepsize [ds_min ds_max]
%           (default [1e-3 10])
%      -> tgi_ds: initial step for obtaining a guess of the solution
%           tangent (default 1e-5*ds0)
%      -> init_corr: if true correct initial solution guess before taking 
%           any pseudo-arclength steps (default true)
%      -> stab_eval: if true evaluate the stability of all found orbits
%           (default true, can be turned of for reduced calculation times)
%    -> stop: stopping conditions for the continuation run
%      -> n_step: maximum number of continuation steps in both directions
%           [n_step_m n_step_p] (default [100 100])
%      -> p_lim: limits of the continuation paramters p(pind) [p_min p_max]
%           can use separate limits in two parameter continuation (2x2)
%           (default [-inf inf])
%      -> conv: if true stop when Newton iteration fails to converge
%           (default true)
%      -> Tneg: stop when a negative segment lenght is encountered 
%           (default true)
%      -> stab_ch: if true stop on changes in stability (default false)
%      -> ext_gr: stop if external grazing is detected (default true)
%      -> int_gr: stop if internal grazing is detected (default true)
%      -> slide: stop if a sliding event is detected (default true)
%      -> user_q: stop on sign changes of the monitor function values 
%           (default false)
%      -> va_tol: tolerance for aborting continuation runs when segment 
%           lengths become too small (default 1e-3)
%      -> gr_tol: tolerance for grazing detection (default 1e-10)
%    -> nr: parameters of the employed Newton iteration (optional field)
%      -> logs: if true print progress of Newton iteration (default false)
%      -> reltol: iteration stopping condition in step norm (default 1e-7)
%      -> abstol: iteration stopping condition in error norm (default 1e-10)
%      -> maxiter: maximum number of iteration steps  (default 10)
%      -> plots: if true plot progress of error function (default false)

% Newton iteration options
opts = corr_opts();

% Bifurcation parameter index
opts.pi = pind;

% Pseudo arclength stepping options
opts.psa = struct();
opts.psa.ds0 = 0.1; % default psa stepsize
opts.psa.ds_lim = [1e-3 10]; % limits of the psa stepsize
opts.psa.tgi_ds = 1e-5; % default relative length of zeroth step
opts.psa.init_corr = true; % correct initial guess before ps steps
opts.psa.stab_eval = true; % evaluate the stability of the found orbits

% Stopping conditions
opts.stop.n_step = [100 100]; % number of allowed steps in both direction
opts.stop.p_lim = [-inf inf]; % by default no limit on the continuation paramters
opts.stop.conv = false; % dont stop at non-convergent solutions
opts.stop.Tneg = true; % stop if negative segment lengths are encountered
opts.stop.stab_ch = false; % dont stop at stability changes
opts.stop.int_gr = true; % stop at interior grazing
opts.stop.ext_gr = true; % stop at exterior grazing
opts.stop.slide = true; % stop at sliding
opts.stop.user_q = false; % no stopping by default at monitor function zero crossings
opts.stop.va_tol = 1e-3; % default segment length tolerance

% Additional options for non-adaptive runs
if nargin<2
    opts.ds = 0.1; % default arclength stepsize
else
    opts.ds = ds; % user defined arclength stepsize
end
if nargin<3
    opts.np = 100; % default number of continuation steps
else
    opts.np = np; % user defined number of continuation steps
end

end

