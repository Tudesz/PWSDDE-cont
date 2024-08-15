function [type,orb,err] = br12_bound_corr(y0,orb0,opts,sys,bifs)
%BR12_BOUND_CORR Detect parameter boundary crossings and correct the last
%point to be just on the boundary in the output of % "br12_cont_adapt.m",
% "br12_cont_fix.m", or "br12_cont_nat.m"
% Input:
%   y0: previous state vector with its bifurcation parameter [u0; Ti0; pi0]
%   orb0: default periodic orbit data structure
%    -> sig: solution signature (event list)
%    -> U: solution vectors (M*N*n)
%    -> T: solution segment lengths (N)
%    -> p: system parameter vectors (l)
%   sys: names of the functions that define the system
%    -> f: vector field and its Jacobians
%    -> e: event function, map and corresponding Jacobians
%    -> q: monitor function to track during continuation (optional)
%    -> tau: time delay and its parameter Jacobian
%    -> mode_no: number of distinct vector field modes
%    -> event_no: number of disctinc events
%    -> tau_no: number of distinct time delays
%   opts: numerical method parameters
%    -> pi: indicies of continuation parameters (length of 1 or 2)
%    -> psa: pseudo-arclength method parameters 
%    -> c_logs: log iterations in initial correction steps (default true)
%      -> ds0: default arclength stepsize (default 0.1)
%      -> ds_lim: minimum and maximum allowed stepsize [ds_min ds_max]
%           (default [1e-3 10])
%      -> tgi_ds: initial step for obtaining a guess of the solution
%           tangent (default 1e-5*ds0)
%      -> init_corr: if true correct initial solution guess before taking 
%           any pseudo-arclength steps (default true)
%    -> stop: stopping conditions for the continuation run 
%      -> n_step: maximum number of continuation steps in both directions
%           [n_step_m n_step_p] (default [100 100])
%      -> p_lim: limits of the continuation paramters p(pind) [p_min p_max]
%           can use separate limits in two parameter continuation (2x2)
%           (default [-inf inf])
%    -> nr: parameters of the employed Newton iteration 
%      -> abstol: iteration stopping condition in error norm (default 1e-10)
%   bifs: extra data for two parameter continuation (optional input)
%    -> type: 1) grazing 2) sliding
%    -> ind: index of bifurcation event in the solution signature

% No logging by default
opts.c_logs = false;
opts.nr.logs = false;

% Unpack orbit data structure
N = length(orb0.sig);   % number of segments
pind = opts.pi;         % bifurcation parameter index
lp = length(pind);      % number of continuation parameters
type = 0;               % solution bifurcation flag (empty by default)

% Set parameter domain boundaries
if lp > size(opts.stop.p_lim,1)
    p_lims = repmat(opts.stop.p_lim,2,1); % same limits used for all parameters unless specified otherwise
else
    p_lims = opts.stop.p_lim; % user defined continuation limits
end
err = [];

% Fill up a temporary orbit data structure
orb = orb0;
orb.U = y0(1:end-N-lp);
orb.T = y0(end-N-lp+1:end-lp); 
orb.p(pind) = y0(end-lp+1:end);

% Check for leaving the prescribed bifurcation parameter domain
p_diff = [orb.p(pind).'-p_lims(:,1); -orb.p(pind).' + p_lims(:,2)];
if any(p_diff < 0,'all')
    type = -1; % bifurcation flag for parameter domain boundary crossing
    % Correct final point to be exactly on the boundary
    [~,vi] = min(p_diff);
    if lp > 1
        % At bifurcation points
        bifs_temp = bifs;
        switch vi
            case 1
                orb.p(pind(1)) = p_lims(1,1);
                bifs_temp.pi = pind(2);
            case 2
                orb.p(pind(2)) = p_lims(2,1);
                bifs_temp.pi = pind(1);
            case 3
                orb.p(pind(1)) = p_lims(1,2);
                bifs_temp.pi = pind(2);
            case 4
                orb.p(pind(2)) = p_lims(2,2);
                bifs_temp.pi = pind(1);
        end
        try
            [orb,err] = orb_corr(orb,sys,opts,bifs_temp);
        catch
            warning('Search for a bounary point failed in br12_bound_corr.m!');
            type = -2; % flag for failed solutions
        end
    else
        % At regular points
        orb.p(pind) = p_lims(vi);
        try
            [orb,err] = orb_corr(orb,sys,opts);
        catch
            warning('Search for a bounary point failed in br12_bound_corr.m!');
            type = -2; % flag for failed solutions
        end
    end
end

% Evaluate the stability of the new orbit
if type > -2 && norm(err) < opts.nr.abstol*1e3
    [orb.mu, orb.mu_crit, ~] = orb_stab(orb,sys);
else
    orb.mu = []; orb.mu_crit = [];
end

end

