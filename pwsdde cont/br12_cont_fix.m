function branch = br12_cont_fix(orb,sys,opts,bifs)
%BR12_CONT_FIX Continuation of periodic orbits in an automous 
% hybrid/non-smooth delay differential equation with a fixed step size 
% pseudo arclength method (1 and 2 parameter continuation are supported)
% Input:
%   orb: starting periodic orbit data structure
%    -> sig: solution signature (event list)
%    -> U: state variable vector (M*N*n)
%    -> T: segment lengths (N)
%    -> p: parameter vector
%    -> n: number of degrees of freedom
%    -> M: Chebyshev mesh resolution
%   sys: names of the functions that define the system
%    -> f: vector field and its Jacobians
%    -> e: event function, map and corresponding Jacobians
%    -> q: monitor function to track during continuation (optional)
%    -> tau: time delay and its parameter Jacobian
%    -> mode_no: number of distinct vector field modes
%    -> event_no: number of disctinc events
%    -> tau_no: number of distinct time delays
%    -> aut: if true the sytem is truly autonomous, thus the vector field 
%       condition for the final point of the orbit is included 
%       during mondoromy matrix formulation (default true)
%   opts: numerical method parameters
%    -> pi: indicies of continuation parameters (length of 1 or 2)
%    -> ds: arclength stepsize (default 0.1)
%    -> np: number of continuation steps (default 100)
%    -> c_logs: log iterations in initial correction steps (default true)
%    -> psa: pseudo-arclength method parameters 
%      -> tgi_ds: initial step for obtaining a guess of the solution
%           tangent (default 1e-5*ds0)
%      -> init_corr: if true correct initial solution guess before taking 
%           any pseudo-arclength steps (default true)
%      -> stab_eval: if true evaluate the stability of all found orbits
%           (default true, can be turned of for reduced calculation times)
%    -> stop: stopping conditions for the continuation run 
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
%    -> nr: parameters of the employed Newton iteration 
%      -> logs: if true print progress of Newton iteration (default false)
%      -> reltol: iteration stopping condition in step norm (default 1e-7)
%      -> abstol: iteration stopping condition in error norm (default 1e-10)
%      -> maxiter: maximum number of iteration steps  (default 10)
%      -> plots: if true plot progress of error function (default false)
%   bifs: extra data for two parameter continuation (optional input)
%    -> type: 1) grazing 2) sliding 3) user defined bifurcation
%    -> ind: index of bifurcation event in the solution signature 
%    -> f: function describing the user defined zero condition and 
%           its Jacobains (only required if bifs.type==3)
% Output:
%   branch: continuation run output structure
%    -> bif_p: values of bifurcation parameters (lp)
%    -> bif_type: bifurcation type if applicable, (possible types: 
%       stability change, grazing, sliding, vanishing segment, fold point)
%    -> error: error of governing system of nonlinear equations
%    -> mu_crit: critical multiplier
%    -> tg: solution tangents in norm (lp+1)
%    -> U: solution vectors (M*n)
%    -> T: solution segment lengths (N)
%    -> p: system parameter vectors (l)
%    -> sig: solution signature (n)
%    -> n: number of degrees of freedom
%    -> M: Chebyshev mesh resolution
%    -> mu: Floquet multipliers
%    -> q: auxiliary monitor function evaluations

% check solution signature
check_sig(orb,sys);

% for backwards compatibility
opts.stop.pi = opts.pi;
opts.psa.pi = opts.pi;
opts.psa.ds = opts.ds;

% Solver data initialization
sig = orb.sig;      % solution signature
N = length(sig);    % number of smooth segments
pind = opts.pi;     % indices of continuation parameters
np = opts.np;       % number of continuation steps
lp = length(pind);  % number of continuation parameters
if nargin > 3
    bifs.pi = pind(1); % match the free parameters
end

% Function set initialization
if lp == 1
    % One parameter continuation
    f_psa = @(x) orb_mpbvp(x,orb,sys,pind);
elseif lp == 2
    % Two parameter continuation
    if nargin < 4 || ~isfield(bifs,'type') || ~isfield(bifs,'ind')
        error('Insufficient bifurcation data for two parameter continuation')
    end
    f_psa = @(x) orb_mpbvp(x,orb,sys,pind,bifs);
else
    error('Only one and two parameter continuation is supported');
end

% Output initialization
br = struct('bif_p',[],'bif_type',[],'error',[],'mu_crit',[],'tg',[],...
    'U',[],'T',[],'p',orb.p,'sig',orb.sig,'M',orb.M,'n',orb.n,'mu',[],'q',[]);
branch = repmat(br,np+1,1);

% Correct initial guess (first continuation point) and find an approximate
% solution tangent vector
if nargin < 4
    [orb1,err,dy] = br12_init_corr(orb,sys,opts);
else
    [orb1,err,dy] = br12_init_corr(orb,sys,opts,bifs);
end

% Save the data of the first orbit
branch(1).U = orb1.U;
branch(1).T = orb1.T;
branch(1).p = orb1.p;
branch(1).bif_p = orb1.p(pind);
branch(1).error = norm(err);
if opts.psa.stab_eval
    [branch(1).mu, branch(1).mu_crit,~] = orb_stab(branch(1),sys);
end
branch(1).bif_type = 'Starting point';

% Initialize continuation run
y0 = [orb1.U; orb1.T; orb1.p(pind).'];  % starting point
dy = 1/norm(dy)*dy; % starting tangent vector
rt0 = tic;
if isfield(sys,'q')
    branch(1).q = feval(sys.q,y0,orb,sys,pind);
    if isfield(orb,'q')
        orb = rmfield(orb,'q');
    end
end
if lp > size(opts.stop.p_lim,1)
    p_lims = repmat(opts.stop.p_lim,2,1); % same limits used for all parameters unless specified otherwise
else
    p_lims = opts.stop.p_lim; % user defined continuation limits
end

% Continuation with Pseudo-Arclength method
fprintf('\nContinuation of periodic orbits\n')
for i = 2:np+1
    % Make a Pseudo-Arclength step
    try
        [y1,dy,tg,err] = psa_step(y0,dy,f_psa,opts);
        bif_type = 0; % bifurcation flag (inactive by default)
    catch
        warning('Pseudo arclength method failed at step %i', i);
        bif_type = -2; % bifurcation flag (non convergent solution)
    end

    % Terminal error detection (negative segments, non convergence)
    if bif_type > -2
        [bif_type,~,branch(i).bif_type] = ...
            br12_term_err(i,y1,err,orb,opts);
    end

    % Terminate the continuation run on parameter boundaries
    p_diff = [y1(end-lp+1:end) - p_lims(:,1);...
        -y1(end-lp+1:end) + p_lims(:,2)];
    if bif_type == 0 && any(p_diff < 0,'all')
        branch(i).bif_type = sprintf('Parameter domain boundary');
        fprintf('   -> Parameter domain boundary reached at step %i\n',i);
        bif_type = -1;
    end

    % Run bifurcation detection routines
    if bif_type == 0
        if nargin < 4
            [bif_type,orb_temp] = bif_ev_detect(i,y0,y1,branch(i-1),...
                sys,opts);
        else
            [bif_type,orb_temp] = bif_ev_detect(i,y0,y1,branch(i-1),...
                sys,opts,bifs);
        end
        branch(i).bif_type = orb_temp.bif_type;
        % Check for fold points
        if i > 2 && tg(2:end).'*branch(i-2).tg(2:end) < 0
            branch(i-1).bif_type = [branch(i-1).bif_type ' (Fold point)'];
            fprintf('   -> Fold point detected at step %i\n',i);
        end
    elseif bif_type > -2
        % Temporary data structure if no bifurcation detection is executed
        orb_temp = orb;
        orb_temp.U = y1(1:end-N-lp);
        orb_temp.T = y1(end-N-lp+1:end-lp); 
        orb_temp.p(pind) = y1(end-lp+1:end);
        if opts.psa.stab_eval
            [orb_temp.mu, orb_temp.mu_crit, ~] = orb_stab(orb_temp,sys);
        else
            orb_temp.mu = []; orb_temp.mu_crit = [];
        end
    end

    % Evaluate the user defined monitor function if it has not been done already
    if isfield(sys,'q') && (~isfield(orb_temp,'q') || isempty(orb_temp.q))
        orb_temp.q = feval(sys.q,y1,orb,sys,pind);
    end

    % Save orbit data
    if bif_type > -2
        branch(i).U = orb_temp.U;
        branch(i).T = orb_temp.T;
        branch(i).p = orb_temp.p;
        branch(i).bif_p = orb_temp.p(pind);
        branch(i).error = norm(err);
        branch(i-1).tg = tg;
        branch(i).mu = orb_temp.mu;
        branch(i).mu_crit = orb_temp.mu_crit;
        if isfield(sys,'q')
            branch(ii).q = orb_temp.q;
        end
    end

    % Terminate the continuation run if necessary
    if bif_type == -2
        branch = branch(1:i-1); % omit last failed step
    end
    if bif_type ~= 0
        break
    end

    % Update last solution point
    y0 = y1;
    
    % Log progress
    if mod(i-1,ceil(np/10))==0
        fprintf('  Step %i/%i, time: %0.3f s\n',i-1, np, toc(rt0));
    end
end

% Omit failed steps (default empty data)
i_max = min([length(branch) i]);
branch = branch(1:i_max);

fprintf('Total continuation time: %0.3f s\n',toc(rt0));

end

