function branch = br12_cont_adapt(orb,sys,opts,bifs)
%BR12_CONT_ADAPT Continuation of periodic orbits in an autonomous 
% hybrid/non-smooth delay differential equation using an adaptive 
% pseudo-arclength method (supports 1 and 2 parameter continuation)
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
%    -> c_logs: log iterations in initial correction steps (default true)
%    -> psa: pseudo-arclength method parameters 
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
%    -> ds: used pseudo arclength stepsize
%    -> tg: solution tangents in norm (lp+1)
%    -> U: solution vectors (M*N*n)
%    -> T: solution segment lengths (N)
%    -> p: system parameter vectors (l)
%    -> sig: solution signature (n)
%    -> n: number of degrees of freedom
%    -> M: Chebyshev mesh resolution
%    -> mu: Floquet multipliers
%    -> q: auxiliary monitor function evaluations

% check solution signature
check_sig(orb,sys);

% Solver data initialization
pind = opts.pi;     % bifurcation parameter index
lp = length(pind);  % number of continuation parameters
N = length(orb.sig);% number of segments
opts.stop.pi = pind;% for backwards compatibility
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
br = struct('bif_p',[],'bif_type',[],'error',[],'mu_crit',[],'ds',[],'tg',[],...
    'U',[],'T',[],'p',[],'sig',orb.sig,'M',orb.M,'n',orb.n,'mu',[],'q',[]);
branch = repmat(br,sum(opts.stop.n_step)+1,1);

% Correct initial guess (first continuation point) and find an approximate
% solution tangent vector
if nargin < 4
    [orb1,err,dy] = br12_init_corr(orb,sys,opts);
else
    [orb1,err,dy] = br12_init_corr(orb,sys,opts,bifs);
end

% Save the data of the first orbit
i0 = opts.stop.n_step(1)+1;
branch(i0).U = orb1.U;
branch(i0).T = orb1.T;
branch(i0).p = orb1.p;
branch(i0).bif_p = orb1.p(pind);
branch(i0).error = norm(err);
[branch(i0).mu, branch(i0).mu_crit,~] = orb_stab(branch(i0),sys);
branch(i0).bif_type = 'Starting point';

% Initialize the continuation run
y_start = [orb1.U; orb1.T; orb1.p(pind).'];  % starting point
dy_start = 1/norm(dy)*dy; % starting tangent vector
rt0 = tic;
if isfield(sys,'q') % user defined monitor functions
    branch(i0).q = feval(sys.q,y_start,orb,sys,pind);
end

% Cover both continuation directions
np = opts.stop.n_step; % maximum steps allowed
for dir = 1:2
    fprintf('\nContinuation of periodic orbits in direction %i\n',dir);
    y0 = y_start; dy = dy_start; % identical starting points
    if dir == 1
        ds = -opts.psa.ds0; % negative continuation direction
    else
        ds = opts.psa.ds0; % positive continuation direction
    end
    
    % Continuation loops
    i = 1;
    while i <= np(dir)
        if dir == 1
            ii = i0 - i;
            ii0 = i0 - i + 1;
            iim = i0 - i + 2;
        else
            ii = i0 + i; 
            ii0 = i0 + i - 1;
            iim = i0 + i - 2;
        end

        % Make a Pseudo-Arclength step
        try
            [y1,dy,tg,ds0,ds1,err] = psa_step_adapt(y0,dy,f_psa,ds,opts);
            bif_type = 0; % bifurcation flag (inactive by default)
        catch
            warning('Pseudo arclength method failed at step %i', i);
            bif_type = -2; % bifurcation flag (non convergent solution)
        end
        
        % Terminal error detection (negative segments, non convergence)
        if bif_type > -2
            [bif_type,conv,branch(ii).bif_type] = ...
                br12_term_err(i,y1,err,orb,opts);
        end

        % Terminate the continuation run on parameter boundaries
        if bif_type == 0
            if nargin < 4
                [bif_type,orb_temp,errb] = br12_bound_corr(y1,orb,...
                    opts,sys);
            else
                [bif_type,orb_temp,errb] = br12_bound_corr(y1,orb,...
                    opts,sys,bifs);
            end
            if bif_type == -1
                fprintf('   -> Parameter domain boundary reached at step %i\n',i);
                branch(ii).bif_type = sprintf('Parameter domain boundary');
                err = errb;
            end
        else
            orb_temp.U = y1(1:end-N-lp);
            orb_temp.T = y1(end-N-lp+1:end-lp); 
            orb_temp.p = orb.p;
            orb_temp.p(pind) = y1(end-lp+1:end);
        end

        % Run bifurcation detection routines
        if bif_type == 0 && i > 1
            if nargin < 4
                [bif_type,orb_temp] = bif_ev_detect(i,y0,y1,branch(ii0),...
                    sys,opts);
            else
                [bif_type,orb_temp] = bif_ev_detect(i,y0,y1,branch(ii0),...
                    sys,opts,bifs);
            end
            branch(ii).bif_type = orb_temp.bif_type;
            % Check for fold points
            if tg(2:end).'*branch(iim).tg(2:end) < 0
                branch(ii0).bif_type = [branch(ii0).bif_type ' (Fold point)'];
                fprintf('   -> Fold point detected at step %i\n',i);
            end
        elseif bif_type > -2
            % Evaluate orbit stabilty without bifurcation search
            [orb_temp.mu, orb_temp.mu_crit, ~] = orb_stab(orb_temp,sys);
        end

        % Do a bisection search for the bifurcation point if the run is
        % about to terminated
        if bif_type > 0 && conv
            try
                if nargin < 4
                    [orbb,dsb,errb] = bif_loc_bisec(branch(ii0),dy,ds0,...
                        f_psa,sys,bif_type,opts);
                else
                    [orbb,dsb,errb] = bif_loc_bisec(branch(ii0),dy,ds0,...
                        f_psa,sys,bif_type,opts,bifs);
                end
            catch
                errb = inf;
            end
            if norm(errb) > 10*opts.nr.abstol || any(orbb.T<0,'all')
                warning('Bisection search for bifurcation point failed!')
            else
                % If successful overwrite the last continuation point
                [orbb.mu, orbb.mu_crit, ~] = orb_stab(orbb,sys);
                orb_temp = orbb; err = errb; ds0 = dsb;
            end

        end

        % Evaluate the user defined monitor function if it has not been done already
        if isfield(sys,'q') && ~isfield(orb_temp,'q')
            orb_temp.q = feval(sys.q,y1,orb,sys,pind);
        end

        % Save orbit data
        if bif_type > -2
            branch(ii).U = orb_temp.U;
            branch(ii).T = orb_temp.T;
            branch(ii).p = orb_temp.p;
            branch(ii).bif_p = orb_temp.p(pind);
            branch(ii).error = norm(err);
            branch(ii0).tg = tg;
            branch(ii).ds = ds0;
            branch(ii).mu = orb_temp.mu;
            branch(ii).mu_crit = orb_temp.mu_crit;
            if isfield(sys,'q')
                branch(ii).q = orb_temp.q;
            end
        else
            np(dir) = i-1; % omit the last failed solution point
        end

        % Stop the continuation run at terminal events
        if bif_type ~= 0 && bif_type ~= -2
            np(dir) = i;
        end
            
        % Update last solution point
        y0 = y1;
        ds = ds1;
        
        % Log progress
        if mod(i,ceil(np(dir)/10))==0
            fprintf('  Step %i/%i, time: %0.3f s\n',i, np(dir), toc(rt0));
        end

        % Increase step counter
        i = i + 1; 
    end

end

% Omit failed steps (default empty data)
branch = branch(i0-np(1):i0+np(2));

fprintf('Total continuation time: %0.3f s\n',toc(rt0));

end

