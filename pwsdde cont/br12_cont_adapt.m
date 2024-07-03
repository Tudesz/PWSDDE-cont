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
%    -> psa: pseudo-arclength method parameters (optional field)
%      -> ds0: default arclength stepsize (default 0.1)
%      -> ds_lim: minimum and maximum allowed stepsize [ds_min ds_max]
%           (default [1e-3 10])
%      -> tgi_ds: initial step for obtaining a guess of the solution
%           tangent (default 1e-5*ds0)
%      -> init_corr: if true correct initial solution guess before taking 
%           any pseudo-arclength steps (default true)
%    -> stop: stopping conditions for the continuation run (optional field)
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

% Default solver options
if ~isfield(opts,'nr') || ~isfield(opts.nr,'abstol')
    opts.nr.abstol = 1e-10; % default solver tolerance
end

% Pseudo arclength stepping options
if ~isfield(opts,'psa')
    opts.psa = struct();
end
if ~isfield(opts.psa,'ds0')
    opts.psa.ds0 = 0.1; % default psa stepsize
end
if ~isfield(opts.psa,'ds_lim')
    opts.psa.ds_lim = [1e-3 10]; % limits of the psa stepsize
end
if ~isfield(opts.psa,'tgi_ds')
    opts.psa.tgi_ds = 1e-5; % default relative length of zeroth step
end
if ~isfield(opts.psa,'init_corr')
    opts.psa.init_corr = true; % correct initial guess before ps steps
end

% Stopping conditions
if ~isfield(opts,'stop')
    opts.stop = struct();
end
if ~isfield(opts.stop,'n_step')
    opts.stop.n_step = [100 100]; % number of allowed steps in both direction
end
if ~isfield(opts.stop,'p_lim')
    opts.stop.p_lim = [-inf inf]; % by default no limit on the continuation paramters
end
if ~isfield(opts.stop,'conv')
    opts.stop.conv = false; % dont stop at non-convergent solutions
end
if ~isfield(opts.stop,'Tneg')
    opts.stop.Tneg = true; % stop if negative segment lengths are encountered
end
if ~isfield(opts.stop,'stab_ch')
    opts.stop.stab_ch = false; % dont stop at stability changes
end
if ~isfield(opts.stop,'int_gr')
    opts.stop.int_gr = true; % stop at interior grazing
end
if ~isfield(opts.stop,'ext_gr')
    opts.stop.ext_gr = true; % stop at exterior grazing
end
if ~isfield(opts.stop,'slide')
    opts.stop.slide = true; % stop at sliding
end
if ~isfield(opts.stop,'user_q')
    opts.stop.user_q = false; % no stopping by default at monitor function zero crossings
end
if ~isfield(opts.stop,'va_tol')
    opts.stop.va_tol = 1e-3; % default segment length tolerance
end
opts.stop.pi = opts.pi; % for backwards compatibility

% Solver data initialization
pind = opts.pi;     % bifurcation parameter index
p0 = orb.p;         % default parameter vector
sig = orb.sig;      % solution signature
N = length(sig);    % number of smooth segments
lp = length(pind);  % number of continuation parameters

% Function set initialization
if lp == 1
    % One parameter continuation
    f = @(x,par) mpbvp(x(1:end-N),x(end-N+1:end),par,orb,sys);
    Jx = @(x,par) mpbvp_Ju(x(1:end-N),x(end-N+1:end),par,orb,sys);
    Jp = @(x,par,pi) mpbvp_Jp(x(1:end-N),x(end-N+1:end),par,orb,sys,pi);
elseif lp == 2
    % Two parameter continuation
    if nargin < 4 || ~isfield(bifs,'type') || ~isfield(bifs,'ind')
        error('Insufficient bifurcation data for two parameter continuation')
    else
        bifs.pi = pind(1);
    end
    switch bifs.type
         case 1
        % Add a grazing bifurcation condition
         f = @(x,par) [mpbvp(x(1:end-N),x(end-N+1:end),par,orb,sys);...
             mpbvp_gr(x(1:end-N),x(end-N+1:end),par,orb,sys,bifs.ind)];
         Jx = @(x,par) [mpbvp_Ju(x(1:end-N),x(end-N+1:end),par,orb,sys);...
             mpbvp_gr_Ju(x(1:end-N),x(end-N+1:end),par,orb,sys,bifs.ind)];
         Jp = @(x,par,pi) [mpbvp_Jp(x(1:end-N),x(end-N+1:end),par,orb,sys,pi);...
             mpbvp_gr_Jp(x(1:end-N),x(end-N+1:end),par,orb,sys,bifs.ind,pi)];
         case 2
        % Add a sliding bifurcation condition
         f = @(x,par) [mpbvp(x(1:end-N),x(end-N+1:end),par,orb,sys);...
             mpbvp_sl(x(1:end-N),x(end-N+1:end),par,orb,sys,bifs.ind)];
         Jx = @(x,par) [mpbvp_Ju(x(1:end-N),x(end-N+1:end),par,orb,sys);...
             mpbvp_sl_Ju(x(1:end-N),x(end-N+1:end),par,orb,sys,bifs.ind)];
         Jp = @(x,par,pi) [mpbvp_Jp(x(1:end-N),x(end-N+1:end),par,orb,sys,pi);...
             mpbvp_sl_Jp(x(1:end-N),x(end-N+1:end),par,orb,sys,bifs.ind,pi)];
         case 3
        % Add a user defined bifurcation condition
         f = @(x,par) [mpbvp(x(1:end-N),x(end-N+1:end),par,orb,sys);...
             bifs.f(x(1:end-N),x(end-N+1:end),par,orb,sys,bifs.ind,[],1)];
         Jx = @(x,par) [mpbvp_Ju(x(1:end-N),x(end-N+1:end),par,orb,sys);...
             bifs.f(x(1:end-N),x(end-N+1:end),par,orb,sys,bifs.ind,[],2)];
         Jp = @(x,par,pi) [mpbvp_Jp(x(1:end-N),x(end-N+1:end),par,orb,sys,pi);...
             bifs.f(x(1:end-N),x(end-N+1:end),par,orb,sys,bifs.ind,pi,3)];
    end
else
    error('Only one and two parameter continuation is supported');
end
funcs = struct('f',f,'Jx',Jx,'Jp',Jp); % for bif_loc_bisec.m

% Output initialization
br = struct('bif_p',[],'bif_type',[],'error',[],'mu_crit',[],'ds',[],'tg',[],...
    'U',[],'T',[],'p',[],'sig',orb.sig,'M',orb.M,'n',orb.n,'mu',[],'q',[]);
branch = repmat(br,sum(opts.stop.n_step)+1,1);

% Correct initial guess (first continuation point)
if ~isfield(opts.nr,'logs') || opts.nr.logs 
    fprintf('\nCorrect initial point\n')
end
if opts.psa.init_corr
    if lp > 1
        % At bifurcation points
        [orb1,err] = orb_corr(orb,sys,opts,bifs);
    else
        % At regular points
        [orb1, err] = orb_corr(orb,sys,opts);
    end
else
    orb1 = orb;
    err = f([orb.U; orb.T],orb.p);
end
x1 = [orb1.U; orb1.T];
i0 = opts.stop.n_step(1)+1;
branch(i0).U = orb1.U;
branch(i0).T = orb1.T;
branch(i0).p = orb1.p;
branch(i0).bif_p = orb1.p(pind);
branch(i0).error = norm(err);
[branch(i0).mu, branch(i0).mu_crit,~] = orb_stab(branch(i0),sys);
branch(i0).bif_type = 'Starting point';

% Find second point to initialize tangent in pseudo-arclength method
if ~isfield(opts,'nr') || ~isfield(opts.nr,'logs') || opts.nr.logs 
    fprintf('Find a second point\n')
end
if opts.psa.init_corr
    ds0 = opts.psa.tgi_ds*abs(opts.psa.ds0); %initial step
    if lp>1
        % In case of two parameter continuation
        orb2 = orb1;
        orb2.p(pind(2)) = orb2.p(pind(2))+ds0; % perturb parameter vector
        orb2 = orb_corr(orb2,sys,opts,bifs);
    else
        % In case of one parameter continuation
        orb2 = orb1;
        orb2.p(pind) = orb2.p(pind)+ds0; % perturb parameter vector   
        orb2 = orb_corr(orb2,sys,opts);
    end
    dy = [orb2.U;orb2.T;orb2.p(pind).']-[x1;orb1.p(pind).'];
else
    dy = [zeros(length(orb.U)+length(orb.T),1); ones(lp,1)];
end

% Initialize continuation run
y_start = [x1; orb1.p(pind).'];  % starting point
dy_start = 1/norm(dy)*dy;        % starting tangent vector
rt0 = tic;
if isfield(sys,'q')
    branch(i0).q = feval(sys.q,y_start,orb,sys,pind);
end
if lp>size(opts.stop.p_lim,1)
    p_lims = repmat(opts.stop.p_lim,2,1); % same limits used for all parameters unless specified otherwise
else
    p_lims = opts.stop.p_lim; % user defined continuation limits
end

% Cover both continuation directions
np = opts.stop.n_step; % maximum steps allowed
for dir = 1:2
    fprintf('\nContinuation of periodic orbits in direction %i\n',dir);
    y0 = y_start; dy = dy_start; % identical starting points
    if dir == 1
        ds = -opts.psa.ds0; % positive continuation direction
    else
        ds = opts.psa.ds0; % negative continuation direction
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
            [y1,dy,tg,ds0,ds1] = ps_arc_step_adapt(y0,dy,p0,f,Jx,Jp,ds,opts);
        catch
            warning('Pseudo arclength method failed at step %i', i);
            np(dir) = i-1;
            break
        end

        % Fill up a temporary orbit data structure
        p1 = pl_insert(p0,y1(end-lp+1:end),pind);
        err = norm(f(y1(1:end-lp),p1));
        orb_ii = branch(ii);
        orb_ii.U = y1(1:end-N-lp);
        orb_ii.T = y1(end-N-lp+1:end-lp);
        orb_ii.p = p1;
        bif_type = 0; % solution bifurcation flag
    
        % Warn/stop in case of nonconvergent solutions
        if err > opts.nr.abstol*1e3
            warning('Solution in br12_cont_adapt did not converge at step %i',i);
            if opts.stop.conv
                bif_type = -1;
                np(dir) = i; % omit last found orbit
            end
        end

        % Warn/stop in case of negative segment lengths
        if any(orb_ii.T<0,'all') || any(isnan(orb_ii.T),'all')
            warning('Negative segment length encountered at step %i',i);
            [~,vi] = min(orb_ii.T);
            branch(ii).bif_type = sprintf('Negative segment at index %i',vi);
            if opts.stop.Tneg
                bif_type = -1;
                np(dir) = i;  % omit last found orbit
            end
        end

        % Check for leaving the prescribed bifurcation parameter domain
        p_diff = [orb_ii.p(pind).'-p_lims(:,1);...
            -orb_ii.p(pind).' + p_lims(:,2)];
        if any(p_diff < 0,'all')
            branch(ii).bif_type = sprintf('Parameter domain boundary');
            fprintf('   -> Parameter domain boundary reached at step %i\n',i);
            % Correct final point to be exactly on the boundary
            [~,vi] = min(abs(p_diff));
            opts_temp = opts; opts_temp.nr.logs = false;
            if lp > 1
                % At bifurcation points
                bifs_temp = bifs;
                switch vi
                    case 1
                        orb_ii.p(pind(1)) = p_lims(1,1);
                        bifs_temp.pi = pind(2);
                    case 2
                        orb_ii.p(pind(2)) = p_lims(2,1);
                        bifs_temp.pi = pind(1);
                    case 3
                        orb_ii.p(pind(1)) = p_lims(1,2);
                        bifs_temp.pi = pind(2);
                    case 4
                        orb_ii.p(pind(2)) = p_lims(2,2);
                        bifs_temp.pi = pind(1);
                end
                [orb_ii,errb] = orb_corr(orb_ii,sys,opts_temp,bifs_temp);
            else
                % At regular points
                orb_ii.p(pind) = p_lims(vi);
                [orb_ii,errb] = orb_corr(orb_ii,sys,opts_temp);
            end
            err = norm(errb);
            np(dir) = i;
        end

        % Evaluate orbit stability
        if bif_type >= 0 && i < np(dir)
            [mu_ii, mu_crit_ii, ~] = orb_stab(orb_ii,sys);
            sc_cond = (1+opts.nr.abstol-abs(branch(ii0).mu_crit))*...
                (1+opts.nr.abstol-abs(mu_crit_ii))<0;
        else
            mu_ii = []; mu_crit_ii = [];
            sc_cond = false;
        end
        
        % Mark changes in stability
        if sc_cond
            if abs(imag(mu_crit_ii))>opts.nr.abstol
                branch(ii).bif_type = 'Stability change (Hopf)';
            elseif real(mu_crit_ii)>0
                branch(ii).bif_type = 'Stability change (Saddle node)';
            else
                branch(ii).bif_type = 'Stability change (Period doubling)';
            end
            fprintf('   -> Change in stability detected at step %i\n',i);
            % Stop continuation if necessary
            if opts.stop.stab_ch
                bif_type = 5;
                np(dir) = i;
            end
        end

        % Auxiliary monitor functions
        if i < np(dir) && isfield(sys,'q')
            branch(ii).q = feval(sys.q,y1,orb,sys,pind);
            if any(branch(ii).q .* branch(ii0).q<0)
                branch(ii).bif_type = 'Sign change in user defined monitor function';
                fprintf('   -> User defined zero crossing detected at step %i\n',i);
                if opts.stop.user_q
                    bif_type = 4;
                    np(dir) = i;
                end
            end
        end
        
        % Check for grazing events
        if nargin < 4
            [graze,g_int,g_ind] = po_graze_det(y0,y1,orb,sys,opts.stop);
        else
            [graze,g_int,g_ind] = po_graze_det(y0,y1,orb,sys,opts.stop,bifs);
        end
        if i < np(dir) && graze
            if g_int
                branch(ii).bif_type = sprintf('Interior grazing (in event condition %i)',g_ind);
            else
                branch(ii).bif_type = sprintf('Exterior grazing (at event %i)',g_ind);
            end
            fprintf('   -> Grazing bifurcation detected at step %i\n',i);
            if (opts.stop.int_gr && g_int) || (opts.stop.ext_gr && ~g_int)
                bif_type = 1;
                np(dir) = i;
            end
        end
        
        % Check for sliding events
        if nargin < 4
            [slide,s_ind] = po_slide_det(y1,orb,sys,opts.stop);
        else
            [slide,s_ind] = po_slide_det(y1,orb,sys,opts.stop,bifs);
        end
        if i < np(dir) && slide
            branch(ii).bif_type = sprintf('Boundary sliding (at event %i)',s_ind);
            fprintf('   -> Sliding bifurcation detected at step %i\n',i);
            if opts.stop.slide
                bif_type = 2;
                np(dir) = i;
            end
        end
        
        % Check orbits for vanishing segments
        vi = find(branch(ii).T<opts.stop.va_tol*max(branch(ii).T)); % index of vanishing segments
        if i < np(dir) && ~isempty(vi)
            branch(ii).bif_type = sprintf('Vanishing segment at index %i',vi(1));
            fprintf('   -> Vanishing segment detected at step %i\n',i);
            bif_type = 4;
            np(dir) = i;
        end

        % Check for fold points
        if i>1 % only from the second point on
            if tg(2:end).'*branch(iim).tg(2:end) < 0
                branch(ii0).bif_type = 'Fold point';
                fprintf('   -> Fold point detected at step %i\n',i);
            end
        end

        % Do a bisection search for the bifurcation point if the run is
        % about to terminated
        if np(dir) == i && bif_type > 0
            try
                if nargin < 4
                    [orbb,dsb,errb,~] = bif_loc_bisec(branch(ii0),dy,ds0,...
                        funcs,sys,bif_type,opts);
                else
                    [orbb,dsb,errb,~] = bif_loc_bisec(branch(ii0),dy,ds0,...
                        funcs,sys,bif_type,opts,bifs);
                end
            catch
                errb = inf;
            end
            if errb > 10*opts.nr.abstol || any(orbb.T<0,'all')
                warning('Bisection search for bifurcation point failed!')
            else
                % If successful overwrite the last continuation point
                orb_ii = orbb;
                err = errb;
                ds0 = dsb;
                [mu_ii, mu_crit_ii, ~] = orb_stab(orb_ii,sys);
            end

        end

        % Save orbit data
        branch(ii).U = orb_ii.U;
        branch(ii).T = orb_ii.T;
        branch(ii).p = orb_ii.p;
        branch(ii).bif_p = orb_ii.p(pind);
        branch(ii0).tg = tg;
        branch(ii).ds = ds0;
        branch(ii).error = err;
        branch(ii).mu = mu_ii;
        branch(ii).mu_crit = mu_crit_ii;
            
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

