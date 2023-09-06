function branch = br12_cont(orb,sys,opts,bifs)
%BR12_CONT Continuation of periodic orbits in an automous hybrid/non-smooth
% delayd differential equations (1 and 2 parameter continuation is
% supported)
% Input:
%   orb: starting periodic orbit data structure
%    -> sig: sloution signature (event list)
%    -> U: state variable vector (M*N*n)
%    -> T: segment lengths (N)
%    -> p: parameter vector
%    -> n: number of degrees of freedom
%    -> M: Chebyshev mesh resolution
%   sys: names of the functions that define the system
%    -> f: vector field and its Jacobians
%    -> e: event function, map and corresponding Jacobians
%    -> tau: time delay and its parameter Jacobian
%    -> mode_no: number of distinct vector field modes
%    -> event_no: number of disctinc events
%   opts: numerical method parameters
%    -> psa: pseudo-arclength method parameters
%      -> ds: arclength stepsize
%      -> pi: indicies of continuation parameters
%      -> np: number of continuation steps
%      -> sc_stop: if true stop on changes in stability (default false)
%      -> egr_stop: stop if external grazing is detected (default true)
%      -> igr_stop: stop if internal grazing is detected (default true)
%      -> sl_stop: stop if a sliding event is detected (default true)
%      -> va_tol: tolerance for aborting continuation runs when segment 
%           lengths become too small (default 1e-3)
%      -> gr_tol: tolerance for grazing detection (default 1e-10)
%      -> init_corr: if true correct initial solution guess before taking 
%           any pseudo-arclength steps (default true)
%    -> nr: parameters of the employed newton Rhapson iteration
%      -> logs: if true print progress of Newton iteration (default false)
%      -> reltol: iteration stopping condition in step norm (default 1e-7)
%      -> abstol: iteration stopping condition in error norm (default 1e-10)
%      -> maxiter: maximum number of iteration steps  (default 10)
%      -> plot: if true plot progress of error function (default false)
%   bifs: extra data for two parameter continuation (optional input)
%    -> type: 1) grazing 2) sliding
%    -> ind: index of bifurcation event in the solution signature
% Output:
%   branch: continuation run output structure
%    -> U: solution vectors (M*n)
%    -> T: solution segment lengths (N)
%    -> p: system parameter vectors (l)
%    -> sig: solution signature (n)
%    -> n: number of degrees of freedom
%    -> M: Chebyshev mesh resolution
%    -> tg: solution tangents in norm
%    -> mu: Floquet multipliers
%    -> mu_crit: critical multiplier
%    -> error: error of governing system of nonlinear equations
%    -> bif_p: values of bifurcation parameters
%    -> bif_type: bifurcation type if applicable, (possible types: 
%       stability change, grazing, sliding, vanishing segment)

% cehck solution signature
check_sig(orb,sys);

fprintf('\nContinuation of periodic orbits\n')

% Solver data initialization
p0 = orb.p;         % default parameter vector
sig = orb.sig;      % solution signature
N = length(sig);    % number of smooth segments
pind = opts.psa.pi; % indices of continuation parameters
np = opts.psa.np;   % number of continuation steps
lp = length(pind);  % number of continuation parameters
if ~isfield(opts.psa,'init_corr')
    opts.psa.init_corr = true; % correct initial guess before ps steps
end
if ~isfield(opts.psa,'sc_stop')
    opts.psa.sc_stop = false; % dont stop at stability changes
end
if ~isfield(opts.psa,'igr_stop')
    opts.psa.igr_stop = true; % stop at interior grazing
end
if ~isfield(opts.psa,'egr_stop')
    opts.psa.egr_stop = true; % stop at exterior grazing
end
if ~isfield(opts.psa,'sl_stop')
    opts.psa.sl_stop = true; % stop at sliding
end
if ~isfield(opts.psa,'va_tol')
    opts.psa.va_tol = 1e-3; % default segment length tolerance
end
if ~isfield(opts,'nr') || ~isfield(opts.nr,'abs_tol')
    opts.nr.abstol = 1e-10; % default solver tolerance
end

% Function set initialization
 Th = @(x,par) mpbvp_mon(x(1:end-N),x(end-N+1:end),par,orb,sys);
 if lp == 1
    % One parameter continuation
    f = @(x,par) mpbvp(x(1:end-N),x(end-N+1:end),par,orb,sys);
    Jx = @(x,par) mpbvp_Ju(x(1:end-N),x(end-N+1:end),par,orb,sys);
    Jp = @(x,par,pi) mpbvp_Jp(x(1:end-N),x(end-N+1:end),par,orb,sys,pi);

elseif lp == 2
    % Two parameter continuation
    if nargin<4 || ~isfield(bifs,'type') || ~isfield(bifs,'ind')
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
    end
else
    error('Only one and two parameter continuation is supported');
 end

% Output initialization
br.bif_p = NaN(1,length(pind));
br.mu_crit = NaN; 
br.U = NaN(length(orb.U),1);
br.T = NaN(N,1);
br.p = NaN(size(p0));
br.sig = orb.sig;
br.M = orb.M;
br.n = orb.n;
br.mu = NaN; 
br.tg = NaN(2,1);
br.error = NaN;
br.bif_type = [];
branch = repmat(br,np+1,1);

% Correct initial guess (first continuation point)
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
branch(1).U = orb1.U;
branch(1).T = orb1.T;
branch(1).p = orb1.p;
branch(1).bif_p = orb1.p(pind);
branch(1).error = norm(err);
branch(1).mu = eig(Th(x1,orb1.p));
[~,muc_i] = max(abs(branch(1).mu));
branch(1).mu_crit = branch(1).mu(muc_i);

% Find second point to initialize tangent in pseudo-arclength method
if opts.psa.init_corr
    ds0 = 1e-3*abs(opts.psa.ds); %initial step
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
    dy = [zeros(length(orb.U)+length(orb.T),1);1];
end

% Continuation with Pseudo-Arclength method
y0 = [x1; orb1.p(pind).'];  % starting point
dy = 1/norm(dy)*dy;         % starting tangent vector
rt0 = tic;
for i = 2:np+1
    [y1,dy,tg] = ps_arc_step(y0,dy,p0,f,Jx,Jp,opts);
    branch(i).U = y1(1:end-N-lp);
    branch(i).T = y1(end-N-lp+1:end-lp);
    branch(i).p = pl_insert(p0,y1(end-lp+1:end),pind);
    branch(i).bif_p = y1(end-lp+1:end);
    branch(i-1).tg = tg;
    branch(i).error = norm(f(y1(1:end-lp),branch(i).p));
    branch(i).mu = eig(Th(y1(1:end-lp),branch(i).p));
    [~,muc_i] = max(abs(branch(i).mu));
    branch(i).mu_crit = branch(i).mu(muc_i);

    % Warn in case of negative segment lengths
    if any(branch(i).T<0,'all')
        warning('Negative segment length encountered at step %i',i);
        % break
    end

    % Warn in case of nonconvergent solutions
    if branch(i).error>1e2*opts.nr.abstol
        warning('Solution in br12_cont did not converge at step %i',i);
        % break
    end
    
    % Check for grazing events
    if nargin<4
        [graze,g_int,g_ind] = po_graze_det(y0,y1,orb,sys,opts.psa);
    else
        [graze,g_int,g_ind] = po_graze_det(y0,y1,orb,sys,opts.psa,bifs);
    end
    if graze
        if g_int
            branch(i).bif_type = sprintf('Interior grazing sig(%i)',g_ind);
        else
            branch(i).bif_type = sprintf('Exterior grazing h(%i)',g_ind);
        end
        % Stop continuation if neccessary
        if (opts.psa.igr_stop && g_int) || (opts.psa.egr_stop && ~g_int)
            fprintf('   -> Grazing bifurcation detected at step %i\n',i);
            break
        end
    end
    
    % Check for sliding events
    if nargin<4
        [slide,s_ind] = po_slide_det(y1,orb,sys,opts.psa);
    else
        [slide,s_ind] = po_slide_det(y1,orb,sys,opts.psa,bifs);
    end
    if slide
        branch(i).bif_type = sprintf('Boundary sliding sig(%i)',s_ind);
        % Stop continuation if neccessary
        if opts.psa.sl_stop
            fprintf('   -> Sliding bifurcation detected at step %i\n',i);
            break
        end
    end
    
    % Check for changes in stability
    if (1+opts.nr.abstol-abs(branch(i-1).mu_crit))*...
            (1+opts.nr.abstol-abs(branch(i).mu_crit))<0
        branch(i).bif_type = 'Stability change';
        % Stop continuation if neccessary
        if opts.psa.sc_stop
            fprintf('   -> Change in stability detected at step %i\n',i);
            break
        end
    end
    
    % Check orbits for vanishing segments
    vi = find(branch(i).T<opts.psa.va_tol*max(branch(i).T)); % index of vanishing segments
    if ~isempty(vi)
        branch(i).bif_type = "Vanishing segment at index " + vi(1);
        % Stop continuation if neccessary
        fprintf('   -> Vanishing segment detected at step %i\n',i);
        break
    end
        
    % Update last solution point
    y0 = y1;
    
    % Log progress
    if mod(i-1,ceil(np/10))==0
        fprintf('  Step %i/%i, time: %0.3f s\n',i-1, np, toc(rt0));
    end
end

% Omit failed steps (default NaN data)
branch = branch(1:i);

fprintf('Total continuation time: %0.3f s\n',toc(rt0));

end

