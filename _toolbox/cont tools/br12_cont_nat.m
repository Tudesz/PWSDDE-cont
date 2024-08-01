function branch = br12_cont_nat(orb,sys,opts,bifs)
%BR12_CONT_NAT Continuation of periodic orbits in an automous 
% hybrid/non-smooth delay differential equation using natural parameter 
% continuation (1 and 2 parameter continuation are supported)
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
%    -> ds: stepsize in the first continuation parameter
%    -> np: number of continuation steps
%    -> jac: if false use fsolve for mpbvp solution (default true)
%    -> c_logs: log iterations in initial correction steps (default true)
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
%   -> nr: parameters of the employed Newton iteration 
%      -> logs: if true print progress of Newton iteration (default false)
%      -> reltol: iteration stopping condition in step norm (default 1e-7)
%      -> abstol: iteration stopping condition in error norm (default 1e-10)
%      -> maxiter: maximum number of iteration steps  (default 10)
%      -> plots: if true plot progress of error function (default false)
%   bifs: extra data for two parameter continuation (optional input)
%    -> type: 1) grazing 2) sliding 3) user defined bifurcation
%    -> ind: index of bifurcation event in the solution signature 
%    -> f: function for user defined zero condition and its Jacobains 
%       (only required if type==3)
% Output:
%   branch: continuation run output structure
%    -> bif_p: values of bifurcation parameters
%    -> bif_type: bifurcation type if applicable, (possible types: 
%       stability change, grazing, sliding, vanishing segment)
%    -> error: error of governing system of nonlinear equations
%    -> mu_crit: critical multiplier
%    -> tg: dummy tangent vector
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

% Solver data initialization
pind = opts.pi;     % indices of continuation parameters
ds = opts.ds;       % used fixed stepsize
np = opts.np;       % number of continuation steps
lp = length(pind);  % number of continuation parameters
opts.stop.pi = pind;% for backwards compatibility

% Output initialization
br = struct('bif_p',[],'bif_type',[],'error',[],'mu_crit',[],'tg',[],...
    'U',[],'T',[],'p',[],'sig',orb.sig,'M',orb.M,'n',orb.n,'mu',[],'q',[]);
branch = repmat(br,np+1,1);

% Correct initial guess
if lp > 1
    % At bifurcation points
    bifs.pi = pind(2); % consider second continuation parameter as part of the state
    [orb0,err] = orb_corr(orb,sys,opts,bifs);
else
    % At regular points
    [orb0, err] = orb_corr(orb,sys,opts);
end

% Save first orbit data
branch(1).U = orb0.U;
branch(1).T = orb0.T;
branch(1).p = orb0.p;
branch(1).bif_p = orb0.p(pind);
branch(1).error = norm(err);
[branch(1).mu, branch(1).mu_crit,~] = orb_stab(branch(1),sys);

% Initialize continuation run
bif_p = orb0.p(pind(1)) + (0:np)*ds; % bifurcation parameter along the branch
rt0 = tic;
if isfield(sys,'q')
    y0 = [branch(1).U; branch(1).T; branch(1).p(pind).'];
    branch(1).q = feval(sys.q,y0,orb,sys,pind);
end
if lp>size(opts.stop.p_lim,1)
    p_lims = repmat(opts.stop.p_lim,2,1); % same limits used for all parameters unless specified otherwise
else
    p_lims = opts.stop.p_lim; % user defined continuation limits
end

% Continuation with natural parameter continuation
fprintf('\nContinuation of periodic orbits\n')
opts.c_logs = false;
for i = 2:np+1
    orb0.p(pind(1)) = bif_p(i);
    % Correct solution points
    try
        if lp > 1
            % At bifurcation points
            [orb1,err] = orb_corr(orb0,sys,opts,bifs);
        else
            % At regular points
            [orb1,err] = orb_corr(orb0,sys,opts);
        end
    catch
        branch = branch(1:i-1);
        warning('Orbit correction failed at step %i', i);
        break
    end

    % Save orbit data
    branch(i).U = orb1.U;
    branch(i).T = orb1.T;
    branch(i).p = orb1.p;
    branch(i).bif_p = orb1.p(pind);
    branch(i).error = norm(err);
    y1 = [branch(i).U; branch(i).T; branch(i).p(pind).'];
    y0 = [branch(i-1).U; branch(i-1).T; branch(i-1).p(pind).'];

    % dummy solution tangent
    norm_diff = norm([branch(i).U; branch(i).T]) ...
        - norm([branch(i-1).U; branch(i-1).T]);
    p_diff = branch(i).p(pind) - branch(i-1).p(pind);
    branch(i-1).tg = [norm_diff; p_diff.'];

    % Warn/stop in case of nonconvergent solutions
    if norm(err) > 1e2*opts.nr.abstol
        warning('Solution in br12_cont did not converge at step %i',i);
        if opts.stop.conv
            branch = branch(1:i-1);
            break
        end
    end

    % Warn/stop in case of negative segment lengths
    if any(branch(i).T<0,'all')
        warning('Negative segment length encountered at step %i',i);
        vi = find(branch(i).T<0);
        branch(i).bif_type = sprintf('Negative segment at index %i',vi(1));
        if opts.stop.Tneg
            break
        end
    end

    % Check for leaving the prescribed bifurcation parameter domain
    p_diff = [branch(i).bif_p.'-p_lims(:,1);...
        -branch(i).bif_p.' + p_lims(:,2)];
    if any(p_diff < 0,'all')
        branch(i).bif_type = sprintf('Parameter domain boundary');
        fprintf('   -> Parameter domain boundary reached at step %i\n',i);
        break
    end

    % Stop in case of NAN solutions
    if any(isnan(y1),'all')
        warning('NAN solution detected at step %i',i);
        break
    end

    % Evaluate orbit stability
    [branch(i).mu, branch(i).mu_crit,~] = orb_stab(branch(i),sys);

    % Check for changes in stability
    if (1+opts.nr.abstol-abs(branch(i-1).mu_crit))*...
            (1+opts.nr.abstol-abs(branch(i).mu_crit))<0
        if abs(imag(branch(i).mu_crit))>opts.nr.abstol
            branch(i).bif_type = 'Stability change (Hopf)';
        elseif real(branch(i).mu_crit)>0
            branch(i).bif_type = 'Stability change (Saddle node)';
        else
            branch(i).bif_type = 'Stability change (Period doubling)';
        end
        fprintf('   -> Change in stability detected at step %i\n',i);
        % Stop continuation if necessary
        if opts.stop.stab_ch
            break
        end
    end

    % Auxiliary monitor functions
    if isfield(sys,'q')
        branch(i).q = feval(sys.q,y1,orb,sys,pind);
        if any(branch(i).q.*branch(i-1).q<0)
            branch(i).bif_type = 'Sign change in user defined monitor function';
            fprintf('   -> User defined zero crossing detected at step %i\n',i);
            if opts.stop.user_q
                break % stop continuation if necessary
            end
        end
    end
    
    % Check for grazing events
    if nargin<4
        [graze,g_int,g_ind] = po_graze_det(y0,y1,orb,sys,opts.stop);
    else
        [graze,g_int,g_ind] = po_graze_det(y0,y1,orb,sys,opts.stop,bifs);
    end
    if graze
        if g_int
            branch(i).bif_type = sprintf('Interior grazing (in event condition %i)',g_ind);
        else
            branch(i).bif_type = sprintf('Exterior grazing (at event %i)',g_ind);
        end
        fprintf('   -> Grazing bifurcation detected at step %i\n',i);
        if (opts.stop.int_gr && g_int) || (opts.stop.ext_gr && ~g_int)
            break % stop continuation if necessary
        end
    end
    
    % Check for sliding events
    if nargin<4
        [slide,s_ind] = po_slide_det(y1,orb,sys,opts);
    else
        [slide,s_ind] = po_slide_det(y1,orb,sys,opts,bifs);
    end
    if slide
        branch(i).bif_type = sprintf('Boundary sliding (at event %i)',s_ind);
        fprintf('   -> Sliding bifurcation detected at step %i\n',i);
        if opts.stop.slide
            break % stop continuation if necessary
        end
    end
    
    % Check orbits for vanishing segments
    vi = find(branch(i).T<opts.stop.va_tol*max(branch(i).T)); % index of vanishing segments
    if ~isempty(vi)
        branch(i).bif_type = sprintf('Vanishing segment at index %i',vi(1));
        fprintf('   -> Vanishing segment detected at step %i\n',i);
        break
    end
        
    % Update last solution point
    orb0 = orb1;
    
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

