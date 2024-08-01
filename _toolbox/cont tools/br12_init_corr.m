function [orbc,err,dy] = br12_init_corr(orb,sys,opts,bifs)
%BR12_INIT_CORR Correct initial solution guess and find a starting value
%for the solution tangent in "br12_cont_adapt.m" or "br12_cont_fix.m"
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
%   opts: numerical method parameters (optional input)
%    -> c_logs: log iterations in initial correction steps (default true)
%    -> jac: if true use provided Jacobian, else use fsolve (default true)
%    -> psa: pseudo-arclength method parameters 
%      -> ds0: default arclength stepsize (default 0.1)
%      -> tgi_ds: initial step for obtaining a guess of the solution
%           tangent (default 1e-5*ds0)
%      -> init_corr: if true correct initial solution guess before taking 
%           any pseudo-arclength steps (default true)
%    -> nr: parameters of the employed newton Rhapson iteration
%      -> logs: if true print progress of Newton iteration (default true)
%      -> reltol: iteration stopping condition in step norm (default 1e-7)
%      -> abstol: iteration stopping condition in error norm (default 1e-10)
%      -> maxiter: maximum number of iteration steps  (default 10)
%      -> plots: if true plot progress of error function (default false)
%   bifs: extra data for finding bifurcation points (optional input)
%    -> type: 1) grazing 2) sliding 3) user defined bifurcation
%    -> ind: index of the bifurcation event in the solution signature
%    -> pi: index of a free system parameter to be corrected
%    -> f: function describing the user defined zero condition and 
%       its Jacobains (only required if bifs.type==3)
% Output:
%   orbc: data structure of the corrected periodic orbit
%    -> sig: solution signature (event list)
%    -> U: state variable vector (M*N*n)
%    -> T: segment lengths (N)
%    -> p: parameter vector
%    -> n: number of degrees of freedom
%    -> M: Chebyshev mesh resolution
%   err: MP-BVP error at its last evaluation
%   dy: initial guess for the solution tangent

% Disable correction if requested
if ~opts.psa.init_corr
    opts.nr.maxiter = 1;
    opts.nr.logs = false;
    opts.c_logs = false;
end

% Correct initial guess (first continuation point)
if opts.c_logs
    fprintf('\n -> Correct initial point:')
end
if nargin > 3
    % At bifurcation points
    [orb1,err] = orb_corr(orb,sys,opts,bifs);
else
    % At regular points
    [orb1, err] = orb_corr(orb,sys,opts);
end
if opts.psa.init_corr
    orbc = orb1; % update the starting point
else
    orbc = orb; % keep the original orbit
end

% Find a second point to initialize the tangent in the pseudo-arclength method
if opts.c_logs
    fprintf('-> Find a second point:')
end
if opts.psa.init_corr
    ds0 = opts.psa.tgi_ds*abs(opts.psa.ds0); %initial step
    if length(opts.pi) > 1
        % In case of two parameter continuation
        orb2 = orb1;
        orb2.p(opts.pi(2)) = orb2.p(opts.pi(2))+ds0; % perturb parameter vector
        orb2 = orb_corr(orb2,sys,opts,bifs);
    else
        % In case of one parameter continuation
        orb2 = orb1;
        orb2.p(opts.pi) = orb2.p(opts.pi)+ds0; % perturb parameter vector   
        orb2 = orb_corr(orb2,sys,opts);
    end
    dy = [orb2.U;orb2.T;orb2.p(opts.pi).'] - ...
        [orb1.U; orb1.T;orb1.p(opts.pi).'];
else
    dy = [zeros(length(orb.U)+length(orb.T),1); ones(lp,1)];
end


end

