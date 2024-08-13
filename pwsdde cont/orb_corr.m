function [orb_c,err] = orb_corr(orb,sys,opts,bifs)
%ORB_CORR Solve boundary value problem of periodic orbits with Newton
%iteration
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
%    -> tau_no: number of distinct time delays
%   opts: numerical method parameters (optional input)
%    -> c_logs: log iterations in initial correction steps (default true)
%    -> jac: if true use provided Jacobian, else use fsolve (default true)
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
%   orb_c: data structure of the corrected periodic orbit
%    -> sig: sloution signature (event list)
%    -> U: state variable vector (M*N*n)
%    -> T: segment lengths (N)
%    -> p: parameter vector
%    -> n: number of degrees of freedom
%    -> M: Chebyshev mesh resolution
%   err: MP-BVP error at its last evaluation

% check solution signature
check_sig(orb,sys);

% default solver options
if nargin < 3 || isempty(opts)
    opts = corr_opts();
end

% Set initial values and define the corresponding NAE
if nargin<4
    % Regular solution points
    x0 = [orb.U; orb.T];
    func = @(x) orb_mpbvp(x,orb,sys);
else
    % Bifurcation points
    x0 = [orb.U; orb.T; orb.p(bifs.pi)];
    func = @(x) orb_mpbvp(x,orb,sys,bifs.pi,bifs);
end

% Newton iteration
if opts.c_logs
    fprintf('\nCorrect solution guess\n')
    opts.nr.logs = true;
end
if opts.c_jac
    % exploiting the analytic Jacobian
    [x1,err] = newton_iter(x0, func, opts.nr);
else
    % relying on fsolve without a Jacobian
    if ~opts.c_logs
        fopts = optimoptions('fsolve','Display','none');
    elseif ~opts.nr.plots
        fopts = optimoptions('fsolve','Display','iter');
    else
        fopts = optimoptions('fsolve','Display','iter',...
            'PlotFcn','optimplotfval');
    end
    [x1,err] = fsolve(@(x) fsolve_trunc(x,func), x0, fopts);
end

% Output structure
orb_c = orb;
N = length(orb.sig);
if nargin<4
    orb_c.U = x1(1:end-N);
    orb_c.T = x1(end-N+1:end);
else
    orb_c.U = x1(1:end-N-1);
    orb_c.T = x1(end-N:end-1);
    orb_c.p(bifs.pi) = x1(end);
end

% Check for invalid output
if any(orb_c.T<0,'all')
    warning('Negative segment length encountered in orb_corr');
end

end

% Auxiliary function for fsolve
function err = fsolve_trunc(x,func)
    [err,~] = func(x);
end