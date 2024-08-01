function opts = corr_opts(logs,jac)
%CORR_OPTS Set default options for numerical BVP solutions cunducted with
% "orb_corr.m", "orb_corr_psa.m", and for using "newton_iter.m" in general 
% Input:
%   logs: log newton iteration progress (default true)
%   jac: if true use fsolve without analytic jacobians (default false)
% Output:
%   opts: numerical method parameters for using Newton's method
%    -> c_logs: log progress during orbit correction (default true)
%    -> c_jac: if true use provided Jacobian, else use fsolve (default true)
%    -> nr: parameters of the employed newton Rhapson iteration
%      -> logs: if true print progress of Newton iteration (default true)
%      -> reltol: iteration stopping condition in step norm (default 1e-7)
%      -> abstol: iteration stopping condition in error norm (default 1e-10)
%      -> maxiter: maximum number of iteration steps  (default 10)
%      -> plots: if true plot progress of error function (default false)

if nargin < 1
    opts.c_logs = true; % log by default during orbit correction
else
    opts.c_logs = logs;
end
if nargin < 2
    opts.c_jac = true;
else
    opts.c_jac = jac; % ony use fsolve if requested
end
opts.nr.logs = false;      % no logging by default during continuation
opts.nr.reltol = 1e-7;     % default for reltol>norm(x_i-x_{i-1})
opts.nr.abstol = 1e-10;    % default for abstoltol>norm(err(x_i))
opts.nr.maxiter = 10;      % defualt max iteration steps
opts.nr.plots = false;     % no plotting by default


end

