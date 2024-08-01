function [x1,dx1,tg,ds0,ds1,err] = psa_step_adapt(x0,dx0,funcs,ds,opts)
%PSA_STEP_ADAPT Adaptive cContinuation step with the pseudo-arclength
%method
% Input:
%   x0: initial extended state vector [u; p]
%   dx0: guess for solution tangent
%   funcs: a function, which returns the governing NAE and its Jacobian wrt
%       state variables and continuation parameters: [F(x), JF(x)]
%   ds: active stepsize
%   opts: continuation run options
%    -> pi: indicies of continuation parameters
%    -> psa: pseudo-arclength method parameters
%      -> ds_lim: minimum and maximum allowed stepsize [ds_min ds_max]
%    -> nr: parameters of the employed newton Rhapson iteration
%      -> logs: if true print progress of Newton iteration (default false)
%      -> reltol: iteration stopping condition in step norm (default 1e-7)
%      -> abstol: iteration stopping condition in error norm (default 1e-10)
%      -> maxiter: maximum number of iteration steps  (default 10)
%      -> plot: if true plot progress of error function (default false)
% Output:
%   x1: new solution point [u; p]
%   dx1: corrected solution tangent
%   tg: solution tangent in solution norm
%   ds0: stepsize used (NaN if solution did not converge)
%   ds1: suggested stepsize for the next psa step
%   err: solver error at the last finished correction step

% initialization
pind = opts.pi;
lp = length(pind);
if ~isfield(opts.nr,'maxiter')
    opts.nr.maxiter = 10;
end
if abs(ds) < opts.psa.ds_lim(1)
    ds = sign(ds)*opts.psa.ds_lim(1); % apply stepsize lower bound
end
if abs(ds) > opts.psa.ds_lim(2)
    ds = sign(ds)*opts.psa.ds_lim(2); % apply stepsize upper bound
end

% Solution tangent vector at x0
[~,JF0] = funcs(x0);
sol_tan = @(x) deal([JF0*x; sum(x.^2) - 1], [JF0; 2*x.']);
dx1 = newton_iter(dx0,sol_tan,opts.nr);

% Obtaining a new solution with automatic stepsize reduction if needed
while abs(ds) >= opts.psa.ds_lim(1)

    % Predictor step
    xp = x0 + ds * dx1; 
    tg = [norm(xp(1:end-lp)) - norm(x0(1:end-lp)); ...
        xp(end-lp+1:end) - x0(end-lp+1:end)];
    
    % Corrector step
    corr = @(x) sol_corr(x,funcs,x0,dx1,ds);
    [x1,err,it] = newton_iter(xp,corr,opts.nr);

    if norm(err) > opts.nr.abstol
        ds = ds/2;
    else
        break
    end
end
ds0 = ds; % save the used timestep

% Suggest a new stepsize based on the number of used iterations
if it < opts.nr.maxiter/3
    ds1 = ds0*2; % increase stepsize
elseif it > opts.nr.maxiter/2
    ds1 = ds0/2; % decrease stepsize
else
    ds1 = ds0;
end

% Enforce limits on the new stepsize
if abs(ds1) < opts.psa.ds_lim(1)
    ds1 = sign(ds1)*opts.psa.ds_lim(1);
elseif abs(ds1) > opts.psa.ds_lim(2)
    ds1 = sign(ds1)*opts.psa.ds_lim(2);
end

end

% Auxiliary function for the corrector step
function [F,JF] = sol_corr(x,funcs,x0,dx,ds)
    [F0, JF0] = funcs(x); % base system without psa condition
    F = [F0; dx.'*(x-x0)-ds];
    JF = [JF0; dx.'];
end
