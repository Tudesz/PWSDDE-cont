function [x1,dx1,tg,ds0,ds1] = ps_arc_step_adapt(x0,dx0,par,f,Ju,Jp,ds,opts)
%PS_ARC_STEP_ADAPT Adaptive cContinuation step with the pseudo-arclength
%method
% Input:
%   x0: initial extended state vector [u; p]
%   dx0: guess for solution tangent
%   par: parameter vector
%   f(x,p): governing system of nonlinear equations
%   Ju(x,p): Jacobian of f wrt state variables (u)
%   Jp(x,p,pi): Jacobian of f wrt bifurcation parameter
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

 % inserting p in the system parameter vector
par_p = @(x) pl_insert(par,x(end-lp+1:end),pind);

% Solution tangent vector at x0
Ju0 = Ju(x0(1:end-lp),par_p(x0));
Jp0 = Jp(x0(1:end-lp),par_p(x0),pind);
dfe = @(x) [Ju0 * x(1:end-lp) + Jp0 * x(end-lp+1:end); sum(x.^2) - 1];
dJe = @(x) [Ju0, Jp0; 2*x.'];
dx1 = newton_iter(dx0,dfe,dJe,opts.nr);

% Obtaining a new solution with automatic stepsize reduction if needed
while abs(ds) >= opts.psa.ds_lim(1)

    % Predictor step
    xp = x0 + ds * dx1; 
    tg = [norm(xp(1:end-lp)) - norm(x0(1:end-lp)); ...
        xp(end-lp+1:end) - x0(end-lp+1:end)];
    
    % Corrector step
    fe = @(x) [f(x(1:end-lp),par_p(x)); dx1.'*(x-x0)-ds];
    Je = @(x) [Ju(x(1:end-lp),par_p(x)), ...
        Jp(x(1:end-lp),par_p(x),pind); dx1.'];
    [x1,err,it] = newton_iter(xp,fe,Je,opts.nr); % new solution

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
