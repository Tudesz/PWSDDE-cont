function [x1,dx1,tg] = ps_arc_step(x0,dx0,par,f,Ju,Jp,opts)
%PS_ARC_STEP Continuation step with the pseudo-arclength method
% Input:
%   x0: initial extended state vector [u; p]
%   dx0: guess for solution tangent
%   par: parameter vector
%   f(x,p): governing system of nonlinear equations
%   Ju(x,p): Jacobian of f wrt state variables (u)
%   Jp(x,p,pi): Jacobian of f wrt bifurcation parameter
%    -> psa: pseudo-arclength method parameters
%      -> ds: arclength stepsize
%      -> pi: indicies of continuation parameters
%      -> np: number of continuation steps
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

% initialization
ds = opts.psa.ds;
pind = opts.psa.pi;
lp = length(pind);

 % inserting p in the system parameter vector
par_p = @(x) pl_insert(par,x(end-lp+1:end),pind);

% Solution tangent vector at x0
Ju0 = Ju(x0(1:end-lp),par_p(x0));
Jp0 = Jp(x0(1:end-lp),par_p(x0),pind);
dfe = @(x) [Ju0 * x(1:end-lp) + Jp0 * x(end-lp+1:end); sum(x.^2) - 1];
dJe = @(x) [Ju0, Jp0; 2*x.'];
dx1 = newton_iter(dx0,dfe,dJe,opts.nr);

% Predictor step
xp = x0 + ds * dx1; 
tg = [norm(xp(1:end-lp)) - norm(x0(1:end-lp)); ...
    xp(end-lp+1:end) - x0(end-lp+1:end)];

% Corrector step
fe = @(x) [f(x(1:end-lp),par_p(x)); dx1.'*(x-x0)-ds];
Je = @(x) [Ju(x(1:end-lp),par_p(x)), ...
    Jp(x(1:end-lp),par_p(x),pind); dx1.'];
x1 = newton_iter(xp,fe,Je,opts.nr); % new solution

end
