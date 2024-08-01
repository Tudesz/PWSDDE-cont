function [x1,dx1,tg,err] = psa_step(x0,dx0,funcs,opts)
%PS_ARC_STEP Continuation step with the pseudo-arclength method
% Input:
%   x0: initial extended state vector [u; p]
%   dx0: guess for solution tangent
%   funcs: a function, which returns the governing NAE and its Jacobian wrt
%       state variables and continuation parameters: [F(x), JF(x)]
%   opts: continuation run options
%    -> psa: pseudo-arclength method parameters
%      -> ds: arclength stepsize
%      -> pi: indicies of continuation parameters
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
%   err: solver error at the last finished correction step

% initialization
ds = opts.psa.ds;
pind = opts.psa.pi;
lp = length(pind);

% Solution tangent vector at x0
[~,JF0] = funcs(x0);
sol_tan = @(x) deal([JF0*x; sum(x.^2) - 1], [JF0; 2*x.']);
dx1 = newton_iter(dx0,sol_tan,opts.nr);

% Predictor step
xp = x0 + ds * dx1; 
tg = [norm(xp(1:end-lp)) - norm(x0(1:end-lp)); ...
    xp(end-lp+1:end) - x0(end-lp+1:end)];

% Corrector step
corr = @(x) sol_corr(x,funcs,x0,dx1,ds);
[x1,err,~] = newton_iter(xp,corr,opts.nr);

end

% Auxiliary function for the corrector step
function [F,JF] = sol_corr(x,funcs,x0,dx,ds)
    [F0, JF0] = funcs(x); % base system without psa condition
    F = [F0; dx.'*(x-x0)-ds];
    JF = [JF0; dx.'];
end