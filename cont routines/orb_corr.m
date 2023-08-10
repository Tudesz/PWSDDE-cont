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
%   opts: numerical method parameters
%    -> jac: if true use provided Jacobian, else use fsolve (default true)
%    -> nr: parameters of the employed newton Rhapson iteration
%      -> logs: if true print progress of Newton iteration (default false)
%      -> reltol: iteration stopping condition in step norm (default 1e-7)
%      -> abstol: iteration stopping condition in error norm (default 1e-10)
%      -> maxiter: maximum number of iteration steps  (default 10)
%      -> plot: if true plot progress of error function (default false)
%   bifs: extra data for finding bifurcation points (optional input)
%    -> type: 1) grazing 2) sliding
%    -> ind: index of bifurcation event in the solution signature
%    -> pi: index of a free system parameter
% Output:
%   x1: corrected state vector [u1; Ti1] (M*n+N)
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

% Function set initialization
N = length(orb.sig); % number of smooth segments

if nargin<4
    % Regular solution points
    f = @(x) mpbvp(x(1:end-N),x(end-N+1:end),orb.p,orb,sys);
    Jx = @(x) mpbvp_Ju(x(1:end-N),x(end-N+1:end),orb.p,orb,sys);
    x0 = [orb.U; orb.T];

else
    switch bifs.type
        % Grazing bifurcation points
        case 1
             f = @(x,p) [mpbvp(x(1:end-N),x(end-N+1:end),p,orb,sys);...
                 mpbvp_gr(x(1:end-N),x(end-N+1:end),p,orb,sys,bifs.ind)];
             Jx = @(x,p) [mpbvp_Ju(x(1:end-N),x(end-N+1:end),p,orb,sys) ...
                 mpbvp_Jp(x(1:end-N),x(end-N+1:end),p,orb,sys,bifs.pi);...
                 mpbvp_gr_Ju(x(1:end-N),x(end-N+1:end),p,orb,sys,bifs.ind) ...
                 mpbvp_gr_Jp(x(1:end-N),x(end-N+1:end),p,orb,sys,bifs.ind,bifs.pi)];       
         % Sliding bifurcation point
        case 2
             f = @(x,p) [mpbvp(x(1:end-N),x(end-N+1:end),p,orb,sys);...
                 mpbvp_sl(x(1:end-N),x(end-N+1:end),p,orb,sys,bifs.ind)];
             Jx = @(x,p) [mpbvp_Ju(x(1:end-N),x(end-N+1:end),p,orb,sys) ...
                 mpbvp_Jp(x(1:end-N),x(end-N+1:end),p,orb,sys,bifs.pi);...
                 mpbvp_sl_Ju(x(1:end-N),x(end-N+1:end),p,orb,sys,bifs.ind) ...
                 mpbvp_sl_Jp(x(1:end-N),x(end-N+1:end),p,orb,sys,bifs.ind,bifs.pi)];
    end
    x0 = [orb.U; orb.T; orb.p(bifs.pi)];
end


% Newton iteration
if nargin<4
    % Simple solution point
    if ~isfield(opts,'jac') || opts.jac
        [x1,err] = newton_iter(x0, f, Jx, opts.nr);
    else
        [x1,err] = fsolve(f, x0);
    end
else
    % Bifurcation point
    par_p = @(x) pl_insert(orb.p,x(end),bifs.pi);
    if ~isfield(opts,'jac') || opts.jac   
        [x1,err] = newton_iter(x0, @(x)f(x(1:end-1),par_p(x)),...
            @(x)Jx(x(1:end-1),par_p(x)),opts.nr);
    else
        [x1,err] = fsolve(@(x)f(x,par_p(x)),x0);
    end
end

% Output structure
orb_c = orb;
if nargin<4
    orb_c.U = x1(1:end-N);
    orb_c.T = x1(end-N+1:end);
else
    orb_c.U = x1(1:end-N-1);
    orb_c.T = x1(end-N:end-1);
    orb_c.p(bifs.pi) = x1(end);
end

end

