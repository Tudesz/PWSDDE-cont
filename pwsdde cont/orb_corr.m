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

if nargin<3 || ~isfield(opts,'nr') || ~isfield(opts.nr,'logs')
    opts.nr.logs = true; % log iteration progress by default
end

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
        case 3
        % A user defined bifurcation condition
             f = @(x,p) [mpbvp(x(1:end-N),x(end-N+1:end),p,orb,sys);...
                 bifs.f(x(1:end-N),x(end-N+1:end),p,orb,sys,bifs.ind,[],1)];
             Jx = @(x,p) [mpbvp_Ju(x(1:end-N),x(end-N+1:end),p,orb,sys) ...
                 mpbvp_Jp(x(1:end-N),x(end-N+1:end),p,orb,sys,bifs.pi);...
                 bifs.f(x(1:end-N),x(end-N+1:end),p,orb,sys,bifs.ind,[],2) ...
                 bifs.f(x(1:end-N),x(end-N+1:end),p,orb,sys,bifs.ind,bifs.pi,3)];
    end
    x0 = [orb.U; orb.T; orb.p(bifs.pi)];
end

% Newton iteration
if opts.nr.logs && ~isfield(opts,'psa')
    fprintf('\nCorrect solution guess\n')
end
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
        [x1,err] = fsolve(@(x)f(x(1:end-1),par_p(x)),x0);
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

% Check for invalid output
if any(orb_c.T<0,'all')
    warning('Negative segment length encountered in orb_corr');
end

end

