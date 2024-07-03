function [orb_c,err] = orb_corr_psa(orb,sys,opts)
%ORB_CORR_PSA: Correct solution guess with a single pseudo-archlength step
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
%   opts: numerical method parameters
%    -> psa: pseudo-arclength method parameters
%      -> ds: arclength stepsize
%      -> pi: indicies of continuation parameters
%    -> nr: parameters of the employed newton Rhapson iteration
%      -> logs: if true print progress of Newton iteration (default true)
%      -> reltol: iteration stopping condition in step norm (default 1e-7)
%      -> abstol: iteration stopping condition in error norm (default 1e-10)
%      -> maxiter: maximum number of iteration steps  (default 10)
%      -> plots: if true plot progress of error function (default false) 
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
pind = opts.psa.pi;     % continuation parameter index
N = length(orb.sig);    % number of smooth segments
f = @(x,par) mpbvp(x(1:end-N),x(end-N+1:end),par,orb,sys);
Jx = @(x,par) mpbvp_Ju(x(1:end-N),x(end-N+1:end),par,orb,sys);
Jp = @(x,par,pi) mpbvp_Jp(x(1:end-N),x(end-N+1:end),par,orb,sys,pi);

% Make a Pseudo-Arclength step
y0 = [orb.U; orb.T; orb.p(pind).'];
dy = zeros(length(y0),1); dy(end) = 1;       
[y1,~,~] = ps_arc_step(y0,dy,orb.p,f,Jx,Jp,opts);

% Output structure
orb_c = orb;
orb_c.U = y1(1:end-N-1);
orb_c.T = y1(end-N-1+1:end-1);
orb_c.p = pl_insert(orb.p,y1(end),pind);
err = norm(f(y1(1:end-1),orb_c.p));

% Check for invalid output
if any(orb_c.T<0,'all')
    warning('Negative segment length encountered in orb_corr_psa');
end

end

