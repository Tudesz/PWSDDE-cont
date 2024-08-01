function [orb_c,err] = orb_corr_psa(orb,sys,opts,bifs)
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

% for backwards compatibility
opts.psa.pi = opts.pi;
opts.psa.ds = opts.ds;

% Function set initialization
pind = opts.psa.pi; % continuation parameter index
lp = length(pind);  % number of continuation parameters
N = length(orb.sig);% number of smooth segments
if nargin < 4
    % Regular solution points
    func = @(x) orb_mpbvp(x,orb,sys,pind);
else
    % Bifurcation points
    func = @(x) orb_mpbvp(x,orb,sys,pind,bifs);
end

% Make a Pseudo-Arclength step
if opts.c_logs
    fprintf('\nCorrect solution guess with a psa step\n')
    opts.nr.logs = true;
end
y0 = [orb.U; orb.T; orb.p(pind).'];
dy = zeros(length(y0),1); dy(end) = 1;       
[y1,~,~,err] = psa_step(y0,dy,func,opts);

% Output structure
orb_c = orb;
orb_c.U = y1(1:end-N-lp);
orb_c.T = y1(end-N-lp+1:end-lp); 
orb_c.p(pind) = y1(end-lp+1:end);

% Check for invalid output
if any(orb_c.T<0,'all')
    warning('Negative segment length encountered in orb_corr_psa');
end

end

