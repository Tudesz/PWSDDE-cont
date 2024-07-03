function [mu,mu_crit,Th,vcrit] = orb_stab(orb,sys)
%ORB_STAB Evaluate the stability of periodic orbits by formulating the
%corresponding monodromy matrix
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
%    -> aut: if true the sytem is truly autonomous, thus the vector field 
%       condition for the final point of the orbit is included 
%       during mondoromy matrix formulation (default true)
% Output:
%   mu: Floquet multipliers of the periodic orbit
%   mu_crit: critical Floquet multiplier with the highest absolute value
%   Th: monodromy matrix
%   vcrit: eigenvector corresponding to mu_crit

Th = mpbvp_mon(orb.U,orb.T,orb.p,orb,sys);
[V,D] = eig(Th); 
mu = diag(D);
[~,muc_i] = max(abs(mu));
mu_crit = mu(muc_i);
vcrit = V(:,muc_i);

end

