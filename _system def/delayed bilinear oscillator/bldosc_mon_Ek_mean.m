function em = bldosc_mon_Ek_mean(y,orb,sys,pind)
%BLDOSC_EK_MEAN Auxiliary monitor function for tracking the average of the
%kinetic energy in the system q = mean(1/2*v(t)^2)
% Input:
%   y: state vector with its bifurcation parameter [u0; Ti0; pi0]
%   orb: periodic orbit data structure (only the metadata part is used)
%    -> sig: sloution signature (event list)
%    -> n: number of degrees of freedom
%    -> M: Chebyshev mesh resolution
%   sys: names of the functions that define the system
%    -> f: vector field and its Jacobians
%    -> e: event function, map and corresponding Jacobians
%    -> tau: time delay and its parameter Jacobian
%   pind: index of continuation parameter
%   func: function to evaluate at all points (in x, xd and p)
% Output:
%   em: average kinetic energy in the orbit

% Initialization
M = orb.M;              % mesh resolution
N = length(orb.sig);    % number of segments

% Unpack solution vector
lp = length(pind);
U0 = y(1:end-N-lp);
Ti0 = y(end-N-lp+1:end-lp);

% Evaluate the mean kinetic energy
e = @(x) 1/2*x(2,:).^2;             % evaluate 1/2 * x'(t)^2
em = po_int(e,U0,Ti0,M)/sum(Ti0);   % average energy via integration

end

