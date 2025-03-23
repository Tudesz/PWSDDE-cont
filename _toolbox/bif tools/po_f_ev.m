function X = po_f_ev(y,orb,sys,pind,func)
%PO_F_EV Auxiliary monitor function for collocation based continuation
% tracking the values of func(x,xd,p)
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
%   X: evaluation of f for all points of the orbit

% Initialization
M = orb.M;              % mesh resolution
N = length(orb.sig);    % number of segments

% Unpack solution vector
lp = length(pind);
p0 = orb.p; 
p0(pind) = y(end-lp+1:end);
Ti0 = y(end-N-lp+1:end-lp);

% Evaluate current and delayed terms
[~,x0] = bvp2sig(y(1:end-N-lp),Ti0,M);
del0 = po_delay_interp(y(1:end-N-lp),Ti0,p0,M,sys);
x0_tau = del0.ud;

% Evaluate monitor function
X = func(x0,x0_tau,p0);

end

