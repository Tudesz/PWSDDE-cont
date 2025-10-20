function X = po_f_ev_bound(y,orb,sys,pind,func,ei)
%PO_F_EV_BOUND Auxiliary monitor function for collocation based continuation
% tracking the values of func(x,xd,p) evaluated at boundary points
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
%   rei: indicies of events to evaluate the monitor function at 
%       (default all)
% Output:
%   X: evaluation of f for selected boundary points of the orbit

% Initialization
M = orb.M;              % mesh resolution
N = length(orb.sig);    % number of segments
if nargin < 6
    ei = 1:N; % all boundaries monitored by default
end

% Unpack solution vector
lp = length(pind);
p0 = orb.p; 
p0(pind) = y(end-lp+1:end);
Ti0 = y(end-N-lp+1:end-lp);

% Evaluate current and delayed terms
[~,x0] = bvp2sig(y(1:end-N-lp),Ti0,M);
del0 = po_delay_interp(y(1:end-N-lp),Ti0,p0,M,sys);
x0_tau = del0.ud;

% Evaluate user defined function at selected boundaries
X = zeros(length(ei),1);
for i = 1:length(ei)
    X(i) = func(x0(:,ei(i)*M),squeeze(x0_tau(:,:,ei(i)*M)),p0);
end

end

