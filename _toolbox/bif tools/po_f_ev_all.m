function X = po_f_ev_all(y,orb,sys,pind,func,res)
%PO_F_EV_ALL Auxiliary monitor function for collocation based continuation
% tracking the values of func(x,xd,p) evaluated for all points
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
%   res: interpolation resolution (default orb.M)
% Output:
%   X: evaluation of f for all points of the orbit

% Initialization
M = orb.M;              % mesh resolution
N = length(orb.sig);    % number of segments

% Unpack solution vector
lp = length(pind);
p0 = pl_insert(orb.p,y(end-lp+1:end),pind);
Ti0 = y(end-N-lp+1:end-lp);

% Evaluate current and delayed terms
if nargin<6
    res = M;
    [~,x0] = bvp2sig(y(1:end-N-lp),Ti0,M);
    del0 = po_delay_interp(y(1:end-N-lp),Ti0,p0,M,sys);
else
    [~,x0] = bvp2sig(y(1:end-N-lp),Ti0,M,res);
    del0 = po_delay_interp(y(1:end-N-lp),Ti0,p0,M,sys,res);
end
x0_tau = del0.ud;

X = zeros(res*N,1);
for i = 1:res*N
    X(i) = func(x0(:,i),squeeze(x0_tau(:,:,i)),p0);
end

end

