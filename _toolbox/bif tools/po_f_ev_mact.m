function X = po_f_ev_mact(y,orb,sys,pind,func,mi,res)
%PO_F_EV_MACT Auxiliary monitor function for collocation based continuation
% tracking the values of func(x,xd,p) evaluated for all points in modes mi
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
%   mi: indices of active modes to evaluate func(x,xd,p) in
%   res: interpolation resolution (default orb.M)
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
if nargin<7
    res = M;
    [~,x0] = bvp2sig(y(1:end-N-lp),Ti0,M);
    del0 = po_delay_interp(y(1:end-N-lp),Ti0,p0,M,sys);
else
    [~,x0] = bvp2sig(y(1:end-N-lp),Ti0,M,res);
    del0 = po_delay_interp(y(1:end-N-lp),Ti0,p0,M,sys,res);
end
x0_tau = del0.ud;

% Evaluate user defined function at active points
X = zeros(res*N,1);
for i = 1:N
    ii = (i-1)*res+1:i*res;
    m_act = feval(sys.e,[],[],[],orb.sig(i),7,1);  % active mode of dde rhs
    if any(m_act == mi,'all')
        for j = 1:res
            X(ii(j)) = func(x0(:,ii(j)),squeeze(x0_tau(:,:,ii(j))),p0);
        end
    else
        X(ii) = NaN;
    end
end
X(isnan(X)) = [];

end

