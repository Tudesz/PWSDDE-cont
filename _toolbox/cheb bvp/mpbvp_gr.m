function gr = mpbvp_gr(U,T,p,orb,sys,del,gr_ind)
%MPBVP_GR Grazing condition for a DDE piecewise-smooth periodic orbit
% Input:
%   U: free state variables (N*M*n)
%   T: free segment lengths (N)
%   p: parameter vector
%   orb: periodic orbit data structure (only the metadata part is used)
%    -> sig: sloution signature (event list)
%    -> n: number of degrees of freedom
%    -> M: Chebyshev mesh resolution
%   sys: names of the functions that define the system
%    -> f: vector field and its Jacobians
%    -> e: event function, map and corresponding Jacobians
%    -> tau: time delay and its parameter Jacobian
%    -> tau_no: number of distinct time delays
%   del: data structure of delayed term evaluations
%    -> ud: state at t-tau(k) (n x nt x M*N*n) (ACCOUNTING FOR NEUTRAL DELAYS!)
%   gr_ind: index of grazing event
% Output:
%   gr: error of governing grazing condition (1)

% Initialization
M = orb.M;              % mesh resolution
h = zeros(M,1);         % event condition evaluations

% Function definitions
h_ej = @(j,x,xd) feval(sys.e,x,xd,p,orb.sig(j),1,0);  % event condition at ej

% Find current and delayed terms
[~,us] = bvp2sig(U,T,M); % signal form of state vector
us_tau = del.ud; % extract the interpolations of delayed terms

% Unpack relevant solution segment and evaluate event conditions
ii = (gr_ind-1)*M+1:gr_ind*M; % in us and ts
for k = 1:M
    ut = us(:,ii(k));
    utau = squeeze(us_tau(:,:,ii(k)));
    h(k) = h_ej(gr_ind,ut,utau);
end

% Take derivative of h at the boundary
DM = 1/T(gr_ind) * cheb_diff(M);
gr = DM(end,:)*h;

end

