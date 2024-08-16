function Jgr = mpbvp_gr_Jp(U,T,p,orb,sys,del,gr_ind)
%MPBVP_GR_JP Grazing condition parameter Jacobian for a DDE 
%piecewise-smooth periodic orbit
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
%    -> sd_delay: if true state dependent delay definition is considered 
%       (default false)
%   del: data structure of delayed term evaluations
%    -> ud: state at t-tau(k) (n x nt x M*N*n) (ACCOUNTING FOR NEUTRAL DELAYS!)
%    -> id: segment index of t-tau(k) mapped back between 0 and T (M*n*N x nt)
%    -> fi: lagrange interpolation coefficients (M*n*N x nt x M)
%   gr_ind: index of grazing event
% Output:
%   Jgr: Jacobian of governing grazing condition wrt parameters (l x 1)

if ~isfield(sys,'sd_delay')
    sys.sd_delay = false; % by default use fix point delays
end

% Initialization
M = orb.M;              % mesh resolution
n = orb.n;              % system dimension
l = length(p);          % number of system parameters
nt = sys.tau_no;        % number of time delays
Jh = zeros(M,l);        % event condition parameter Jacobian evaluations
D0 = cheb_diff(M);      % base differentiation matrix

% Function definitions
Jhp_ej = @(j,x,xd) feval(sys.e,x,xd,p,orb.sig(j),3,0);  % event condition at ej
Jh_ej = @(j,x,xd,i) feval(sys.e,x,xd,p,orb.sig(j),2,i);  % event condition at ej
if ~sys.sd_delay
    dp_tau = @(~,i) feval(sys.tau,p,i,3); % parameter jacobian of tau_i (fixed)
else
    dp_tau = @(x,i) feval(sys.tau,x,p,i,3); % parameter jacobian of tau_i (state dependent)
end

% Find current and delayed terms
[~,us] = bvp2sig(U,T,M); % signal form of state vector
us_tau = del.ud; % interpolations of the delayed terms
i_tau = del.id; % containing segment indices
fi_tau = del.fi; % Lagrange interpolation coefficients

% Unpack relevant solution segment and evaluate event condition parameter Jacobians
ii = (gr_ind-1)*M+1:gr_ind*M; % in us and ts

for k = 1:M
    ut = us(:,ii(k));
    utau = squeeze(us_tau(:,:,ii(k)));
    Jh(k,:) = Jh(k,:) + Jhp_ej(gr_ind,ut,utau);
    
    % Effect on time delay
    for i = 1:nt
        % Interpolation coefficients
        fi_jk = squeeze(fi_tau(ii(k),i,:)); % interpolation coefficients
        dfi_jk = zeros(M);
        dfi_jk(k,:) = fi_jk.'*D0; % time derivative of interpolation coefficients
        ijk = (i_tau(ii(k),i)-1)*M*n+1:i_tau(ii(k),i)*M*n; % indicies of interpolation segment elements
        dtau_jk = -1/T(i_tau(ii(k),i))*dp_tau(ut,i);
        % Derivatives of the boundary conditions
        dh_jk = Jh_ej(gr_ind,ut,utau,i); % Jacobian of h wrt u(t-tau(k))
        Jh(k,:) = Jh(k,:) + kron(dh_jk,dfi_jk(k,:))*U(ijk)*dtau_jk;
    end
end

% Take derivative at the boundary
DM = 1/T(gr_ind) * cheb_diff(M);
Jgr = DM(end,:) * Jh;
end

