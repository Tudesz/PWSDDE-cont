function Jgr = mpbvp_gr_Ju(U,T,p,orb,sys,gr_ind)
%NSDDE_GR_BVP_JU Grazing condition Jacobian for a DDE piecewise-smooth 
%periodic orbit
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
%   gr_ind: index of grazing event
% Output:
%   Jgr: Jacobian of governing grazing condition (N*M*n+N x 1)

% Initialization
M = orb.M;              % mesh resolution
n = orb.n;              % system dimension
N = length(orb.sig);    % number of segments
D0 = cheb_diff(M);      % base differentiation matrix
h = zeros(M,1);         % event condition evaluations
Jh = zeros(M,N*M*n);    % event condition Jacobian evaluations
JTh = zeros(M,N);       % event condition Jacobian evaluations wrt Ti

% Function definitions
tau = feval(sys.tau,p,1); % point time delays
h_ej = @(j,x,xd) feval(sys.e,x,xd,p,orb.sig(j),1);  % event condition at ej
Jh_ej = @(j,x,xd,i) feval(sys.e,x,xd,p,orb.sig(j),2,i); % event condition at ej


% Evaluate current and delayed terms
[~,us] = bvp2sig(U,T,M); % signal form of state vector
[us_tau,i_tau,fi_tau,~,dT] = po_delay_interp(U,T,M,tau); % interpolation of delayed terms

% Unpack relevant solution segment and evaluate event condition Jacobians
ij = (gr_ind-1)*M*n+1:gr_ind*M*n; % in D and fg
ii = (gr_ind-1)*M+1:gr_ind*M; % in us and ts
DM = 1/T(gr_ind) * D0;
for k = 1:M
    % Continuous terms
    ut = us(:,ii(k));
    utau = squeeze(us_tau(:,:,ii(k)));
    h(k) = h_ej(gr_ind,ut,utau);
    Jh(k,ij(k:M:end)) = Jh(k,ij(k:M:end)) + Jh_ej(gr_ind,ut,utau,0);
    
    % Delayed terms
    for i = 1:length(tau)
        % Interpolation coefficients
        fi_jk = zeros(M);
        fi_jk(k,:) = squeeze(fi_tau(ii(k),i,:)); % interpolation coefficients
        dfi_jk = zeros(M);
        dfi_jk(k,:) = (D0.'*fi_jk(k,:).').'; % time derivative of interpolation coefficients
        ijk = (i_tau(ii(k),i)-1)*M*n+1:i_tau(ii(k),i)*M*n; % indicies of interpolation segment elements
        dT_jk = squeeze(sum(dT(ii(k),i,:,:),4)).';
        % Derivatives of the boundary conditions
        dh_jk = Jh_ej(gr_ind,ut,utau,i); % Jacobian of h wrt u(t-tau(k))
        Jh(k,ijk) = Jh(k,ijk) + kron(dh_jk,fi_jk(k,:));
        JTh(k,:) = JTh(k,:) + kron(dh_jk,dfi_jk(k,:))*U(ijk)*dT_jk;
    end
end

% Take derivative at the boundary
Jgr = DM(end,:) * [Jh JTh];

% Derivative wrt Ti
Jgr(end-N+gr_ind) = Jgr(end-N+gr_ind) - 1/T(gr_ind)*DM(end,:)*h;

end

