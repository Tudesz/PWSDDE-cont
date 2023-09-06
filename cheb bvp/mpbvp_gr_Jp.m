function Jgr = mpbvp_gr_Jp(U,T,p,orb,sys,gr_ind,lp)
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
%   gr_ind: index of grazing event
%   lp: bifurcation parameter index (default all)
% Output:
%   Jgr: Jacobian of governing grazing condition wrt parameters (l x 1)

% Initialization
M = orb.M;              % mesh resolution
n = orb.n;              % system dimension
l = length(p);          % number of system parameters
if nargin<7
    lp = 1:l;           % indicies of bifurcation parameters
end
Jh = zeros(M,l);        % event condition parameter Jacobian evaluations
D0 = cheb_diff(M);      % base differentiation matrix

% Function definitions
tau = feval(sys.tau,p,1); % point time delays
dtau = feval(sys.tau,p,2); % derivatives of point time delays
Jhp_ej = @(j,x,xd) feval(sys.e,x,xd,p,orb.sig(j),3,0);  % event condition at ej
Jh_ej = @(j,x,xd,i) feval(sys.e,x,xd,p,orb.sig(j),2,i);  % event condition at ej


% Evaluate current and delayed terms
[~,us] = bvp2sig(U,T,M); % signal form of state vector
[us_tau,i_tau,fi_tau,~,~] = po_delay_interp(U,T,M,tau); % interpolation of delayed terms

% Unpack relevant solution segment and evaluate event condition parameter Jacobians
ii = (gr_ind-1)*M+1:gr_ind*M; % in us and ts

for k = 1:M
    ut = us(:,ii(k));
    utau = squeeze(us_tau(:,:,ii(k)));
    Jh(k,:) = Jh(k,:) + Jhp_ej(gr_ind,ut,utau);
    
    % Effect on time delay
    for i = 1:length(tau)
        % Interpolation coefficients
        fi_jk = squeeze(fi_tau(ii(k),i,:)); % interpolation coefficients
        dfi_jk = zeros(M);
        dfi_jk(k,:) = fi_jk.'*D0; % time derivative of interpolation coefficients
        ijk = (i_tau(ii(k),i)-1)*M*n+1:i_tau(ii(k),i)*M*n; % indicies of interpolation segment elements
        dtau_jk = -1/T(i_tau(ii(k),i))*dtau(i,:);
        % Derivatives of the boundary conditions
        dh_jk = Jh_ej(gr_ind,ut,utau,i); % Jacobian of h wrt u(t-tau(k))
        Jh(k,:) = Jh(k,:) + kron(dh_jk,dfi_jk(k,:))*U(ijk)*dtau_jk;
    end
end

% Take derivative at the boundary
DM = 1/T(gr_ind) * cheb_diff(M);
Jgr = DM(end,:) * Jh(:,lp);
end

