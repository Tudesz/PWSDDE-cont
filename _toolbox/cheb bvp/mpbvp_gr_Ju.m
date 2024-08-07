function Jgr = mpbvp_gr_Ju(U,T,p,orb,sys,del,gr_ind)
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
%    -> tau_no: number of distinct time delays
%   del: data structure of delayed term evaluations
%    -> ud: state at t-tau(k) (n x nt x M*N*n) (ACCOUNTING FOR NEUTRAL DELAYS!)
%    -> id: segment index of t-tau(k) mapped back between 0 and T (M*n*N x nt)
%    -> fi: lagrange interpolation coefficients (M*n*N x nt x M)
%    -> dT: derivative of the querry point wrt Ti without looping around 
%       (M*n*N x nt x N x n_tau)
%   gr_ind: index of grazing event
% Output:
%   Jgr: Jacobian of governing grazing condition (N*M*n+N x 1)

% Initialization
M = orb.M;              % mesh resolution
n = orb.n;              % system dimension
N = length(orb.sig);    % number of segments
nt = sys.tau_no;        % number of time delays
D0 = cheb_diff(M);      % base differentiation matrix
h = zeros(M,1);         % event condition evaluations
Jh = zeros(M,N*M*n);    % event condition Jacobian evaluations
JTh = zeros(M,N);       % event condition Jacobian evaluations wrt Ti

% Function definitions
h_ej = @(j,x,xd) feval(sys.e,x,xd,p,orb.sig(j),1,0);  % event condition at ej
Jh_ej = @(j,x,xd,i) feval(sys.e,x,xd,p,orb.sig(j),2,i); % event condition at ej
td_type = @(i) feval(sys.tau,[],i,1);                    % delay type identifier

% Find current and delayed terms
[~,us] = bvp2sig(U,T,M); % signal form of state vector
us_tau = del.ud; % interpolations of the delayed terms
i_tau = del.id; % containing segment indices
fi_tau = del.fi; % Lagrange interpolation coefficients
dT = del.dT; % derivatives wrt segment lengths

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
    for i = 1:nt

        % Interpolation
        i_k = i_tau(ii(k),i);  % index of segment containing t-tau_i
        ik = (i_k-1)*M*n+1:i_k*M*n; % indicies of interpolation segment elements
        fi_ki = squeeze(fi_tau(ii(k),i,:)).'; % interpolation coefficients
        dfi_ki = fi_ki*D0; % time derivative of interpolation coefficients
        dT_jk = squeeze(sum(dT(ii(k),i,:,:),4)).';

        % Derivatives of the event condition
        dh_jk = Jh_ej(gr_ind,ut,utau,i); % Jacobian of h wrt u(t-tau(k))
        Jh(k,ik) = Jh(k,ik) + kron(dh_jk,fi_ki);
        JTh(k,:) = JTh(k,:) + kron(dh_jk,dfi_ki)*U(ik)*dT_jk;
        if td_type(i) == 2 % extra term for neutral delays
            JTh(k,i_k) = JTh(k,i_k) - kron(dh_jk,fi_ki)*U(ik)/T(i_k);
        end
    end
end

% Take derivative at the boundary
Jgr = DM(end,:) * [Jh JTh];

% Derivative wrt Ti
Jgr(end-N+gr_ind) = Jgr(end-N+gr_ind) - 1/T(gr_ind)*DM(end,:)*h;

end

