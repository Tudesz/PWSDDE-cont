function Jsl = mpbvp_sl_Jp(U,T,p,orb,sys,del,sl_ind)
%MPBVP_SL_JP Sliding condition Jacobian for a DDE piecewise-smooth 
%periodic orbit
% IMPORTANT: in the current version, these assumptions must hold:
%   -> dh/dx is independet of p
%   -> g(x) = x at the sliding event
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
% Output:
%   Jsl: Jacobian of governing sliding condition wrt parameters (l x 1)

if ~isfield(sys,'sd_delay')
    sys.sd_delay = false; % by default use fix point delays
end

% Initialization
M = orb.M;              % mesh resolution
n = orb.n;              % system dimension
l = length(p);          % number of system parameters
nt = sys.tau_no;        % number of time delays
D0 = cheb_diff(M);      % base differentiation matrix

% Function definitions
Jf_mj = @(mj,x,xd,i) feval(sys.f,x,xd,p,mj,2,i); % vector field Jacobian in mode mj
if ~sys.sd_delay
    dp_tau = @(~,i) feval(sys.tau,p,i,3); % parameter jacobian of tau_i (fixed)
else
    dp_tau = @(x,i) feval(sys.tau,x,p,i,3); % parameter jacobian of tau_i (state dependent)
end

% unpack solution vector
[~,x0] = bvp2sig(U,T,M);       
m_in = feval(sys.e,[],[],[],orb.sig(sl_ind),7,1);     % mode before event
m_out = feval(sys.e,[],[],[],orb.sig(sl_ind),7,2);    % mode after event
u_sl = x0(:,sl_ind*M);          % state at sliding event

% Find delayed terms
us_tau = del.ud; % interpolations of the delayed terms
i_tau = del.id; % containing segment indices
fi_tau = del.fi; % Lagrange interpolation coefficients
ud_sl = squeeze(us_tau(:,:,sl_ind*M));

% Evaluate relevant vector fields and event conditions
f_in = feval(sys.f,u_sl,ud_sl,p,m_in,1,0);  % vector field before event
f_out = feval(sys.f,u_sl,ud_sl,p,m_out,1,0); % vector field after event
dh = feval(sys.e,u_sl,ud_sl,p,orb.sig(sl_ind),2,0); % event surface jacobian
Jf_in = feval(sys.f,u_sl,ud_sl,p,m_in,3,0); % vector field Jacobian before event
Jf_out = feval(sys.f,u_sl,ud_sl,p,m_out,3,0); % vector field Jacobian after event

% Jacobian of sliding condition
h_sl = (dh*f_in + dh*f_out)/(dh*f_in - dh*f_out);
dh_sl = 1/(dh*f_in - dh*f_out)*(dh*(Jf_in+Jf_out))...
    -(dh*f_in + dh*f_out)/(dh*f_in - dh*f_out)^2*(dh*(Jf_in-Jf_out));

% Effect on time delay
dht_sl = zeros(1,length(p));
for k = 1:nt
    % Interpolation
    fi_k = squeeze(fi_tau(sl_ind*M,k,:)).'; % interpolation coefficients
    dfi_k = (D0.'*fi_k.').'; % time derivative of interpolation coefficients
    i_k = (i_tau(sl_ind*M,k)-1)*M*n+1:i_tau(sl_ind*M,k)*M*n; % indicies of interpolation segment elements
    dtau_k = -1/T(i_tau(sl_ind*M,k))*dp_tau(u_sl,k);
    % Corresponding derivatives
    Jfd_in = Jf_mj(m_in,u_sl,ud_sl,k); % vector field Jacobian before event
    Jfd_out = Jf_mj(m_out,u_sl,ud_sl,k); % vector field Jacobian after event
    dhd_sl = 1/(dh*f_in - dh*f_out)*(dh*(Jfd_in+Jfd_out))...
    -(dh*f_in + dh*f_out)/(dh*f_in - dh*f_out)^2*...
    (dh*(Jfd_in-Jfd_out));
    % Place elements in Jacobian matrix
    dht_sl = dht_sl + kron(dhd_sl,dfi_k)*U(i_k)*dtau_k;
end

% Pick relevant elements
Jsl = sign(h_sl)*(dh_sl+dht_sl);

end

