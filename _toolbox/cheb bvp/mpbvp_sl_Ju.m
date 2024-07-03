function Jsl = mpbvp_sl_Ju(U,T,p,orb,sys,sl_ind)
%MPBVP_SL_JU Sliding condition Jacobian for a DDE piecewise-smooth 
%periodic orbit
% IMPORTANT: in the current version, these assumptions must hold:
%   -> dh/dx is independet of x and xd
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
%   sl_ind: index of sliding event
% Output:
%   Jsl: Jacobian of governing sliding condition (N*M*n+N x 1)

% Initialization
M = orb.M;              % mesh resolution
n = orb.n;              % system dimension
N = length(orb.sig);    % number of segments
nt = sys.tau_no;        % number of time delays
D0 = cheb_diff(M);      % base differentiation matrix

% Function definitions
Jf_mj = @(mj,x,xd,i) feval(sys.f,x,xd,p,mj,2,i); % vector field Jacobian in mode mj
td_type = @(i) feval(sys.tau,[],i,1);    % delay type identifier

% Initialization
[~,x0] = bvp2sig(U,T,M);        % unpack solution vector
m_in = feval(sys.e,[],[],[],orb.sig(sl_ind),7,1);     % mode before event
m_out = feval(sys.e,[],[],[],orb.sig(sl_ind),7,2);    % mode after event
u_sl = x0(:,sl_ind*M);          % state at sliding event

% Interpolation of delayed terms
[x0_tau,i_tau,fi_tau,~,dT] = po_delay_interp(U,T,p,M,sys);
ud_sl = squeeze(x0_tau(:,:,sl_ind*M));

% Evaluate relevant vector fields and event conditions
f_in = feval(sys.f,u_sl,ud_sl,p,m_in,1,0);  % vector field before event
f_out = feval(sys.f,u_sl,ud_sl,p,m_out,1,0); % vector field after event
dh = feval(sys.e,u_sl,ud_sl,p,orb.sig(sl_ind),2,0); % event surface jacobian
Jf_in = Jf_mj(m_in,u_sl,ud_sl,0); % vector field Jacobian before event
Jf_out = Jf_mj(m_out,u_sl,ud_sl,0); % vector field Jacobian after event

% Jacobian of sliding condition
h_sl = (dh*f_in + dh*f_out)/(dh*f_in - dh*f_out);
dh_sl = 1/(dh*f_in - dh*f_out)*(dh*(Jf_in+Jf_out))...
    -(dh*f_in + dh*f_out)/(dh*f_in - dh*f_out)^2*(dh*(Jf_in-Jf_out));

% Place values at relevant positions
Jsl = zeros(1,n*N*M+N);
Jsl(1,(sl_ind-1)*M*n+M:M:sl_ind*M*n) = sign(h_sl)*dh_sl;

% Add derivatives from delayed terms
for k = 1:nt
    % Interpolation
    fi_k = squeeze(fi_tau(sl_ind*M,k,:)).'; % interpolation coefficients
    dfi_k = fi_k*D0; % time derivative of interpolation coefficients
    i_k = i_tau(sl_ind*M,k); % index of segment containing t-tau_k
    ik = (i_k-1)*M*n+1:i_k*M*n; % indicies of interpolation segment elements
    dT_k = squeeze(sum(dT(sl_ind*M,k,:,:),4)).';

    % Corresponding derivatives
    Jfd_in = Jf_mj(m_in,u_sl,ud_sl,k); % vector field Jacobian before event
    Jfd_out = Jf_mj(m_out,u_sl,ud_sl,k); % vector field Jacobian after event
    dhd_sl = 1/(dh*f_in - dh*f_out)*(dh*(Jfd_in+Jfd_out))...
    -(dh*f_in + dh*f_out)/(dh*f_in - dh*f_out)^2*...
    (dh*(Jfd_in-Jfd_out));

    % Place elements in Jacobian matrix
    Jsl(1,ik) = Jsl(1,ik) + sign(h_sl)*kron(dhd_sl,fi_k);
    Jsl(1,end-N+1:end) = Jsl(1,end-N+1:end) + ...
        sign(h_sl)*kron(dhd_sl,dfi_k)*U(ik)*dT_k;
    if td_type(k) == 2 % extra term for neutral delays
        Jsl(1,end-N+i_k) = Jsl(1,end-N+i_k) - ...
            sign(h_sl)*kron(dhd_sl,fi_k)*U(ik)/T(i_k);
    end
end

end

