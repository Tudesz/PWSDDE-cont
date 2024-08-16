function Ju = mpbvp_Ju(U,T,p,orb,sys,del)
%MPBVP_JU Jacobian for the discretized multiple-point boundary value problem 
%of autonomous non-smooth periodic orbits in DDEs with respect to the 
%state vector U
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
%    -> dT: derivative of the querry point wrt Ti without looping around 
%       (M*n*N x nt x N x n_tau)
% Output:
%   Ju: Jacobian of governing system of nonlinear equations wrt U

if ~isfield(sys,'sd_delay')
    sys.sd_delay = false; % by default use fix point delays
end

% Initialization
M = orb.M;              % mesh resolution
n = orb.n;              % system dimension
N = length(orb.sig);    % number of segments
nt = sys.tau_no;        % number of time delays
D0 = cheb_diff(M);      % base differentiation matrix
Jfg = zeros(M*N*n);     % Jacobian matrix of rhs and interfaces
Jh = zeros(N,M*N*n);    % Jacobian of event conditions
JTfg = zeros(M*N*n,N);  % Jacobian of f and g with respect to T
JTh = zeros(N,N);       % Jacobian of h with respect to T

% Function definitions
pi_ej = @(ej) feval(sys.e,[],[],[],orb.sig(ej),7,1);    % incoming modes at events
Jf_mj = @(mj,x,xd,i) feval(sys.f,x,xd,p,mj,2,i);        % vector field Jacobian in mode mj
Jh_ej = @(j,x,xd,i) feval(sys.e,x,xd,p,orb.sig(j),2,i); % event condition at ej
Jg_ej = @(j,x,xd,i) feval(sys.e,x,xd,p,orb.sig(j),5,i); % event map at ej
if ~sys.sd_delay
    td_type = @(i) feval(sys.tau,[],i,1); % delay type identifier
else
    td_dx = @(x,i) feval(sys.tau,x,p,i,4);   % time delay jacobian wrt x
    td_type = @(i) feval(sys.tau,[],[],i,1); % delay type identifier
end

% Find current and delayed terms
[~,us] = bvp2sig(U,T,M); % signal form of state vector
us_tau = del.ud; % interpolations of the delayed terms
i_tau = del.id; % containing segment indices
fi_tau = del.fi; % Lagrange interpolation coefficients
dT = del.dT; % derivatives wrt segment lengths

% Go segment by segment (M*n+1 rows each)
for j = 1:N
    % Indicies of corresponding solution segment
    ij = (j-1)*M*n+1:j*M*n; % in D and fg
    ii = (j-1)*M+1:j*M; % in us and ts
    
    % Current vector field mode
    mj = pi_ej(j);
    
    % Spectral differentiation matrix
    DM = 1/T(j) * D0; % rescale with interval length
    DM(end,:) = 0; % make room for interface conditions
    D = kron(eye(n),DM);    
    
    % Place interface condition counterparts in off diagonal blocks
    ij_M = (j-1)*M*n+M*(1:n); % row index of the current continuity condition
    if j==N
        g_ij = 1+M*(0:n-1); % map to the first segment
    else
        g_ij = j*M*n+1+M*(0:n-1); % map to the next segment
    end
    Jfg(ij_M,g_ij) = Jfg(ij_M,g_ij) + eye(n);
    
    % Derivative wrt Ti
    JTfg(ij,j) = JTfg(ij,j) - 1/T(j)*D*U(ij); % [-1/Ti^2*DM]*u0
    
    % Evalutation of DDE RHS Jacobian
    df = zeros(M,n*n);
    for k = 1:M
        ut = us(:,ii(k));
        utau = squeeze(us_tau(:,:,ii(k)));
        ij_k = (j-1)*M*n+k+M*(0:n-1); % row index of the current vector field or continuity condition
        
        % Continuous terms
        if k<M
            % Inner points
            df(k,:) = reshape(Jf_mj(mj,ut,utau,0),1,[]); % DDE RHS
        else
            % Boundary points
            df(k,:) = reshape(Jg_ej(j,ut,utau,0),1,[]); % event map
            Jh(j,ij_k) = Jh(j,ij_k) + Jh_ej(j,ut,utau,0);  % event condition
        end

        % Delayed terms wrt u(t-tau)
        for i = 1:nt

            i_jk = i_tau(ii(k),i);  % index of segment containing t-tau_i
            ijk = (i_jk-1)*M*n+1:i_jk*M*n; % indicies of interpolation segment elements
            fi_jk = squeeze(fi_tau(ii(k),i,:)).'; % Lagrange coefficients
            dfi_jk = (D0.'*fi_jk.').'; % Derivatives of Lagrange coefficients
            dT_jk = squeeze(sum(dT(ii(k),i,:,:),4)).'; % derivatives of t-tau_i wrt Ti
            if sys.sd_delay
                dtau_dx = -1/T(i_jk)*td_dx(ut,i); % derivative of tau wrt x
            end
            if k<M
                % Interior points
                df_jk = Jf_mj(mj,ut,utau,i); % Jacobian of f wrt u(t-tau(k))
                Jfg(ij_k,ijk) = Jfg(ij_k,ijk) - kron(df_jk,fi_jk);
                JTfg(ij_k,:) = JTfg(ij_k,:) - kron(df_jk,dfi_jk)*U(ijk)*dT_jk;
                if td_type(i) == 2 % extra term for neutral delays
                    JTfg(ij_k,i_jk) = JTfg(ij_k,i_jk) + kron(df_jk,fi_jk)*U(ijk)/T(i_jk);
                end
                if sys.sd_delay % extra terms for state dependent delays
                    Jfg(ij_k,ij_k) = Jfg(ij_k,ij_k) - ...
                        kron(df_jk,dfi_jk)*U(ijk)*dtau_dx;
                end
            else
                % Boundary points
                dg_jk = Jg_ej(j,ut,utau,i); % Jacobian of g wrt u(t-tau(k))
                Jfg(ij_k,ijk) = Jfg(ij_k,ijk) - kron(dg_jk,fi_jk);
                JTfg(ij_k,:) = JTfg(ij_k,:) - kron(dg_jk,dfi_jk)*U(ijk)*dT_jk;
                if td_type(i) == 2 % extra term for neutral delays
                    JTfg(ij_k,i_jk) = JTfg(ij_k,i_jk) + kron(dg_jk,fi_jk)*U(ijk)/T(i_jk);
                end
                if sys.sd_delay % extra terms for state dependent delays
                    Jfg(ij_k,ij_k) = Jfg(ij_k,ij_k) - ...
                        kron(dg_jk,dfi_jk)*U(ijk)*dtau_dx;
                end
                dh_jk = Jh_ej(j,ut,utau,i); % Jacobian of h wrt u(t-tau(k))
                Jh(j,ijk) = Jh(j,ijk) + kron(dh_jk,fi_jk);
                JTh(j,:) = JTh(j,:) + kron(dh_jk,dfi_jk)*U(ijk)*dT_jk;
                if td_type(i) == 2 % extra term for neutral delays
                    JTh(j,i_jk) = JTh(j,i_jk) - kron(dh_jk,fi_jk)*U(ijk)/T(i_jk);
                end
                if sys.sd_delay % extra terms for state dependent delays
                    Jh(j,ij_k) = Jh(j,ij_k) + ...
                        kron(dh_jk,dfi_jk)*U(ijk)*dtau_dx;
                end
            end            
        end
    end
    
    % Evaluate Jacobians
    Jfg(ij,ij) = Jfg(ij,ij) + D - block_diag(df); % df in block diagonal form

end

% Concatenate Jacobian matricies
Ju = [[Jfg; Jh] [JTfg; JTh]];

end

