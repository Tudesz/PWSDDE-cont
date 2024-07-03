function Ju = mpbvp_Ju(U,T,p,orb,sys)
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
% Output:
%   Ju: Jacobian of governing system of nonlinear equations wrt U

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
td_type = @(i) feval(sys.tau,[],i,1);                   % delay type identifier

% Evaluate current and delayed terms
[~,us] = bvp2sig(U,T,M); % signal form of state vector
[us_tau,i_tau,fi_tau,~,dT] = po_delay_interp(U,T,p,M,sys); % interpolation of delayed terms

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
    De = zeros(M); De(end,1) = 1;
    if j==N
        Jfg(ij,1:M*n) = Jfg(ij,1:M*n) + kron(eye(n),De);
    else
        Jfg(ij,ij+M*n) = Jfg(ij,ij+M*n) + kron(eye(n),De);
    end
    
    % Derivative wrt Ti
    JTfg(ij,j) = JTfg(ij,j) - 1/T(j)*D*U(ij); % [-1/Ti^2*DM]*u0
    
    % Evalutation of DDE RHS Jacobian
    df = zeros(M,n*n);
    for k = 1:M
        ut = us(:,ii(k));
        utau = squeeze(us_tau(:,:,ii(k)));
        
        % Continuous terms
        if k<M
            % Inner points
            df(k,:) = reshape(Jf_mj(mj,ut,utau,0),1,[]); % DDE RHS
        else
            % Boundary points
            df(k,:) = reshape(Jg_ej(j,ut,utau,0),1,[]); % event map
            Jh(j,ij(M:M:end)) = Jh(j,ij(M:M:end)) + Jh_ej(j,ut,utau,0);  % event condition
        end

        % Delayed terms wrt u(t-tau)
        for i = 1:nt

            i_jk = i_tau(ii(k),i);  % index of segment containing t-tau_i
            ijk = (i_jk-1)*M*n+1:i_jk*M*n; % indicies of interpolation segment elements
            fi_jk = zeros(M);  % Lagrange coefficients
            dfi_jk = zeros(M); % Derivatives of Lagrange coefficients
            fi_jk(k,:) = squeeze(fi_tau(ii(k),i,:));
            dfi_jk(k,:) = (D0.'*fi_jk(k,:).').';
            dT_jk = squeeze(sum(dT(ii(k),i,:,:),4)).'; % derivatives of t-tau_i wrt Ti
            
            if k<M
                % Interior points
                df_jk = Jf_mj(mj,ut,utau,i); % Jacobian of f wrt u(t-tau(k))
                Jfg(ij,ijk) = Jfg(ij,ijk) - kron(df_jk,fi_jk);
                JTfg(ij,:) = JTfg(ij,:) - kron(df_jk,dfi_jk)*U(ijk)*dT_jk;
                if td_type(i) == 2 % extra term for neutral delays
                    JTfg(ij,i_jk) = JTfg(ij,i_jk) + kron(df_jk,fi_jk)*U(ijk)/T(i_jk);
                end
            else
                % Boundary points
                dg_jk = Jg_ej(j,ut,utau,i); % Jacobian of g wrt u(t-tau(k))
                Jfg(ij,ijk) = Jfg(ij,ijk) - kron(dg_jk,fi_jk);
                JTfg(ij,:) = JTfg(ij,:) - kron(dg_jk,dfi_jk)*U(ijk)*dT_jk;
                if td_type(i) == 2 % extra term for neutral delays
                    JTfg(ij,i_jk) = JTfg(ij,i_jk) + kron(dg_jk,fi_jk)*U(ijk)/T(i_jk);
                end
                dh_jk = Jh_ej(j,ut,utau,i); % Jacobian of h wrt u(t-tau(k))
                Jh(j,ijk) = Jh(j,ijk) + kron(dh_jk,fi_jk(k,:));
                JTh(j,:) = JTh(j,:) + kron(dh_jk,dfi_jk(k,:))*U(ijk)*dT_jk;
                if td_type(i) == 2 % extra term for neutral delays
                    JTh(j,i_jk) = JTh(j,i_jk) - kron(dh_jk,fi_jk(k,:))*U(ijk)/T(i_jk);
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

