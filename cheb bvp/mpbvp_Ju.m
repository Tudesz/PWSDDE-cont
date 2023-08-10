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
% Output:
%   Ju: Jacobian of governing system of nonlinear equations wrt U

% Initialization
M = orb.M;              % mesh resolution
n = orb.n;              % system dimension
N = length(orb.sig);    % number of segments
D0 = cheb_diff(M);      % base differentiation matrix
Jfg = zeros(M*N*n);     % Jacobian matrix of rhs and interfaces
Jh = zeros(N,M*N*n);    % Jacobian of event conditions
JTfg = zeros(M*N*n,N);  % Jacobian of f and g with respect to T
JTh = zeros(N,N);       % Jacobian of h with respect to T

% Function definitions
tau = feval(sys.tau,p,1); % point time delays
pi_ej = @(ej) feval(sys.e,[],[],[],orb.sig(ej),7,1);    % incoming modes at events
Jf_mj = @(mj,x,xd,i) feval(sys.f,x,xd,p,mj,2,i);        % vector field Jacobian in mode mj
Jh_ej = @(j,x,xd,i) feval(sys.e,x,xd,p,orb.sig(j),2,i); % event condition at ej
Jg_ej = @(j,x,xd,i) feval(sys.e,x,xd,p,orb.sig(j),5,i); % event map at ej


% Evaluate current and delayed terms
[~,us] = bvp2sig(U,T,M); % signal form of state vector
[us_tau,i_tau,fi_tau,~,dT] = po_delay_interp(U,T,M,tau); % interpolation of delayed terms

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
        for i = 1:length(tau)
            fi_jk = zeros(M);
            fi_jk(k,:) = squeeze(fi_tau(ii(k),i,:)); % interpolation coefficients
            dfi_jk = zeros(M);
            dfi_jk(k,:) = (D0.'*fi_jk(k,:).').'; % time derivative of interpolation coefficients
            ijk = (i_tau(ii(k),i)-1)*M*n+1:i_tau(ii(k),i)*M*n; % indicies of interpolation segment elements
            dT_jk = squeeze(sum(dT(ii(k),i,:,:),4)).'; % derivatives wrt Ti
            if k<M
                % Interior points
                df_jk = Jf_mj(mj,ut,utau,i); % Jacobian of f wrt u(t-tau(k))
                Jfg(ij,ijk) = Jfg(ij,ijk) - kron(df_jk,fi_jk);
                JTfg(ij,:) = JTfg(ij,:) - kron(df_jk,dfi_jk)*U(ijk)*dT_jk;
            else
                % Boundary points
                dg_jk = Jg_ej(j,ut,utau,i); % Jacobian of g wrt u(t-tau(k))
                Jfg(ij,ijk) = Jfg(ij,ijk) - kron(dg_jk,fi_jk);
                JTfg(ij,:) = JTfg(ij,:) - kron(dg_jk,dfi_jk)*U(ijk)*dT_jk;
                dh_jk = Jh_ej(j,ut,utau,i); % Jacobian of h wrt u(t-tau(k))
                Jh(j,ijk) = Jh(j,ijk) + kron(dh_jk,fi_jk(k,:));
                JTh(j,:) = JTh(j,:) + kron(dh_jk,dfi_jk(k,:))*U(ijk)*dT_jk;
            end            
        end
    end
    
    % Evaluate Jacobians
    Jfg(ij,ij) = Jfg(ij,ij) + D - block_diag(df); % df in block diagonal form

end

% Concatenate Jacobian matricies
Ju = [[Jfg; Jh] [JTfg; JTh]];

end

