function Theta = mpbvp_mon(U,T,p,orb,sys,del)
%MPBVP_MON Formulate monodromy matrix of periodic orbits based on IFT
%for non-smooth autonomous DDEs
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
%    -> aut: if true the sytem is truly autonomous, thus the vector field 
%       condition for the final point of the orbit is included 
%       during mondoromy matrix formulation (default true)
%   del: data structure of delayed term evaluations
%    -> ud: state at t-tau(k) (n x nt x M*N*n) (ACCOUNTING FOR NEUTRAL DELAYS!)
%    -> id: segment index of t-tau(k) mapped back between 0 and T (M*n*N x nt)
%    -> fi: lagrange interpolation coefficients (M*n*N x nt x M)
%    -> dT: derivative of the querry point wrt Ti without looping around 
%       (M*n*N x nt x N x n_tau)
%    -> nd: index of the periodic orbit where the tq is found without 
%       looping around (0 current, 1,2,... past orbits) (M*n*N x nt)
% Output:
%   Theta: monodromy matrix of periodic orbit

if ~isfield(sys,'aut')
    sys.aut = true; % include the Jacobian of the f(x(T)) = x'(T) condition by default
end

% Initialization
M = orb.M;              % mesh resolution
n = orb.n;              % system dimension
N = length(orb.sig);    % number of segments
nt = sys.tau_no;        % number of time delays
D0 = cheb_diff(M);      % base differentiation matrix
jT = M*N*n+(1:N);       % indicies of T and h in the NAE

% Function definitions
pi_ej = @(ej) feval(sys.e,[],[],[],orb.sig(ej),7,1);     % incoming modes at events
Jf_mj = @(mj,x,xd,i) feval(sys.f,x,xd,p,mj,2,i);         % vector field Jacobian in mode mj
Jh_ej = @(j,x,xd,i) feval(sys.e,x,xd,p,orb.sig(j),2,i);  % event condition at ej
Jg_ej = @(j,x,xd,i) feval(sys.e,x,xd,p,orb.sig(j),5,i);  % event map at ej
td_type = @(i) feval(sys.tau,[],i,1);                    % delay type identifier

% Find current and delayed terms
[~,us] = bvp2sig(U,T,M); % signal form of state vector
us_tau = del.ud; % interpolations of the delayed terms
i_tau = del.id; % containing segment indices
fi_tau = del.fi; % Lagrange interpolation coefficients
dT = del.dT; % derivatives wrt segment lengths
nd_tau = del.nd; % index indicating the active periodic orbit instance

% Extend the Jacobian with the solution tail x0 on [-tau_max,0]
if nt>0
    ntau = max(nd_tau,[],'all');   % ntau*T is neccessary to cover tau_max
else
    ntau = 1; % calculation of stability for ODEs
end
Ju = repmat({zeros(n*M*N+N)},ntau+1,1); % Extended Jacobian matrix separated to nt+1 parts

% Go segment by segment (M*n+1 rows each)
for j = 1:N
    % Indicies of corresponding solution segment
    ij = (j-1)*M*n+1:j*M*n; % in D and fg
    ii = (j-1)*M+1:j*M; % in us and ts
    
    % Current vector field mode
    mj = pi_ej(j);
    
    % Spectral differentiation matrix
    DM = 1/T(j) * D0; % rescale with interval length
    if ~(sys.aut && j==N)
        DM(end,:) = 0; % make room for interface conditions 
    end
    D = kron(eye(n),DM);
    
    % Place interface condition counterparts in off diagonal blocks
    ut = us(:,ii(end));
    utau = squeeze(us_tau(:,:,ii(end)));
    dg_i =  Jg_ej(j,ut,utau,0); % event map jacobian shifted back by one instance
    ij_M = (j-1)*M*n+M*(1:n); % row index of the current continuity condition
    if j==N
        % keep only one side of the periodicity condition in Ju{end}
        g_ij = 1+M*(0:n-1); % map to the first segment
        Ju{end}(ij_M,g_ij) = Ju{end}(ij_M,g_ij) + eye(n);
        Ju{end-1}(ij_M,ij_M) = Ju{end-1}(ij_M,ij_M) - dg_i;
    else
        g_ij = j*M*n+1+M*(0:n-1); % map to the next segment
        Ju{end}(ij_M,g_ij) = Ju{end}(ij_M,g_ij) + eye(n);
        Ju{end}(ij_M,ij_M) = Ju{end}(ij_M,ij_M) - dg_i;
    end

    % Derivative wrt Ti
    Ju{end}(ij,jT(j)) = Ju{end}(ij,jT(j)) - 1/T(j)*D*U(ij); % [-1/Ti^2*DM]*u0
    
    % Evalutation of DDE RHS Jacobian
    df = zeros(M,n*n);
    for k = 1:M
        ut = us(:,ii(k));
        utau = squeeze(us_tau(:,:,ii(k)));
        ij_k = (j-1)*M*n+k+M*(0:n-1); % row index of the current vector field or continuity condition

        % Continuous terms
        if k<M || (sys.aut && j==N)
            % Inner points
            df(k,:) = reshape(Jf_mj(mj,ut,utau,0),1,[]); % DDE RHS
        end
        if k==M
            % Boundary points (event map g handled above)
            Ju{end}(jT(j),ij_k) = Ju{end}(jT(j),ij_k) + Jh_ej(j,ut,utau,0);  % event condition (move to Ju{end-1} if i==Ne?)
        end

        % Delayed terms wrt u(t-tau)
        for i = 1:nt

            fi_jk = squeeze(fi_tau(ii(k),i,:)).'; % Lagrange coefficients
            dfi_jk = (D0.'*fi_jk.').'; % Derivatives of Lagrange coefficients
            i_jk = i_tau(ii(k),i);  % index of segment containing t-tau_i
            ijk = (i_jk-1)*M*n+1:i_jk*M*n; % indicies of interpolation segment elements
            n_jk = nd_tau(ii(k),i); % index of orbit which contains t-tau

            if k<M || (sys.aut && j==N)
                % Interior points
                df_jk = Jf_mj(mj,ut,utau,i); % Jacobian of f wrt u(t-tau(k))
                Ju{end-n_jk}(ij_k,ijk) = Ju{end-n_jk}(ij_k,ijk) - kron(df_jk,fi_jk);
                for nk = 0:ntau
                    dT_njk = squeeze(dT(ii(k),i,:,end-nk)).';
                    Ju{end-nk}(ij_k,jT) = Ju{end-nk}(ij_k,jT)...
                        - kron(df_jk,dfi_jk)*U(ijk)*dT_njk;
                end
                if td_type(i) == 2 % extra term for neutral delays
                    Ju{end-n_jk}(ij_k,jT(i_jk)) = Ju{end-n_jk}(ij_k,jT(i_jk))...
                        + kron(df_jk,fi_jk)*U(ijk)/T(i_jk);
                end
            end

            if k==M
                % Boundary points
                dg_jk = Jg_ej(j,ut,utau,i); % Jacobian of g wrt u(t-tau(k))
                dh_jk = Jh_ej(j,ut,utau,i); % Jacobian of h wrt u(t-tau(k))
                if j==N % Moved back by 1 along with g_N
                    Ju{end-n_jk-1}(ij_k,ijk) = Ju{end-n_jk-1}(ij_k,ijk) ...
                        - kron(dg_jk,fi_jk);
                    for nk = 0:ntau-1
                        dT_njk = squeeze(dT(ii(k),i,:,end-nk)).';
                        Ju{end-nk-1}(ij_k,jT) = Ju{end-nk-1}(ij_k,jT)...
                            - kron(dg_jk,dfi_jk)*U(ijk)*dT_njk;
                    end
                    if td_type(i) == 2 % extra terms for neutral delays
                        Ju{end-n_jk-1}(ij_k,jT(i_jk)) = Ju{end-n_jk-1}(ij_k,jT(i_jk))...
                            + kron(dg_jk,fi_jk)*U(ijk)/T(i_jk);
                    end
                else
                    Ju{end-n_jk}(ij_k,ijk) = Ju{end-n_jk}(ij_k,ijk) ...
                        - kron(dg_jk,fi_jk);
                    for nk = 0:ntau
                        dT_njk = squeeze(dT(ii(k),i,:,end-nk)).';
                        Ju{end-nk}(ij_k,jT) = Ju{end-nk}(ij_k,jT)...
                            - kron(dg_jk,dfi_jk)*U(ijk)*dT_njk;
                    end
                    if td_type(i) == 2 % extra terms for neutral delays
                        Ju{end-n_jk}(ij_k,jT(i_jk)) = Ju{end-n_jk}(ij_k,jT(i_jk))...
                            + kron(dg_jk,fi_jk)*U(ijk)/T(i_jk);
                    end
                end
                Ju{end-n_jk}(jT(j),ijk) = Ju{end-n_jk}(jT(j),ijk) ...
                    + kron(dh_jk,fi_jk);
                for nk = 0:ntau
                    dT_njk = squeeze(dT(ii(k),i,:,end-nk)).';
                    Ju{end-nk}(jT(j),jT) = Ju{end-nk}(jT(j),jT)...
                        + kron(dh_jk,dfi_jk)*U(ijk)*dT_njk;
                end
                if td_type(i) == 2 % extra terms for neutral delays
                    Ju{end-n_jk}(jT(j),jT(i_jk)) = Ju{end-n_jk}(jT(j),jT(i_jk))...
                        - kron(dh_jk,fi_jk)*U(ijk)/T(i_jk);
                end
            end            
        end
    end
    
    % Evaluate Jacobians
    Ju{end}(ij,ij) = Ju{end}(ij,ij) + D - block_diag(df); % df in block diagonal form
    
end

% Concatenate Jacobian matricies (with extra continuity conditions)
Znm = zeros((ntau-1)*(n*M*N+N),n*M*N+N);
Enm = eye((ntau-1)*(n*M*N+N));
Ju_p = [Znm.' Ju{end}; -Enm Znm];
Ju_m = [Ju{1:end-1}; Znm Enm];

% Calculate monodromy matrix
Theta = -Ju_p\Ju_m;

end

