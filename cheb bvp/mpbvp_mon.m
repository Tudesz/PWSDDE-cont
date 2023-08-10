function Theta = mpbvp_mon(U,T,p,orb,sys)
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
% Output:
%   Theta: monodromy matrix of periodic orbit

% Initialization
M = orb.M;              % mesh resolution
n = orb.n;              % system dimension
N = length(orb.sig);    % number of segments        
D0 = cheb_diff(M);      % base differentiation matrix
jT = M*N*n+(1:N);       % indicies of T and h in the NAE

% Function definitions
tau = feval(sys.tau,orb.p,1); % point time delays
pi_ej = @(ej) feval(sys.e,[],[],[],orb.sig(ej),7,1);     % incoming modes at events
Jf_mj = @(mj,x,xd,i) feval(sys.f,x,xd,p,mj,2,i);         % vector field Jacobian in mode mj
Jh_ej = @(j,x,xd,i) feval(sys.e,x,xd,p,orb.sig(j),2,i);  % event condition at ej
Jg_ej = @(j,x,xd,i) feval(sys.e,x,xd,p,orb.sig(j),5,i);  % event map at ej


% Extend the Jacobian with the solution tail x0 on [-tau_max,0]
if ~isempty(tau)
    nt = ceil(max(tau)/sum(T));   % nt*T is neccessary to cover tau_max
else
    nt = 1; % calculation of stability for ODEs
end
Ju = repmat({zeros(n*M*N+N)},nt+1,1); % Extended Jacobian matrix separated to nt+1 parts

% Evaluate current and delayed terms
[~,us] = bvp2sig(U,T,M); % signal form of state vector
[us_tau,i_tau,fi_tau,nd_tau,dT] = po_delay_interp(U,T,M,tau); % interpolation of delayed terms

% Go segment by segment (M*n+1 rows each)
for j = 1:N
    % Indicies of corresponding solution segment
    ij = (j-1)*M*n+1:j*M*n; % in D and fg
    ii = (j-1)*M+1:j*M; % in us and ts
    
    % Current vector field mode
    mj = pi_ej(j);
    
    % Spectral differentiation matrix
    DM = 1/T(j) * D0; % rescale with interval length
    if j<N % disregard periodicity condition
        DM(end,:) = 0; % make room for interface conditions 
    end
    D = kron(eye(n),DM);
    
    % Place interface condition counterparts in off diagonal blocks
    ut = us(:,ii(end));
    utau = squeeze(us_tau(:,:,ii(end)));
    g_i =  Jg_ej(j,ut,utau,0); % event map jacobian shifted back by one instance
    De_p = zeros(M); De_p(end,1) = 1;
    De_m = zeros(M); De_m(end,end) = 1;
    if j==N
        % keep only one side of the periodicity condition in Ju{end}
        Ju{end}(ij,1:M*n) = Ju{end}(ij,1:M*n) + kron(eye(n),De_p);
        Ju{end-1}(ij,ij) = Ju{end-1}(ij,ij) - kron(g_i,De_m);
    else
        Ju{end}(ij,ij+M*n) = Ju{end}(ij,ij+M*n) + kron(eye(n),De_p);
        Ju{end}(ij,ij) = Ju{end}(ij,ij) - kron(g_i,De_m);
    end

    % Derivative wrt Ti
    Ju{end}(ij,jT(j)) = Ju{end}(ij,jT(j)) - 1/T(j)*D*U(ij); % [-1/Ti^2*DM]*u0
    
    % Evalutation of DDE RHS Jacobian
    df = zeros(M,n*n);
    for k = 1:M
        ut = us(:,ii(k));
        utau = squeeze(us_tau(:,:,ii(k)));
        % Continuous terms
        if k<M || j==N
            % Inner points
            df(k,:) = reshape(Jf_mj(mj,ut,utau,0),1,[]); % DDE RHS
        end
        if k==M
            % Boundary points (event map g handled above)
            Ju{end}(jT(j),ij(M:M:end)) = Ju{end}(jT(j),ij(M:M:end)) + ...
                Jh_ej(j,ut,utau,0);  % event condition (move to Ju{end-1} if i==Ne?)
        end
        % Delayed terms wrt u(t-tau)
        for i = 1:length(tau)
            fi_jk = zeros(M);
            fi_jk(k,:) = squeeze(fi_tau(ii(k),i,:)); % interpolation coefficients
            dfi_jk = zeros(M);
            dfi_jk(k,:) = (D0.'*fi_jk(k,:).').'; % time derivative of interpolation coefficients
            ijk = (i_tau(ii(k),i)-1)*M*n+1:i_tau(ii(k),i)*M*n; % indicies of interpolation segment elements
            n_jk = nd_tau(ii(k),i); % index of orbit which contains t-tau
            if k<M || j==N
                % Interior points
                df_jk = Jf_mj(mj,ut,utau,i); % Jacobian of f wrt u(t-tau(k))
                Ju{end-n_jk}(ij,ijk) = Ju{end-n_jk}(ij,ijk) - kron(df_jk,fi_jk);
                for nk = 0:nt
                    dT_njk = squeeze(dT(ii(k),i,:,end-nk)).';
                    Ju{end-nk}(ij,jT) = Ju{end-nk}(ij,jT)...
                        - kron(df_jk,dfi_jk)*U(ijk)*dT_njk;
                end
            end
            if k==M
                % Boundary points
                dg_jk = Jg_ej(j,ut,utau,i); % Jacobian of g wrt u(t-tau(k))
                dh_jk = Jh_ej(j,ut,utau,i); % Jacobian of h wrt u(t-tau(k))
                if j==N % Moved back by 1 along with g_N
                    Ju{end-n_jk-1}(ij,ijk) = Ju{end-n_jk-1}(ij,ijk) ...
                        - kron(dg_jk,fi_jk);
                    for nk = 0:nt-1
                        dT_njk = squeeze(dT(ii(k),i,:,end-nk)).';
                        Ju{end-nk-1}(ij,jT) = Ju{end-nk-1}(ij,jT)...
                            - kron(dg_jk,dfi_jk)*U(ijk)*dT_njk;
                    end
                else
                    Ju{end-n_jk}(ij,ijk) = Ju{end-n_jk}(ij,ijk) ...
                        - kron(dg_jk,fi_jk);
                    for nk = 0:nt
                        dT_njk = squeeze(dT(ii(k),i,:,end-nk)).';
                        Ju{end-nk}(ij,jT) = Ju{end-nk}(ij,jT)...
                            - kron(dg_jk,dfi_jk)*U(ijk)*dT_njk;
                    end
                end
                Ju{end-n_jk}(jT(j),ijk) = Ju{end-n_jk}(jT(j),ijk) ...
                    + kron(dh_jk,fi_jk(k,:));
                for nk = 0:nt
                    dT_njk = squeeze(dT(ii(k),i,:,end-nk)).';
                    Ju{end-nk}(jT(j),jT) = Ju{end-nk}(jT(j),jT)...
                        + kron(dh_jk,dfi_jk(k,:))*U(ijk)*dT_njk;
                end
            end            
        end
    end
    
    % Evaluate Jacobians
    Ju{end}(ij,ij) = Ju{end}(ij,ij) + D - block_diag(df); % df in block diagonal form
    
end

% Concatenate Jacobian matricies (with extra continuity conditions)
Znm = zeros((nt-1)*(n*M*N+N),n*M*N+N);
Enm = eye((nt-1)*(n*M*N+N));
Ju_p = [Znm.' Ju{end}; -Enm Znm];
Ju_m = [Ju{1:nt-1}, Ju{nt}; Znm Enm];

% Calculate monodromy matrix
Theta = -Ju_p\Ju_m;

end

