function Jp = mpbvp_Jp(U,T,p,orb,sys,del)
%MPBVP_JP Jacobian for discretized multiple-point boundary value problem 
%of non-smooth periodic orbits in DDEs with respect to parameter vector
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
% Output:
%   Jp: Jacobian of governing system of nonlinear equations

% Initialization
M = orb.M;              % mesh resolution
n = orb.n;              % system dimension
N = length(orb.sig);    % number of segments       
l = length(p);          % number of system parameters
nt = sys.tau_no;        % number of time delays
D0 = cheb_diff(M);      % base differentiation matrix
Jfg = zeros(M*N*n,l);   % Jacobian matrix of rhs and interfaces
Jh = zeros(N,l);        % Jacobian of event conditions

% Function definitions
dtau = zeros(nt,l); % parameter derivatives of point time delays
for i=1:nt
    dtau(i,:) = feval(sys.tau,p,i,3);
end
pi_ej = @(ej) feval(sys.e,[],[],[],orb.sig(ej),7,1);    % incoming modes at events
Jfp_mj = @(mj,x,xd) feval(sys.f,x,xd,p,mj,3,0);         % vector field Jacobian in mode mj
Jhp_ej = @(j,x,xd) feval(sys.e,x,xd,p,orb.sig(j),3,0);  % event condition at ej
Jgp_ej = @(j,x,xd) feval(sys.e,x,xd,p,orb.sig(j),6,0);  % event map at ej
Jf_mj = @(mj,x,xd,i) feval(sys.f,x,xd,p,mj,2,i);         % vector field Jacobian in mode mj
Jh_ej = @(j,x,xd,i) feval(sys.e,x,xd,p,orb.sig(j),2,i);  % event condition at ej
Jg_ej = @(j,x,xd,i) feval(sys.e,x,xd,p,orb.sig(j),5,i);  % event map at ej

% Find current and delayed terms
[~,us] = bvp2sig(U,T,M); % signal form of state vector
us_tau = del.ud; % interpolations of the delayed terms
i_tau = del.id; % containing segment indices
fi_tau = del.fi; % Lagrange interpolation coefficients

% Go segment by segment (M*n rows each)
for j = 1:N
    % Indicies of corresponding solution segment
    ij = (j-1)*M*n+1:j*M*n; % in Jfg
    ii = (j-1)*M+1:j*M; % in us and ts
    
    % Current vector field mode
    mj = pi_ej(j);
    
    % Evalutation of DDE RHS and event map parameter Jacobians
    dfg = zeros(M,n*l);
    for k = 1:M
        ut = us(:,ii(k));
        utau = squeeze(us_tau(:,:,ii(k)));
        ij_k = (j-1)*M*n+k+M*(0:n-1); % row index of the current vector field or continuity condition
        
        if k<M % interior points
            dfg(k,:) = reshape(Jfp_mj(mj,ut,utau),1,[]); % ODE rhs
        else % boundary points
            dfg(k,:) = reshape(Jgp_ej(j,ut,utau),1,[]); % event map
            Jh(j,:) = Jh(j,:) +  Jhp_ej(j,ut,utau); % event condition
        end
        
        % Effect on time delay
         for i = 1:nt
            fi_jk = squeeze(fi_tau(ii(k),i,:)).'; % interpolation coefficients
            dfi_jk = (D0.'*fi_jk.').'; % Derivatives of Lagrange coefficients
            ijk = (i_tau(ii(k),i)-1)*M*n+1:i_tau(ii(k),i)*M*n; % indicies of interpolation segment elements
            dtau_jk = -1/T(i_tau(ii(k),i))*dtau(i,:);
            if k<M
                % Interior points
                df_jk = Jf_mj(mj,ut,utau,i); % Jacobian of f wrt u(t-tau(k))
                Jfg(ij_k,:) = Jfg(ij_k,:) - kron(df_jk,dfi_jk)*U(ijk)*dtau_jk;
            else
                % Boundary points
                dg_jk = Jg_ej(j,ut,utau,i); % Jacobian of g wrt u(t-tau(k))
                Jfg(ij_k,:) = Jfg(ij_k,:) - kron(dg_jk,dfi_jk)*U(ijk)*dtau_jk;
                dh_jk =  Jh_ej(j,ut,utau,i); % Jacobian of h wrt u(t-tau(k))
                Jh(j,:) = Jh(j,:) + kron(dh_jk,dfi_jk)*U(ijk)*dtau_jk;
            end            

         end
    end
    
    % Evaluate Jacobians
    Jfg(ij,:) = Jfg(ij,:) - reshape(dfg,[],l);
    
end

% Evaluate Jacobian matrix of the system of nonlinear equations
Jp = [Jfg; Jh];

end

