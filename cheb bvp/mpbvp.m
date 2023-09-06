function err = mpbvp(U,T,p,orb,sys)
%MPBVP Discretized multiple-point boundary value problem of
%autonomous non-smooth periodic orbits in DDEs
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
%   err: error of governing system of nonlinear equations (M*N*n + N)

% Initialization
M = orb.M;              % mesh resolution
n = orb.n;              % system dimension
N = length(orb.sig);    % number of segments
fg = zeros(M*N*n,1);    % function evaluations
h = zeros(N,1);         % interface conditions
D0 = cheb_diff(M);      % base differentiation matrix
D = zeros(M*N*n);       % full differentiation matrix

% Function definitions
tau = feval(sys.tau,p,1); % point time delays
pi_ej = @(ej) feval(sys.e,[],[],[],orb.sig(ej),7,1);  % incoming modes at events
f_mj = @(mj,x,xd) feval(sys.f,x,xd,p,mj,1,0);         % vector field in mode mj
h_ej = @(j,x,xd) feval(sys.e,x,xd,p,orb.sig(j),1,0);  % event condition at ej
g_ej = @(j,x,xd) feval(sys.e,x,xd,p,orb.sig(j),4,0);  % jump map at ej


% Evaluate current and delayed terms
[~,us] = bvp2sig(U,T,M); % signal form of state vector
[us_tau,~,~,~,~] = po_delay_interp(U,T,M,tau); % interpolation of delayed terms

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
    D(ij,ij) = kron(eye(n),DM);
    
    % Place interface condition counterparts in off diagonal blocks
    De = zeros(M); De(end,1) = 1;
    if j==N
        D(ij,1:M*n) = D(ij,1:M*n) + kron(eye(n),De);
    else
        D(ij,ij+M*n) = D(ij,ij+M*n) + kron(eye(n),De);
    end
    
    % Evalutation of DDE RHS, event maps, and event conditions
    f = zeros(M,n);
    for k = 1:M
        ut = us(:,ii(k));
        utau = squeeze(us_tau(:,:,ii(k)));
        if k<M % inner points
            f(k,:) = f_mj(mj,ut,utau).'; % DDE RHS
        else % endpoints
            f(k,:) = g_ej(j,ut,utau).'; % event map 
            h(j) = h_ej(j,ut,utau); % event condition
        end
    end
    fg(ij) = reshape(f,[],1);
end

% Concatenate error vector
err = [D*U-fg; h];

end

