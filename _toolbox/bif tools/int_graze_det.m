function  gr_cond = int_graze_det(U,T,p,orb,sys)
%INT_GRAZE Interior grazing condition for event detection in periodic
%orbits of PWS-DDEs
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
%    -> event_no: number of disctinc events
% Output:
%   gr_cond: interior grazing condition for each event type (N*(M-2)*n) x event_no)

% Initialization
M = orb.M;              % mesh resolution
N = length(orb.sig);    % number of segments
m = sys.event_no;       % number of contact surfaces
h = zeros(M*N,m);       % event condition evaluations

% Function definitions
h_ei = @(i,x,xd) feval(sys.e,x,xd,p,i,1,0);  % event condition at ej

% Evaluate current and delayed terms
[~,us] = bvp2sig(U,T,M); % signal form of state vector
[us_tau,~,~,~,~] = po_delay_interp(U,T,p,M,sys); % interpolation of delayed terms

% Evaluate event conditions for all points
for k=1:length(us)
    for i=1:m
        h(k,i) = h_ei(i,us(:,k),squeeze(us_tau(:,:,k)));
    end
end

% Differentiation matrix for evaluating the Lie derivative
% D = zeros(M*N);
% DM = cheb_diff(M);
% for j = 1:N
%      ii = (j-1)*M+1:j*M; % indicies in D
%      D(ii,ii) = 1/T(j)*DM;
% end
% dh = D*h; %Lie derivative of the event condtions

% Formulate the interior grazing condition
% gr_cond = h.*(D*h); % with Lie derivative
gr_cond = h;    % without Lie derivative
gr_cond([1:M:N*M M:M:N*M],:) = []; % omit boundary points


end

