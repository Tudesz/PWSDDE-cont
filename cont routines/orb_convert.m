function orb_mod = orb_convert(orb,sys,guess,res)
%ORB_CONVERT Convert orbit state vector to a different signature by 
% manually picking new event locations
% Input:
%   orb: starting periodic orbit data structure
%    -> sig: sloution signature (event list)
%    -> U: state variable vector (M*N*n)
%    -> T: segment lengths (N)
%    -> p: parameter vector
%    -> n: number of degrees of freedom
%    -> M: Chebyshev mesh resolution
%   sys: names of the functions that define the system
%    -> f: vector field and its Jacobians
%    -> e: event function, map and corresponding Jacobians
%    -> tau: time delay and its parameter Jacobian
%    -> mode_no: number of distinct vector field modes
%    -> event_no: number of disctinc events
%   guess: number of events in new solution signature or predefined Ti1 (if
%       0, a single segment is considered, if multiple values are given they 
%       are employed as a first apporximate for the event times)
%   res: interpolation resolution for event selection (default 1000)
% Output:
%   orb_mod: data structure of the new periodic orbit
%    -> sig: sloution signature (event list)
%    -> U: state variable vector (M*N*n)
%    -> T: segment lengths (N)
%    -> p: parameter vector
%    -> n: number of degrees of freedom
%    -> M: Chebyshev mesh resolution
%   err: MP-BVP error at its last evaluation

% Initialization, unpack orbit structure
if nargin<7
    res = 1000; % default search resolution
end
n = orb.n; M = orb.M;
U = orb.U; Tj = orb.T; p = orb.p;
T = sum(Tj);

% Function definitions
tau = feval(sys.tau,p,1); % point time delays
h_ej = @(j,x,xd) feval(sys.e,x,xd,p,j,1);  % event condition at ej

% Create higher resolution orbit for event search
[ts,us] = bvp2sig([U; U],[Tj; Tj],M,ceil(res/(2*length(Tj)))); % signal form of state vector
% interpolation of delayed terms 
[us_tau,~,~,~] = po_delay_interp([U; U],[Tj;Tj],M,...
   tau,ceil(res/(2*length(Tj)))); 

% Evaluate event surfaces in all points
ev_tol = zeros(sys.event_no,length(ts));
for j = 1:sys.event_no
    for i = 1:length(ts)
        ut = us(:,i);
        utau = squeeze(us_tau(:,:,i));
        ev_tol(j,i) = h_ej(j,ut,utau);
    end
end

% Select events based on plots
if length(guess) == 1 && guess ~=0
    figure('Name','Event selection');
    plot(ts,ev_tol); xlabel('$t$'); ylabel('$h_i$'); hold on 
    yline(0,':k'); hold off
    ev_in = ginput(guess);
    ev_t = ev_in(:,1);
    disp(ev_t.');
elseif guess == 0
    ev_t = [0 T];
else
    ev_t = guess; % predefined segment lengths
end
N = length(ev_t)-1;  % number of new segments

% Interpolate solution onto new mesh
[ts_u,ind_u] = unique(ts);
us_u = us(:,ind_u);
x0 = zeros(n,M*N);
t0 = zeros(1,M*N);
sig1 = zeros(1,N);
for j = 1:N
    [~,sig1(j)] = min(abs(ev_tol(:,find(ts>=ev_t(j+1),1)))); % determine which event occured
    ii = (j-1)*M+1:j*M; % indicies of u and t
    t0(ii) = cheb_mesh(M,[ev_t(j) ev_t(j+1)]); % Chebyshev mesh on segments
    for i = 1:n
        x0(i,(j-1)*M+1:j*M) = interp1(ts_u,us_u(i,:),t0(ii)); % spline interploation
    end
end
[U1, T1] = sig2bvp(t0,x0,M);
orb_mod = orb;
orb_mod.U = U1;
orb_mod.T = T1;
orb_mod.sig = sig1;

% Close figure
ob = findobj('type','figure','name','Event selection');
close(ob);
end

