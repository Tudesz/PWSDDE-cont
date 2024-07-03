function sl = mpbvp_sl(U,T,p,orb,sys,sl_ind)
%MPBVP_SL Sliding condition for a DDE piecewise-smooth periodic orbit
% IMPORTANT: in the current version, the following assumptions must hold:
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
%   sl: error of governing sliding condition (1)

% Initialization
M = orb.M;                 % mesh resolution
[~,x0] = bvp2sig(U,T,M);   % unpack solution vector
m_in = feval(sys.e,[],[],[],orb.sig(sl_ind),7,1);     % mode before event
m_out = feval(sys.e,[],[],[],orb.sig(sl_ind),7,2);    % mode after event
u_sl = x0(:,sl_ind*M);     % state at sliding event

% Interpolation of delayed terms
[x0_tau,~,~,~,~] = po_delay_interp(U,T,p,M,sys);
ud_sl = squeeze(x0_tau(:,:,sl_ind*M));

% Evaluate relevant vector fields and event conditions
f_in = feval(sys.f,u_sl,ud_sl,p,m_in,1,0);  % vector field before event
f_out = feval(sys.f,u_sl,ud_sl,p,m_out,1,0); % vector field after event
dh = feval(sys.e,u_sl,ud_sl,p,orb.sig(sl_ind),2,0); % event surface jacobian

% Check sliding condition
h_sl = (dh*f_in + dh*f_out)/(dh*f_in - dh*f_out);
sl = abs(h_sl)-1;
end

