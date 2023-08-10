function sim = sim_stab_test(res,sys,ind,opts,prt)
%SIM_STAB_TEST Test the stability of a periodic orbit with forward time
%simulation
% Input:
%   res: continuation run output or orbit data structure
%    -> u: solution vectors (M*n)
%    -> T: segment lengths (N)
%    -> p: parameter vector (l)
%    -> sig: solution signature (n)
%   sys: names of the functions that define the 
%    -> f: vector field and its Jacobians
%    -> e: event function, map and corresponding Jacobians
%    -> tau: time delay and its parameter Jacobian
%    -> mode_no: number of distinct vector field modes
%    -> event_no: number of disctinc events
%   prt: norm of applied perturbation (default 0)
% Output:
%   sim: DDE solver output structure for troubleshooting and initial
%   searches

% Initialize orbit parameters
if length(ind) ~=1
    ind = length(res);
end
U = res(ind).U;
p0 = res(ind).p;
Tj = res(ind).T;

% apply some perturbation
if nargin>4
    U = U + prt*norm(U)*(1-2*rand(size(U)));
end

% Forward time simulation starting from the periodic orbit
T = sum(Tj);
y0 = @(t) po_interp(ceil(-t/T)*T+t,U,Tj,res(ind).M);
[~,sim] = sim_ns_dde(y0,p0,sys,[],opts);

end

