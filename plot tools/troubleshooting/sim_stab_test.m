function sim = sim_stab_test(orb,sys,opts,prt)
%SIM_STAB_TEST Test the stability of a periodic orbit with forward time
%simulation
% Input:
%   orb: periodic orbit data structure
%    -> U: solution vectors (M*n)
%    -> T: segment lengths (N)
%    -> p: parameter vector (l)
%    -> sig: solution signature (n)
%   sys: names of the functions that define the 
%    -> f: vector field and its Jacobians
%    -> e: event function, map and corresponding Jacobians
%    -> tau: time delay and its parameter Jacobian
%    -> mode_no: number of distinct vector field modes
%    -> event_no: number of disctinc events
%   opts: solver and data processing options
%    -> t_end: end time of simulation (default 1000)
%    -> m0: starting mode of simulation (default 1)
%    -> h_act: indicies of active events (default all)
%    -> i_max: maximum iteration number, event counter (default 1000)
%    -> calc_delayed: if true also include the delayed terms in the output
%       structure res (default true)
%   prt: norm of applied perturbation (default 0)
% Output:
%   sim: DDE solver output structure for troubleshooting and initial
%   searches

% Initialize orbit parameters
U = orb.U;
p0 = orb.p;
Tj = orb.T;

% apply some perturbation
if nargin>3
    U = U + prt*norm(U)*(1-2*rand(size(U)));
end

% Forward time simulation starting from the periodic orbit
T = sum(Tj);
y0 = @(t) po_interp(ceil(-t/T)*T+t,U,Tj,orb.M);
if ~isfield(sys,'ds_delay') || ~sys.sd_delay
    [~,sim] = sim_ns_dde(y0,p0,sys,[],opts);
else
    [~,sim] = sim_ns_sd_dde(y0,p0,sys,[],opts);
end
    
end

