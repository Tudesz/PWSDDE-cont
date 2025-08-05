function [sim, orb_p] = sim_stab_test(orb,sys,opts,prt,orb_eval)
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
%    -> t_end: end time of simulation (default sum(T))
%    -> m0: starting mode of simulation (default exit mode of e_N)
%    -> h_act: indicies of active events (default all)
%    -> i_max: maximum iteration number, event counter (default 1000)
%    -> calc_delayed: if true also include the delayed terms in the output
%       structure res (default true)
%   prt: norm of applied perturbation or a perturbation vector of size M*n
%       or M*n+N (default 0)
%   orb_eval: if true also return a final orbit structure in orb_p
%       (default false)
% Output:
%   sim: DDE solver output structure for troubleshooting and initial
%   searches
%   orb_p: new orbit at the end of the simulation with the same structure as
%   orb (empty by default)

if nargin < 5
    orb_eval = false; % don't evaluate orb_p by default
end

% default simulation options
if nargin < 3 || isempty(opts) || ~isfield(opts,'t_end')
    opts.t_end = sum(orb.T);
end
if ~isfield(opts,'m0')
    opts.m0 = feval(sys.e,[],[],[],orb.sig(end),7,2); % exit mode of the last event
end


% Initialize orbit parameters
U = orb.U;
p0 = orb.p;
Tj = orb.T;

% apply some perturbation
if nargin>3 && ~isempty(prt)
    if isscalar(prt)
        U = U + prt*norm(U)*randn(size(U));
    elseif length(prt) == length(U)
        U = U + prt;
    elseif length(prt) == length(U) + length(Tj)
        U = U + prt(1:length(U));
        Tj = Tj + prt(length(U)+1:end);
    end
end

% final orbit data
if orb_eval
    o_init.M = orb.M;
    o_init.N = length(orb.sig);
    o_init.N0 = o_init.N;
else
    o_init = [];
end

% Forward time simulation starting from the periodic orbit
T = sum(Tj);
y0 = @(t) po_interp(ceil(-t/T)*T+t,U,Tj,orb.M);
if ~isfield(sys,'sd_delay') || ~sys.sd_delay
    [orb_p,sim] = sim_ns_dde(y0,p0,sys,o_init,opts);
else
    [orb_p,sim] = sim_ns_sd_dde(y0,p0,sys,o_init,opts);
end
    
end

