function [orb,res] = sim_ns_dde(y0,p,sys,o_init,opts)
%SIM_NS_DDE Simulation routine for non-smooth DDEs to find periodic
%solutions
% Input:
%   y0: intial histoy function @(t) -> R^n
%   par: system parameter vector
%   sys: names of the functions that define the 
%    -> f: vector field and its Jacobians
%    -> e: event function, map and corresponding Jacobians
%    -> tau: time delay and its parameter Jacobian
%    -> mode_no: number of distinct vector field modes
%    -> event_no: number of disctinc events
%   o_init: supplementary data for periodic orbit search (if empty only res
%   is returned)
%    -> N: number of events in the periodic orbit
%    -> M: Chebyshev mesh resolution
%   opts: solver and data processing options
%    -> t_end: end time of simulation (default 1000)
%    -> m0: starting mode of simulation (default 1)
%    -> h_act: indicies of active events (default all)
%    -> i_max: maximum iteration number, event counter (default 1000)
% Output:
%   orb: periodic orbit data structure (based on the simulation guess)
%    -> sig: sloution signature (event list)
%    -> U: state variable vector (M*N*n)
%    -> T: segment lengths (N)
%    -> p: parameter vector
%    -> n: number of degrees of freedom
%    -> M: Chebyshev mesh resolution
%   res: DDE solver output structure for troubleshooting and initial
%   searches

% Initialization
if nargin<5
    % if none of the solver options are not defined
    opts.i_max = 1000; % default maximum iteration number
    opts.h_act = true(1,sys.event_no); % by default all events are active
    opts.m0 = 1; % by default start with mode 1
    opts.t_end = 1000; % default simulation time
else
    if ~isfield(opts,'i_max')
        opts.i_max = 1000;
    end
    if ~isfield(opts,'h_act')
        opts.h_act = true(1,sys.event_no);
    end
    if ~isfield(opts,'m0')
        opts.m0 = 1;
    end
    if ~isfield(opts,'t_end')
        opts.t_end = 1000;
    end
end
iter = 1; % iteration counter (safety switch)
ti = zeros(opts.i_max,1); % starting times
mi = zeros(opts.i_max,1); % vector field modes
ei = zeros(opts.i_max,1); % encountered events
n = length(y0(0));  % system dimension
mi(1) = opts.m0;    % starting vector field mode
y_hist = @(t) history_func(t,y0(0),y0,[],ti);
rt = tic;
tau = feval(sys.tau,p,1);
if isempty(tau) 
    tau = 1; % introduce a dummy time delay
end

% Run simulation loop
fprintf('\nRun simulation routine for non-smooth DDEs\n')
while ti(iter)<opts.t_end
    % Update solver data
    ddefun = @(t,y,Z) feval(sys.f,y,Z,p,mi(iter),1,0);
    solver_opts = ddeset('Event',@(t,y,Z) ...
        event_func(y,Z,sys,p,mi(iter),opts.h_act));
    
    % Run simulation and extract data
    sol = dde23(ddefun,tau,y_hist,[ti(iter) opts.t_end],solver_opts);
    sol.mode = mi(iter);
    res(iter) = sol; % solution
    
    % Find which event occured
    if isempty(sol.ie) || sol.x(end) == opts.t_end
        break
    else
        h_ind = 1:length(opts.h_act);
        ei_temp = intersect(sol.ie,h_ind(1==opts.h_act),'stable');
        ei(iter+1) = ei_temp(end); 
    end
        
    % Update solution
    ti(iter+1) = sol.x(end);
    if length(sol.x)>2
        % Find delayed terms
        hist_temp = @(t) history_func(t,sol.y(:,end),y0,res,ti);
        y_tau = zeros(length(sol.y(:,end)),length(tau));
        for i = 1:length(tau)
            y_tau(:,i) = hist_temp(sol.x(end)-tau(i));
        end
        yt = sol.y(:,end);
        % Apply event map and mode change
        y1 = feval(sys.e,yt,y_tau,p,ei(iter+1),4,0);
        mi(iter+1) = feval(sys.e,yt,y_tau,p,ei(iter+1),7,2);
    else
        % Continue solution as if nothing had happened
        y1 = sol.y(:,end);
        mi(iter+1) = mi(iter);
    end
    y_hist = @(t) history_func(t,y1,y0,res,ti);
    
    % Safety switch 
    if iter+1 > opts.i_max
        warning('Maximum iteration %i reached',opts.i_max);
        break
    end
    iter = iter+1;
end

if ~isempty(o_init)
    % Find periodic solution
    i0 = iter-2*o_init.N; % index of zeroth event
    t = zeros(o_init.M*o_init.N,1); % time mesh
    u = zeros(n,o_init.M*o_init.N); % solution vectors
    orb.sig = zeros(1,o_init.N); % solution signature
    for j = 1:o_init.N
        jj = (j-1)*o_init.M+1:j*o_init.M;
        t(jj) = cheb_mesh(o_init.M,[ti(i0-1+j) ti(i0+j)]);
        u(:,jj) = deval(res(i0-1+j),t(jj));
        orb.sig(j) = ei(i0+j);
    end
    t = t-t(1); % shift time mesh back to 0
    orb.M = o_init.M; orb.p = p; orb.n = n;
    [orb.U, orb.T] = sig2bvp(t,u,orb.M);
else
    orb = [];
end

fprintf('   Simulation time: %0.3f s\n',toc(rt));

end

% Auxilarry event functions
function [value,isterminal,direction] = event_func(y,Z,sys,par,mode,h_act)
    % Checking event surface conditions
    value = ones(1,sys.event_no);
    for i = 1:sys.event_no
        if feval(sys.e,y,Z,par,i,7,1) == mode % only check currently active event surfaces
            value(i) = feval(sys.e,y,Z,par,i,1,0); % detect zero crossings
        end
    end
    isterminal = h_act; % terminate integration for allowed events
    direction = zeros(1,sys.event_no); % don't check crossing direction
end

% Auxilarry history function
function y = history_func(t,y1,y0,sols,tis)
    % Check which segment contains y(t)
    if t<=tis(1)
        % use original history function
        y = y0(t);
    elseif t==max(tis)
        % return mapped new value for t=ti
         y = y1; 
    else
        % find index of segment containing t
        i0 = find(t<=tis,1); 
        y = deval(sols(i0-1),t);
    end
end

