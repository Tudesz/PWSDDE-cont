function [orb,res] = sim_ns_sd_dde(y0,p,sys,o_init,opts)
%SIM_NS_SD_DDE Simulation routine for non-smooth DDEs/NDDEs to find 
% periodic solutions while considering state dependent delays
% Input:
%   y0: intial histoy function @(t) -> R^n
%   p: system parameter vector
%   sys: names of the functions that define the 
%    -> f: vector field and its Jacobians
%    -> e: event function, map and corresponding Jacobians
%    -> tau: time delay and its parameter Jacobian
%    -> mode_no: number of distinct vector field modes
%    -> event_no: number of disctinc events
%    -> tau_no: number of distinct time delays
%    -> sd_delay: flag for state dependent delay definition
%       (should be true, if false use sim_ns_dde.m instead)
%   o_init: supplementary data for periodic orbit search (if empty only res
%   is returned)
%    -> N: number of events in the periodic orbit
%    -> M: Chebyshev mesh resolution
%    -> N0: index of the zeroth event (default 2*N)
%   opts: solver and data processing options
%    -> t_end: end time of simulation (default 1000)
%    -> m0: starting mode of simulation (default 1)
%    -> h_act: indicies of active events (default all)
%    -> i_max: maximum iteration number, event counter (default 1000)
%    -> calc_delayed: if true also include data on the delayed terms in 
%       the res output structure: feild Z (default true)
% Output:
%   orb: periodic orbit data structure (based on the simulation guess)
%    -> sig: sloution signature (event list)
%    -> U: state variable vector (M*N*n)
%    -> T: segment lengths (N)
%    -> p: parameter vector
%    -> n: number of degrees of freedom
%    -> M: Chebyshev mesh resolution
%   res: DDE solver output structure for troubleshooting and visualization
%    -> mode: active mode of the vector field
%    -> event: executed event at the end of the solution segment
%    -> Z: delayed instances of the solution y(t-tau) at sol.x

% Initialization
if nargin<5
    % if none of the solver options are not defined
    opts.i_max = 1000; % default maximum iteration number
    opts.h_act = true(1,sys.event_no); % by default all events are active
    opts.m0 = 1; % by default start with mode 1
    opts.t_end = 1000; % default simulation time
    opts.calc_delayed = true; % also return the delayed terms in res by default
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
    if ~isfield(opts,'calc_delayed')
        opts.calc_delayed = true;
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

% Sort out time delays
tau_type = zeros(1,sys.tau_no);
for i = 1:sys.tau_no
    tau_type(i) = feval(sys.tau,[],p,i,1); % delay types
end
tau_ind = find(tau_type==1); % index of normal delays
ntau_ind = find(tau_type==2);% index of nbeutral delays

% Run simulation loop
if isempty(ntau_ind)
    res = repmat(struct('solver',[],'history',[],...
        'x',[],'y',[],'xe',[],'ye',[],'ie',[],'stats',[],'yp',[],...
        'mode',[],'event',[],'Z',[],'tau',[]),1,opts.i_max);
else
    res = repmat(struct('solver',[],'history',[],'x',[],'y',[],...
        'xe',[],'ye',[],'ie',[],'stats',[],'yp',[],'IVP',[],...
        'mode',[],'event',[],'Z',[],'tau',[]),1,opts.i_max);
end
fprintf('\nRun simulation routine for non-smooth DDEs\n')
while ti(iter)<opts.t_end
    % Run simulation and extract data
    if isempty(ntau_ind)
        % SD DDE case
        ddefun = @(t,y,Z) feval(sys.f,y,Z,p,mi(iter),1,0);
        dely = @(t,y) delay_calc_sd(t,y,p,sys); % normal delays
        solver_opts = ddeset('Event',@(t,y,Z) ...
        event_func(y,Z,sys,p,mi(iter),opts.h_act));
        sol = ddesd(ddefun,dely,y_hist,[ti(iter) opts.t_end],...
            solver_opts);
    else
        % SD NDDE case
        nddefun = @(t,y,ydel,ypdel) feval(sys.f,y,...
            yd_sort(y,ydel,ypdel,sys,p),p,mi(iter),1,0);
        dely = @(t,y) delay_calc_sd(t,y,p,sys,tau_ind);   % normal delays
        delyp = @(t,y) delay_calc_sd(t,y,p,sys,ntau_ind); % neutral delays
        solver_opts = ddeset('Event',@(t,y,ydel,ypdel) ...
            event_func(y,yd_sort(y,ydel,ypdel,sys,p),...
            sys,p,mi(iter),opts.h_act));
        sol = ddensd(nddefun,dely,delyp,y_hist,...
            [ti(iter) opts.t_end],solver_opts);
    end
    sol.mode = mi(iter);
    sol.Z = [];
    sol.tau = [];
    sol.event = [];
    res(iter) = sol;
    ti(iter+1) = sol.x(end);

    % Include sorted delayed terms in res if needed
    if opts.calc_delayed && sys.tau_no>0
        Z = zeros(size(sol.y,1),sys.tau_no,size(sol.y,2));
        tau_sol = zeros(sys.tau_no,size(Z,3));
        for k = 1:sys.tau_no
            for i = 1:size(Z,3)
                tau_sol(k,i) = feval(sys.tau,sol.y(:,i),p,k,2);
                td = sol.x(i)-tau_sol(k,i);
                if tau_type(k) == 1
                    Z(:,k,i) = history_func(td,sol.y(:,end),y0,res,ti);
                else
                    Z(:,k,i) = yp_history(td,y0,res,ti);
                end
            end
        end
        res(iter).tau = tau_sol;
        res(iter).Z = Z;
    end

    
    % Find which event occured
    if isempty(sol.ie) || sol.x(end) == opts.t_end
        break
    else
        h_ind = 1:length(opts.h_act);
        ei_temp = intersect(sol.ie,h_ind(1==opts.h_act),'stable');
        ei(iter+1) = ei_temp(end);
        res(iter).event = ei(iter+1);
    end
        
    % Update solution
    if length(sol.x)>2
        % Find delayed terms
        y_tau = zeros(length(sol.y(:,end)),sys.tau_no);
        for i = 1:sys.tau_no
            td = sol.x(end) - feval(sys.tau,sol.y(:,end),p,i,2); % sd delay
            if tau_type(i) == 1 % normal delay
                y_tau(:,i) = history_func(td,sol.y(:,end),y0,res,ti);
            else % neutral delay
                y_tau(:,i) = yp_history(td,y0,res,ti);
            end

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

% Omit unused steps
res = res(1:iter);

if ~isempty(o_init)
    % Starting event index
    if ~isfield(o_init,'N0')
        N0 = 2*o_init.N;
    else
        N0 = o_init.N0;
    end
    % Find periodic solution
    i0 = iter-N0; % index of zeroth event
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

% Auxiliary event functions
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

% Auxiliary history function
function y = history_func(t,y1,y0,sols,tis)
    % Check which segment contains y(t)
    if t<=tis(1)
        % use original history function
        y = y0(t);
    elseif t>=max(tis)
        % return mapped new value for t=ti
         y = y1; 
    else
        % find index of segment containing t
        i0 = find(t<=tis,1); 
        y = deval(sols(i0-1),t);
    end
end

% Additional functions for NNDEs and ddensd
function Z = yd_sort(y,ydel,ypdel,sys,p)
    % Sort delayed terms for sys.f
    Z = zeros(length(y),sys.tau_no);
    id = 1; ipd = 1; % counters for normal and neutral delays
    for i = 1:sys.tau_no
        type = feval(sys.tau,[],p,i,1);
        if type == 1
            Z(:,i) = ydel(:,id);
            id = id + 1;
        else
            Z(:,i) = ypdel(:,ipd);
            ipd = ipd + 1;
        end
    end
end
function td = delay_calc_sd(t,y,p,sys,ind)
    % find t-tau_i(x(t))
    if nargin<5
        ind = 1:sys.tau_no;
    end
    % evaluate all requested delays
    tau = zeros(1,length(ind));
    for i = 1:length(ind)
        tau(i) = feval(sys.tau,y,p,ind(i),2);
    end
    td = t - tau;
end
% Auxiliary history function for neutral terms
function yp = yp_history(t,y0,sols,tis)
    % Check which segment contains y(t)
    if t<=tis(1)
        % use original history function
        yp = zeros(size(y0));
    else
        % find index of segment containing t
        i0 = find(t<=tis,1); 
        [~,yp] = deval(sols(i0-1),t);
    end
end

