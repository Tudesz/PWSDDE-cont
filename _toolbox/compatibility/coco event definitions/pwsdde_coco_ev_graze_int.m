function prob = pwsdde_coco_ev_graze_int(prob,sys,ev_types,boundary)
%PWSDDE_COCO_EV_GRAZE_INT Add monitor function and detection for interior
%grazing bifurcations to the COCO compatible problem
% Input:
%   prob: COCO compatible problem structure generated by pwsdde_coco_prob.m
%    -> init_orb: starting periodic orbit data structure
%   sys: names of the functions that define the system
%    -> f: vector field and its Jacobians
%    -> e: event function, map and corresponding Jacobians
%    -> tau: time delay and its parameter Jacobian
%    -> mode_no: number of distinct vector field modes
%    -> event_no: number of disctinc events
%    -> tau_no: number of distinct time delays
%   ev_types: event types to check for sliding (default all)
%   boundary: if true stop the continuation run when a grazing event is
%       detected (default true)
% Output:
%   prob: COCO compatible continuation problem structure with interior
%   grazing detection functions

% initialization
orb0 = prob.init_orb;
fcn = @(f) @(prob,data,u) deal(data, f(u));

if nargin<3 || isempty(ev_types)
    ev_types = 1:sys.event_no;   % check all events by default
end
if nargin<4
    boundary = true; % boundary event by default
end

% neccessary indicies
M = orb0.M;             % Chebyshev mesh resolution
n = orb0.n;             % number of degrees of freedom
N = length(orb0.sig);   % number of smooth segments
L = length(orb0.p);     % number of system parameters
U_idx = 1:M*n*N;        % indicies of the state variable vector
T_idx = M*n*N+(1:N);    % indicies of segment lengths
P_idx = M*n*N+N+(1:L);  % indicies of system parameters

% add a continuation variable corresponding to each monitored event funciton
gr_int = @(u) gr_crit(u(U_idx),u(T_idx),u(P_idx),orb0,sys);
prob = coco_add_func(prob,'graze_int',fcn(gr_int),[],'regular',...
        coco_get_def_par_names('gr.in',ev_types),'uidx',[U_idx,T_idx,P_idx]);

% span all monitored event boundaries
for i = ev_types
    % add event corresponding to reaching the sliding boundary
    pid = sprintf('gr.in(%i)',i);
    if boundary
        prob = coco_add_event(prob,'IGR','boundary', pid,'=',0);
    else
        prob = coco_add_event(prob,'IGR','special point', pid,'=',0);
    end

end

end


% Auxiliary functions

function grc = gr_crit(U,T,p,orb,sys)
    % Critical values of the grazing conditions
    gr_cond = int_graze_det(U,T,p,orb,sys);
    [~,ind] = min(abs(gr_cond),[],'linear');
    grc = gr_cond(ind);
end




