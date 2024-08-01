function prob = pwsdde_coco_ev_graze_bd(prob,sys,ev_ind,boundary)
%PWSDDE_COCO_EV_GRAZE_BD Add monitor function and detection for boundary
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
%   ev_ind: indicies of events to check for grazing (default all)
%   boundary: if true stop the continuation run when a grazing event is
%       detected (default true)
% Output:
%   prob: COCO compatible continuation problem structure with boundary
%   grazing detection functions

% initialization
orb0 = prob.init_orb;
fcn = @(f) @(prob,data,u) deal(data, f(u));

if nargin<3 || isempty(ev_ind)
    ev_ind = 1:length(orb0.sig);   % check all events by default
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

% span all monitored event boundaries
for i = ev_ind
    fid = sprintf('graze_bd%i',i);
    pid = sprintf('gr.bd(%i)',i);
    % add a continuation variable corresponding to the sliding surface
    gr_n = @(u) mpbvp_gr(u(U_idx),u(T_idx),u(P_idx),orb0,sys,...
        po_delay_interp(u(U_idx),u(T_idx),u(P_idx),orb0.M,sys),i);
    Jgr_n = @(u) [mpbvp_gr_Ju(u(U_idx),u(T_idx),u(P_idx),orb0,sys,...
        po_delay_interp(u(U_idx),u(T_idx),u(P_idx),orb0.M,sys),i) ...
        mpbvp_gr_Jp(u(U_idx),u(T_idx),u(P_idx),orb0,sys,...
        po_delay_interp(u(U_idx),u(T_idx),u(P_idx),orb0.M,sys),i)];
    prob = coco_add_func(prob,fid,fcn(gr_n),fcn(Jgr_n),...
        [],'regular',pid,'uidx',[U_idx,T_idx,P_idx]);
    % add event corresponding to reaching the sliding boundary
    if boundary
        prob = coco_add_event(prob,'BGR','boundary', pid,'=',0);
    else
        prob = coco_add_event(prob,'BGR','special point', pid,'=',0);
    end
end

end





