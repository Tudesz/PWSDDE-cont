function prob = pwsdde_coco_ev_slide(prob,sys,ev_ind,boundary)
%PWSDDE_COCO_EV_SLIDE Add monitor function and detection for sliding
%bifurcations to the COCO compatible problem
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
%   ev_ind: indicies of events to check for sliding (default all)
%   boundary: if true stop the continuation run when a sliding event is
%       detected (default true)
% Output:
%   prob: COCO compatible continuation problem structure with sliding
%   detection functions

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
    m_in = feval(sys.e,[],[],[],orb0.sig(i),7,1);
    m_out = feval(sys.e,[],[],[],orb0.sig(i),7,2);
    if m_in~=m_out % dont check events where no mode change is induced
        fid = sprintf('slide%i',i);
        pid = sprintf('sl(%i)',i);
        % add a continuation variable corresponding to the sliding surface
        sl_n = @(u) mpbvp_sl(u(U_idx),u(T_idx),u(P_idx),orb0,sys,i);
        Jsl_n = @(u) [mpbvp_sl_Ju(u(U_idx),u(T_idx),u(P_idx),orb0,sys,i) ...
            mpbvp_sl_Jp(u(U_idx),u(T_idx),u(P_idx),orb0,sys,i)];
        prob = coco_add_func(prob,fid,fcn(sl_n),fcn(Jsl_n),...
            [],'regular',pid,'uidx',[U_idx,T_idx,P_idx]);
        % add event corresponding to reaching the sliding boundary
        if boundary
            prob = coco_add_event(prob,'SL','boundary', pid,'<',0);
        else
            prob = coco_add_event(prob,'SL','special point', pid,'=',0);
        end
    end
end

end





