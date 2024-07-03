function prob = pwsdde_coco_prob(sys,orb0,p_names)
%PWSDDE_COCO_PROB Define COCO compatible continuation problem for PWS-DDEs
% Input:
%   sys: names of the functions that define the system
%    -> f: vector field and its Jacobians
%    -> e: event function, map and corresponding Jacobians
%    -> tau: time delay and its parameter Jacobian
%    -> mode_no: number of distinct vector field modes
%    -> event_no: number of disctinc events
%    -> tau_no: number of distinct time delays
%    -> aut: if true the sytem is truly autonomous, thus the vector field 
%       condition for the final point of the orbit is included 
%       during mondoromy matrix formulation (default true)
%   orb0: starting periodic orbit data structure
%    -> sig: solution signature (event list)
%    -> U: state variable vector (M*N*n)
%    -> T: segment lengths (N)
%    -> p: parameter vector
%    -> n: number of degrees of freedom
%    -> M: Chebyshev mesh resolution
%   p_names: list of parameter names in a cell array
% Output:
%   prob: COCO compatible continuation problem structure
%    -> init_orb: starting periodic orbit data structure

% necessary indicies
M = orb0.M;             % Chebyshev mesh resolution
n = orb0.n;             % number of degrees of freedom
N = length(orb0.sig);   % number of smooth segments
L = length(orb0.p);     % number of system parameters (must be equal to length(p_names))
U_idx = 1:M*n*N;        % indicies of the state variable vector
T_idx = M*n*N+(1:N);    % indicies of segment lengths
P_idx = M*n*N+N+(1:L);  % indicies of system parameters

% initialization
prob = coco_prob;
fcn = @(f) @(prob,data,u) deal(data, f(u));
prob.init_orb = orb0; % append initial orbit data structure

% Inclusion of the governing MP-BVP
f = @(u) mpbvp(u(U_idx),u(T_idx),u(P_idx),orb0,sys);
Jf = @(u) [mpbvp_Ju(u(U_idx),u(T_idx),u(P_idx),orb0,sys) ...
    mpbvp_Jp(u(U_idx),u(T_idx),u(P_idx),orb0,sys)];
u0 = [orb0.U; orb0.T; reshape(orb0.p,[],1)];
prob = coco_add_func(prob,'mp_bvp',fcn(f),fcn(Jf),[],'zero','u0',u0);

% Define continuation parameters
prob = coco_add_pars(prob,'pars',P_idx,p_names);

% Evaluate critical Floquet multiplier
mu_n = @(u) mu_crit(mpbvp_mon(u(U_idx),u(T_idx),u(P_idx),orb0,sys));
prob = coco_add_func(prob,'floq',fcn(mu_n),[],'regular',...
    {'mu_crit'},'uidx',[U_idx,T_idx,P_idx]);

% Add complementary monitor functions
u_n = @(u) norm(po_ampl(u(U_idx),u(T_idx),M)); % norm of solution amplitudes
prob = coco_add_func(prob,'U_norm',fcn(u_n),[],'regular',...
    {'Ampl'},'uidx',[U_idx,T_idx]);
e_n = @(u) norm(mpbvp(u(U_idx),u(T_idx),u(P_idx),orb0,sys)); % solution error
prob = coco_add_func(prob,'f_norm',fcn(e_n),[],'regular',...
    {'|f|'},'uidx',[U_idx,T_idx,P_idx]);

% Add detection routines for stability related bifurcations
% prob = coco_add_event(prob,'SN','special point', 'mu_crit','=',1);
% prob = coco_add_event(prob,'PD','special point', 'mu_crit','=',-1);

end


% Auxiliary functions

function muc = mu_crit(Th)
    % Evaluate critical Floquet-multiplier from monodormy matrix
    mu = eig(Th);
    [~,ind] = max(abs(mu));
    muc = mu(ind);
end

function ampl = po_ampl(U,T,M)
    % Evaluate solution amplitudes
    [~,u] = bvp2sig(U,T,M,100);
    ampl = (max(u,[],2)-min(u,[],2))./2;
end
