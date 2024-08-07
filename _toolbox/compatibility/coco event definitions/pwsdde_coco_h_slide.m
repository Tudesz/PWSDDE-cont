function prob = pwsdde_coco_h_slide(prob,sys,bif_ind)
%PWSDDE_COCO_H_SLIDE Add an extra zero function to the COCO compatible 
%problem for following sliding bifurcations
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
%   bif_ind: index of grazing bifurcation in the solution signature
% Output:
%   prob: COCO compatible continuation problem structure with extra zero
%       function for a grazing event

% initialization
orb0 = prob.init_orb;
fcn = @(f) @(prob,data,u) deal(data, f(u));

% neccessary indicies
M = orb0.M;             % Chebyshev mesh resolution
n = orb0.n;             % number of degrees of freedom
N = length(orb0.sig);   % number of smooth segments
L = length(orb0.p);     % number of system parameters
U_idx = 1:M*n*N;        % indicies of the state variable vector
T_idx = M*n*N+(1:N);    % indicies of segment lengths
P_idx = M*n*N+N+(1:L);  % indicies of system parameters

% Inclusion of the grazing condition in the MP-BVP
f = @(u) mpbvp_sl(u(U_idx),u(T_idx),u(P_idx),orb0,sys,...
    po_delay_interp(u(U_idx),u(T_idx),u(P_idx),orb0.M,sys),bif_ind);
Jf = @(u) [mpbvp_sl_Ju(u(U_idx),u(T_idx),u(P_idx),orb0,sys,...
    po_delay_interp(u(U_idx),u(T_idx),u(P_idx),orb0.M,sys),bif_ind) ...
    mpbvp_sl_Jp(u(U_idx),u(T_idx),u(P_idx),orb0,sys,...
    po_delay_interp(u(U_idx),u(T_idx),u(P_idx),orb0.M,sys),bif_ind)];
prob = coco_add_func(prob,'slide',fcn(f),fcn(Jf),[],'zero',...
    'uidx',[U_idx,T_idx,P_idx]);
end

