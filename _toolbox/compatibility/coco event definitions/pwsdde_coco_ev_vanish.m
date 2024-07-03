function prob = pwsdde_coco_ev_vanish(prob,T_lim,rel,boundary)
%PWSDDE_COCO_EV_VANISH Add monitor function and detection for vanishing
%segments to the COCO compatible problem
% Input:
%   prob: COCO compatible problem structure generated by pwsdde_coco_prob.m
%    -> init_orb: starting periodic orbit data structure
%   T_lim: minimum allowed segment length
%   rel: if true, the vanishing segment event is defined relative to
%       the period length (default false)
%   boundary: if true stop the continuation run when one of the  smooth 
%       segments vanishes (default true)
% Output:
%   prob: COCO compatible continuation problem structure with vanishing
%       segment detection function

if nargin<3 || isempty(rel)
    rel = false; % absolute segment length limit by default
end
if nargin<4
    boundary = true; % boundary event by default
end

% neccessary indicies
M = prob.init_orb.M;             % Chebyshev mesh resolution
n = prob.init_orb.n;             % number of degrees of freedom
N = length(prob.init_orb.sig);   % number of smooth segments
T_idx = M*n*N+(1:N);             % indicies of segment lengths

% add continuation variable corresponding to min(T)
fcn = @(f) @(prob,data,u) deal(data, f(u)); 
if rel
    f_min = @(T) min(T)/sum(T);
else
    f_min = @(T) min(T);
end
prob = coco_add_func(prob,'vanish',fcn(f_min),[],'regular',...
    'T_min','uidx',T_idx);

% add event corresponding to vanishing segments
if boundary
    prob = coco_add_event(prob,'VA','boundary', 'T_min','<',T_lim);
else
    prob = coco_add_event(prob,'VA','special point', 'T_min','=',T_lim);
end

end

