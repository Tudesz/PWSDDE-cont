function sys = struct2func(struct)
%STRUCT2FUNC Convert old struct system definition to new function form
% Input:
%   struct: old system definition data structure
%     -> f{m}: governing autonomous vector spaces R^n -> R^n (DDE rhs)
%     -> h{i}: switching or event surfaces R^n -> R
%     -> g{i}: event maps R^n->R^n
%     -> Ju{m}: vector space Jacobians R^n -> R^(n x n)(DDE rhs)
%     -> Jud{m,k}: vector space Jacobians wrt x(t-tau(k))
%     -> dh{i}: event surface Jacobians R^n -> R^(1 x n)
%     -> dg{i}: event map Jacobians R^n -> R^(n x n)
%     -> dhd{i,k}: event surface Jacobians wrt x(t-tau(k))
%     -> dgd{i,k}: event map Jacobians wrt x(t-tau(k))
%     -> Jp{m}: Jacobian of f R^n -> R^(n x l)
%     -> dhp{i}: Jacobian of h R^n -> R^(1 x l)
%     -> dgp{i}: Jacobian of g R^n -> R^(n x l)
%     -> tau: system time delays
%     -> tau_dp: derivative of time delays wrt system parameters
%     -> pm{i}: mode changes on event surfaces
% Output:
%   sys: functions compatible with the new continuation framework
%     -> f: vector field definitions
%     -> e: event functions and maps
%     -> tau: time delays
%     -> mode_no: number of vector field modes
%     -> event_no: number of possible events

if isfield(struct,'tau')
    sys.tau_no = size(struct.Jud,2);
else
    sys.tau_no = 0;
end
sys.mode_no = length(struct.f);
sys.event_no = length(struct.h);

sys.f = @(x,xd,p,mode,type,l) struct2func_f(struct,x,xd,p,mode,type,l);
sys.e = @(x,xd,p,id,type,l) struct2func_e(struct,x,xd,p,id,type,l);
if sys.tau_no > 0
    sys.tau = @(p,ind,type) struct2func_tau(struct,p,ind,type);
else
    sys.tau = @(p,ind,type) [];
end

end

% Auxiliary functions

function X = struct2func_f(struct,x,xd,p,mode,type,l)
%STRUCT2FUNC_F Convert vector field definitions from old struct data to new
%function definition form
% (rhs of the pws-dde and its Jacobians)
% Input:
%   struct: old system definition data structure
%       -> f{m}: governing autonomous vector spaces R^n -> R^n (DDE rhs)
%       -> Ju{m}: vector space Jacobians R^n -> R^(n x n)(DDE rhs)
%       -> Jud{m,k}: vector space Jacobians wrt x(t-tau(k))
%       -> Jp{m}: Jacobian of f R^n -> R^(n x l)
%   x: state vector at time t (n)
%   xd: delayed state vectors at t-tau_l (n x L)
%   p: system parameter vector
%   mode: mode index of the vector field mode
%   type: flag for which property is requested:
%       1, The function evaluation itself f(x,xd,p)
%       2, Jacobian wrt x or xd J_x(x,xd,p)
%       3, Parameter Jacobians J_p(x,xd,p)
%   l: index of the requested delayed Jacobian J_x
% Output:
%   X: appropriate evaluation of the vector field of one of its Jacobians

if isfield(struct,'Jud')
    switch type
        case 1 % Vector field
            X = struct.f{mode}(x,xd,p);
        case 2 % Jacobian matricies
            if l>0
                X = struct.Jud{mode,l}(x,xd,p);
            else
                X = struct.Ju{mode}(x,xd,p);
            end
        case 3 % Parameter Jacobian
            X = struct.Jp{mode}(x,xd,p);
    end
else
    switch type
        case 1 % Vector field
            X = struct.f{mode}(x,p);
        case 2 % Jacobian matricies
            X = struct.Ju{mode}(x,p);
        case 3 % Parameter Jacobian
            X = struct.Jp{mode}(x,p);
    end
end

end

function X = struct2func_e(struct,x,xd,p,id,type,l)
%F_STRUCT2FUNC Convert event definitions from old struct data to new
%function definition form
% (event function, event map, corresponding jacobians, mode changes)
% Input:
%   struct: old system definition data structure
%       -> h{i}: switching or event surfaces R^n -> R
%       -> g{i}: event maps R^n->R^n
%       -> dh{i}: event surface Jacobians R^n -> R^(1 x n)
%       -> dg{i}: event map Jacobians R^n -> R^(n x n)
%       -> dhd{i,k}: event surface Jacobians wrt x(t-tau(k))
%       -> dgd{i,k}: event map Jacobians wrt x(t-tau(k))
%       -> dhp{i}: Jacobian of h R^n -> R^(1 x l)
%       -> dgp{i}: Jacobian of g R^n -> R^(n x l)
%       -> pm{i}: mode changes on event surfaces
%   x: state vector at time t (n)
%   xd: delayed state vectors at t-tau_l (n x L)
%   p: system parameter vector
%   id: event identifier
%   type: flag for which property is requested:
%       1, evaluation of the event condition h(x,xd,p)
%       2, Jacobian of h wrt x or xd J_h(x,xd,p)
%       3, Parameter Jacobian of h J_hp(x,xd,p)
%       4, evaluation of the jump map g(x,xd,p)
%       5, Jacobian of g wrt x or xd J_g(x,xd,p)
%       6, Parameter Jacobian of g J_gp(x,xd,p)
%       7, Induced mode change pi_e
%   l: index of the requested delayed Jacobian J_h or J_g
% Output:
%   X: appropriate evaluation of the event functions of Jacobians

if isfield(struct,'dhd')
    switch type
        case 1 % Event condition
            X = struct.h{id}(x,xd,p);
        case 2 % Event condition Jacobians
            if l>0
                X = struct.dhd{id,l}(x,xd,p);
            else
                X = struct.dh{id}(x,xd,p);
            end
        case 3 % Parameter Jacobain of h
            X = struct.dhp{id}(x,xd,p);
        case 4 % Event map
            X = struct.g{id}(x,xd,p);
        case 5 % Event map Jacobians
            if l>0
                X = struct.dgd{id,l}(x,xd,p);
            else
                X = struct.dg{id}(x,xd,p);
            end
        case 6 % Parameter Jacobian of 6
            X = struct.dgp{id}(x,xd,p);
        case 7 % Mode change
            X = struct.pm{id}(l);
    end
else
    switch type
        case 1 % Event condition
            X = struct.h{id}(x,p);
        case 2 % Event condition Jacobians
            X = struct.dh{id}(x,p);
        case 3 % Parameter Jacobain of h
            X = struct.dhp{id}(x,p);
        case 4 % Event map
            X = struct.g{id}(x,p);
        case 5 % Event map Jacobians
            X = struct.dg{id}(x,p);
        case 6 % Parameter Jacobian of 6
            X = struct.dgp{id}(x,p);
        case 7 % Mode change
            X = struct.pm{id}(l);
    end
end

end

function X = struct2func_tau(struct,p,ind,type)
%STRUCT2FUNC_TAU Convert time delay definitions from old struct data to new
%function definition form
% Input:
%   struct: old system definition data structure
%       -> tau: system time delays
%       -> tau_dp: derivative of time delays wrt system parameters
%   p: system parameter vector
%   ind: index of requested time delay
%   type: flag for which property is requested:
%       1, type of the delayed term:
%           -> 1, fix point delay (normal)
%           -> 2, fix point delay (neutral)
%       2, value of the time delay tau(p)
%       3, parameter jacobian of the time delay Jp_tau(p)
% Output:
%   X: appropriate evaluation of the delay or its parameter Jacobain

% Select the requested time delay
switch type
    case 1 % Delay type
        X = 1;
    case 2 % Time delays
        tau = struct.tau(p);
        X = tau(ind);
    case 3 % Time delay derivatives
        tau_dp = struct.tau_dp(p);
        X = tau_dp(ind,:);
end

end


