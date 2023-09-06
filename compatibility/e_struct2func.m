function X = e_struct2func(sys,x,xd,p,id,type,l)
%F_STRUCT2FUNC Convert event definitions from old struct data to new
%function definition form
% (event function, event map, corresponding jacobians, mode changes)
% Input:
%   sys: old system definition data structure
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

switch type
    case 1 % Event condition
        X = sys.h{id}(x,xd,p);
    case 2 % Event condition Jacobians
        if l>0
            X = sys.dhd{id,l}(x,xd,p);
        else
            X = sys.dh{id}(x,xd,p);
        end
    case 3 % Parameter Jacobain of h
        X = sys.dhp{id}(x,xd,p);
    case 4 % Event map
        X = sys.g{id}(x,xd,p);
    case 5 % Event map Jacobians
        if l>0
            X = sys.dgd{id,l}(x,xd,p);
        else
            X = sys.dg{id}(x,xd,p);
        end
    case 6 % Parameter Jacobian of 6
        X = sys.dgp{id}(x,xd,p);
    case 7 % Mode change
        X = sys.pm{id}(l);
end

end

