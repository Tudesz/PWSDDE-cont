function X = f_struct2func(sys,x,xd,p,mode,type,l)
%F_STRUCT2FUNC Convert vector field definitions from old struct data to new
%function definition form
% (rhs of the pws-dde and its Jacobians)
% Input:
%   sys: old system definition data structure
%       -> f{m}: governing autonomous vector spaces R^n -> R^n (DDE rhs)
%       -> Ju{m}: vector space Jacobians R^n -> R^(n x n)(DDE rhs)
%       -> Jud{m,k}: vector space Jacobians wrt x(t-tau(k))
%       -> Jp{m}: Jacobian of f R^n -> R^(n x l)
%   x: state vector at time t (n)
%   xd: delayed state vectors at t-tau_l (n x L)
%   p: system parameter vector
%   m: mode index of the vector field mode
%   type: flag for which property is requested:
%       1, The function evaluation itself f(x,xd,p)
%       2, Jacobian wrt x or xd J_x(x,xd,p)
%       3, Parameter Jacobians J_p(x,xd,p)
%   l: index of the requested delayed Jacobian J_x
% Output:
%   X: appropriate evaluation of the vector field of one of its Jacobians

switch type
    case 1 % Vector field
        X = sys.f{mode}(x,xd,p);
    case 2 % Jacobian matricies
        if l>0
            X = sys.Jud{mode,l}(x,xd,p);
        else
            X = sys.Ju{mode}(x,xd,p);
        end
    case 3 % Parameter Jacobian
        X = sys.Jp{mode}(x,xd,p);
end

end

