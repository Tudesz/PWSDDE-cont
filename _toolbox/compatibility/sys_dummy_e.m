function X = sys_dummy_e(x,xd,p,id,type,l,uid)
%SYS_DUMMY_E Dummy event definition for smooth systems
% Input:
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
%   l: index of the requested delayed Jacobian J_h or J_g or in case of
%       pi_e it denotes the incoming and outcoming modes with 1 and 2
%   uid: index of state variable for whose zero crossing the dummy event is
%       defined (default last)
% Output:
%   X: appropriate evaluation of the event functions or their Jacobians

if nargin < 7
    uid = length(x); % last state variable by default
end

switch type
    case 1

        X = x(uid);
    case 2
        X = zeros(1,size(x,1));
        if l == 0
            X(uid) = 1;
        end
    case 3
        X = zeros(1,length(p));
    case 4
        X = x;
    case 5
        if l == 0
            X = eye(size(x,1));
        else
            X = zeros(size(x,1));
        end
    case 6
        X = zeros(size(x,1),length(p));
    case 7
        X = 1;
end

end