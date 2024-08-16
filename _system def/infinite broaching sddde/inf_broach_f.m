function X = inf_broach_f(x,xd,p,mode,type,l)
%INF_BROACH_F Vector field definitions for an example "infinite" broaching
% problem (rhs of the pws-dde and its Jacobians)
% Input:
%   x: state vector at time t (n)
%   xd: delayed state vectors at t-tau_l (n x L)
%   p: system parameter vector
%   mode: index of the vector field mode
%   type: flag for which property is requested:
%       1, The function evaluation itself f(x,xd,p)
%       2, Jacobian wrt x or xd J_x(x,xd,p) or J_xd(x,xd,p)
%       3, Parameter Jacobians J_p(x,xd,p)
%   l: index of the requested delayed Jacobian J_x or J_xd
% Output:
%   X: appropriate evaluation of the vector field or one of its Jacobians

% State variables
%   x(1): y(t)
%   x(2): y'(t)
%   x(3): z(t)
%   x(4): z'(t)
%   x(5): tau(t)
%   x(6): t

% Dimensionless system paramters
% p(1): c1
% p(2): k1
% p(3): c2
% p(4): k2
% p(5): w0
% p(6): xi
% p(7): nt
% p(8): d0

% Determine the number of teeth in contact
if mode == 1
    n = p(7) + 1; % maximum number of teeth cutting
elseif mode == 2
    n = p(7); % minimum number of teeth cutting
end


switch type
    % Vector field
    case 1
        fc = n*p(5)*(1 + xd(1,1) - x(1));
        X = [x(2); -p(1)*x(2) - p(2)*x(1) + p(6)*fc;
            x(4); -p(3)*x(4) - p(4)*x(3) - fc;
            xd(4,1) - x(4); 1];
    % Jacobian matrices
    case 2
        switch l
            % wrt x(t)
            case 0
                X = [0 1 0 0 0 0; -p(2)-p(6)*n*p(5) -p(1) 0 0 0 0;
                    0 0 0 1 0 0; n*p(5) 0 -p(4) -p(3) 0 0;
                    0 0 0 -1 0 0; zeros(1,6)];
            % wrt x(t-tau_1)
            case 1
                X = [zeros(1,6); p(6)*n*p(5) zeros(1,5);
                    zeros(1,6); -n*p(5) zeros(1,5);
                    0 0 0 1 0 0; zeros(1,6)];
        end
    % Parameter Jacobian
    case 3
        X = [zeros(1,8); 
            -x(2) -x(1) 0 0 p(6)*n*(1 + xd(1,1) - x(1)) ...
            p(5)*n*(1 + xd(1,1) - x(1)) p(5)*p(6)*(1 + xd(1,1) - x(1)) 0;
            zeros(1,8); 
            0 0 -x(4) -x(3) -n*(1 + xd(1,1) - x(1)) 0 ...
            -p(5)*(1 + xd(1,1) - x(1)) 0;
            zeros(2,8)];
end


end

