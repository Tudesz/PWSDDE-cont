function X = nlne_rob_f(x,xd,p,mode,type,l)
%NLNE_ROB_F Vector field definitions for an example nonlinear robotic arm 
% (rhs of the pws-dde and its Jacobians)
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
%   x(1): x1(t)
%   x(2): x1'(t)
%   x(3): x2(t)
%   x(4): x2'(t)

% Sytem parameters
%   p(1): chi
%   p(2): r
%   p(3): gamma
%   p(4): mu
%   p(5): k
%   p(6): knl
%   p(7): tau

switch mode

    % 1) No contact mode
    case 1 
        switch type
            % Vector field
            case 1
                X = [x(2);
                -2*p(1)*p(2)*p(3)*(x(2)-x(4))...
                - p(3)^2*p(2)*(x(1)-x(3))-x(1)...
                - p(4)*p(2)*(x(1)-x(3))^3+p(5)*xd(2,1)+p(6)*xd(2,1)^3;
                x(4);
                -2*p(1)*p(3)*(x(4)-x(2))-p(3)^2*(x(3)-x(1))...
                - p(4)*(x(3)-x(1))^3];
            % Jacobian matrices
            case 2
                switch l
                    % wrt x(t)
                    case 0
                        X = [0, 1, 0, 0; ...
                        -p(3)^2*p(2)-1-3*p(2)*p(4)*(x(1)-x(3))^2,...
                        -2*p(1)*p(2)*p(3),...
                        p(2)*p(3)^2+3*p(2)*p(4)*(x(1)-x(3))^2,...
                        2*p(1)*p(2)*p(3); ...
                        0, 0, 0, 1; ...
                        p(3)^2+3*p(4)*(-x(1)+x(3))^2, 2*p(1)*p(3), ...
                        -p(3)^2-3*p(4)*(-x(1)+x(3))^2, -2*p(1)*p(3)];
                    % wrt x(t-tau_1)
                    case 1
                        X = zeros(4);
                        X(2,2)=p(5)+3*p(6)*xd(2,1)^2;
                end
            % Parameter Jacobian
            case 3
                X = [[0, -2*p(2)*p(3)*(x(2)-x(4)), 0,...
                    -2*p(3)*(x(4)-x(2))].'...
                    [0, (-p(3)^2)*(x(1)-x(3))...
                    - p(4)*(x(1)-x(3))^3-2*p(1)*p(3)*(x(2)...
                    -x(4)), 0, 0].'...
                    [0, -2*p(2)*p(3)*(x(1)-x(3))...
                    - 2*p(1)*p(2)*(x(2)-x(4)), 0, ...
                    - 2*p(3)*(-x(1)+x(3))...
                    - 2*p(1)*(-x(2)+x(4))].'...
                    [0, (-p(2))*(x(1)-x(3))^3, ...
                    0, -(-x(1)+x(3))^3].'...
                    [0, xd(2,1), 0 , 0].' ...
                    [0, xd(2,1)^3, 0, 0].'...
                    [0, 0, 0, 0].'];
        end
end

end

