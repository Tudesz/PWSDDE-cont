function X = bldosc_f(x,xd,p,mode,type,l)
%BLDOSC_F Vector field definitions for an example delayed bilinear 
%oscillator problem (rhs of the pws-dde and its Jacobians)
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
%   x(1): x(t)
%   x(2): x'(t)
%   x(3): cos(om*t)
%   x(4): sin(om*t)

% Sytem parameters
%   p(1): tau
%   p(2): k
%   p(3): zeta
%   p(4): beta
%   p(5): f
%   p(6): om


% Select vector field mode
switch mode

    % 1) No contact mode
    case 1 
        switch type
            % Vector field
            case 1
                X = [x(2); 
                    -2*p(3)*x(2) - x(1) + p(5)*x(3) + p(2)*(xd(2,1)-x(2));
                    -p(6)*x(4) + (1-x(3)^2-x(4)^2)*x(3);
                    p(6)*x(3) + (1-x(3)^2-x(4)^2)*x(4)];
            % Jacobian matrices
            case 2
                switch l
                    % wrt x(t)
                    case 0
                        X = [0 1 0 0;
                            -1 -2*p(3)-p(2) p(5) 0;
                            0 0 1-3*x(3)^2-x(4)^2 -p(6)-2*x(3)*x(4);
                            0 0 p(6)-2*x(3)*x(4) 1-x(3)^2-3*x(4)^2];
                    % wrt x(t-tau_1)
                    case 1
                        X = [zeros(1,4); 
                            0 p(2) 0 0;
                            zeros(2,4)];
                end
            % Parameter Jacobian
            case 3
                X = [zeros(1,6);
                0, xd(2,1)-x(2), -2*x(2) 0, x(3), 0;
                zeros(1,5), -x(4);
                zeros(1,5), x(3)];
        end

    % 2) Contact mode
    case 2 
        switch type
            % Vector field
            case 1
                X = [x(2);...
                -2*p(3)*x(2) - x(1) + p(5)*x(3) + ...
                p(2)*(xd(2,1)-x(2)) - p(4)*(x(1)-1);
                -p(6)*x(4) + (1-x(3)^2-x(4)^2)*x(3);
                p(6)*x(3) + (1-x(3)^2-x(4)^2)*x(4)];
            % Jacobian matrices
            case 2
                switch l
                    % wrt x(t)
                    case 0
                        X = [0 1 0 0;
                            -1-p(4) -2*p(3)-p(2) p(5) 0;
                            0 0 1-3*x(3)^2-x(4)^2 -p(6)-2*x(3)*x(4);
                            0 0 p(6)-2*x(3)*x(4) 1-x(3)^2-3*x(4)^2];
                    % wrt x(t-tau_1)
                    case 1
                        X = [zeros(1,4); 
                            0 p(2) 0 0;
                            zeros(2,4)];
                end
            % Parameter Jacobian
            case 3
                X = [zeros(1,6);
                0 xd(2,1)-x(2) -2*x(2) -(x(1)-1) x(3) 0;
                zeros(1,5) -x(4);
                zeros(1,5) x(3)];
        end

end

end