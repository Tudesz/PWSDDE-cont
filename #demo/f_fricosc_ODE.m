function X = f_fricosc_ODE(x,xd,p,mode,type,l)
%F_FRICOSC_ODE Vector field definitions for an example dry friction oscillator 
% problem (rhs of the pws-dde and its Jacobians)
% Input:
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

% State variables
%   x(1): x(t)
%   x(2): x'(t)
%   x(3): phi(t)

% Sytem parameters
%   p(1): zeta
%   p(2): omega
%   p(3): tau
%   p(4): f_0
%   p(5): eta


% Select vector field mode
switch mode

    % 1) Sliding with x'(t)>0
    case 1 
        switch type
            % Vector field
            case 1
                X = [x(2); 
                    -x(1) - 2*p(1)*x(2) + p(4)*cos(x(3)-p(2)*p(3)) - 1;
                    p(2)];
            % Jacobian matrices
            case 2
                switch l
                    % wrt x(t)
                    case 0
                        X = [0, 1, 0;
                            -1, -2*p(1), -p(4)*sin(x(3)-p(2)*p(3));
                            0, 0, 0];
                end
            % Parameter Jacobian
            case 3
                X = [zeros(1,5);
                    -2*x(2), p(3)*p(4)*sin(x(3)-p(2)*p(3)), ...
                    p(2)*p(4)*sin(x(3)-p(2)*p(3)), cos(x(3)-p(2)*p(3)), 0;
                    0, 1, 0, 0, 0];
        end

    % 2) Sliding with x'(t)<0
    case 2 
        switch type
            % Vector field
            case 1
                X = [x(2); 
                    -x(1) - 2*p(1)*x(2) + p(4)*cos(x(3)-p(2)*p(3)) + 1;
                    p(2)];
            % Jacobian matrices
            case 2
                switch l
                    % wrt x(t)
                    case 0
                        X = [0, 1, 0;
                            -1, -2*p(1), -p(4)*sin(x(3)-p(2)*p(3));
                            0, 0, 0];
                end
            % Parameter Jacobian
            case 3
                X = [zeros(1,5);
                    -2*x(2), p(3)*p(4)*sin(x(3)-p(2)*p(3)), ...
                    p(2)*p(4)*sin(x(3)-p(2)*p(3)), cos(x(3)-p(2)*p(3)), 0;
                    0, 1, 0, 0, 0];
        end

    % 3) Sticking with x'(t)=0
    case 3
        switch type
            % Vector field
            case 1
                X = [x(2); 
                    0;
                    p(2)];
            % Jacobian matrices
            case 2
                switch l
                    % wrt x(t)
                    case 0
                        X = [0, 1, 0;
                            zeros(2,3)];
                end
            % Parameter Jacobian
            case 3
                X = [zeros(2,5);
                    0, 1, 0, 0, 0];
        end
end

end