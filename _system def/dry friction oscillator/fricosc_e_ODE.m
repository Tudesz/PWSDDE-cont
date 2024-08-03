function X = fricosc_e_ODE(x,xd,p,id,type,l)
%FRICOSC_E_ODE Event definitions for an example dry friction oscillator
% problem in an ODE form for validation
% (event function, event map, corresponding jacobians, mode changes)
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
% Output:
%   X: appropriate evaluation of the event functions or their Jacobians

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


% Select the event which occured
switch id

    % 1) Velocity reversal while sliding (1->2)
    case 1 
        switch type
            % Event condition
            case 1
                X = x(2);
            % Event condition Jacobian
            case 2
                switch l
                    % wrt x(t)
                    case 0
                        X = [0, 1, 0];
                end
            % Parameter Jacobian of h
            case 3
                X = zeros(1,5);
            % Event map
            case 4
                X = x;
            % Event map Jacobian
            case 5
                switch l
                    % wrt x(t)
                    case 0
                        X = eye(3);
                end
            % Parameter Jacobian of g
            case 6
                X = zeros(3,5);
            % Mode change
            case 7
               switch l
                    % Incoming mode
                    case 1
                        X = 1;
                    % Outgoing mode
                    case 2
                        X = 2;
                end
        end

    % 2) Velocity reversal while sliding (2->1)
    case 2 
        switch type
            % Event condition
            case 1
                X = x(2);
            % Event condition Jacobian
            case 2
                switch l
                    % wrt x(t)
                    case 0
                        X = [0, 1, 0];
                end
            % Parameter Jacobian of h
            case 3
                X = zeros(1,5);
            % Event map
            case 4
                X = x;
            % Event map Jacobian
            case 5
                switch l
                    % wrt x(t)
                    case 0
                        X = eye(3);
                end
            % Parameter Jacobian of g
            case 6
                X = zeros(3,5);
            % Mode change
            case 7
               switch l
                    % Incoming mode
                    case 1
                        X = 2;
                    % Outgoing mode
                    case 2
                        X = 1;
                end
        end

    % 3) Transition to sticking (1->3)
    case 3 
        switch type
            % Event condition
            case 1
                X = x(2);
            % Event condition Jacobian
            case 2
                switch l
                    % wrt x(t)
                    case 0
                        X = [0, 1, 0];
                end
            % Parameter Jacobian of h
            case 3
                X = zeros(1,5);
            % Event map
            case 4
                X = x;
            % Event map Jacobian
            case 5
                switch l
                    % wrt x(t)
                    case 0
                        X = eye(3);
                end
            % Parameter Jacobian of g
            case 6
                X = zeros(3,5);
            % Mode change
            case 7
               switch l
                    % Incoming mode
                    case 1
                        X = 1;
                    % Outgoing mode
                    case 2
                        X = 3;
                end
        end

    % 4) Separation from sticking (3->1)
    case 4
        switch type
            % Event condition
            case 1
                X = -x(1) - 2*p(1)*x(2) + p(4)*cos(x(3)-p(2)*p(3)) - p(5);
            % Event condition Jacobian
            case 2
                switch l
                    % wrt x(t)
                    case 0
                        X = [-1, -2*p(1), -p(4)*sin(x(3)-p(2)*p(3))];
                end
            % Parameter Jacobian of h
            case 3
                X = [-2*x(2), p(4)*sin(x(3)-p(2)*p(3)),...
                    p(4)*sin(x(3)-p(2)*p(3)), cos(x(3)-p(2)*p(3)), -1];
            % Event map
            case 4
                X = x;
            % Event map Jacobian
            case 5
                switch l
                    % wrt x(t)
                    case 0
                        X = eye(3);
                end
            % Parameter Jacobian of g
            case 6
                X = zeros(3,5);
            % Mode change
            case 7
               switch l
                    % Incoming mode
                    case 1
                        X = 3;
                    % Outgoing mode
                    case 2
                        X = 1;
                end
        end

    % 5) Transition to sticking (2->3)
    case 5 
        switch type
            % Event condition
            case 1
                X = x(2);
            % Event condition Jacobian
            case 2
                switch l
                    % wrt x(t)
                    case 0
                        X = [0, 1, 0];
                end
            % Parameter Jacobian of h
            case 3
                X = zeros(1,5);
            % Event map
            case 4
                X = x;
            % Event map Jacobian
            case 5
                switch l
                    % wrt x(t)
                    case 0
                        X = eye(3);
                end
            % Parameter Jacobian of g
            case 6
                X = zeros(3,5);
            % Mode change
            case 7
                switch l
                    % Incoming mode
                    case 1
                        X = 2;
                    % Outgoing mode
                    case 2
                        X = 3;
                end
        end

    % 6) Separation from sticking (3->2)
    case 6
        switch type
            % Event condition
            case 1
                X = -x(1) - 2*p(1)*x(2) + p(4)*cos(x(3)-p(2)*p(3)) + p(5);
            % Event condition Jacobian
            case 2
                switch l
                    % wrt x(t)
                    case 0
                        X = [-1, -2*p(1), -p(4)*sin(x(3)-p(2)*p(3))];
                end
            % Parameter Jacobian of h
            case 3
                X = [-2*x(2), p(4)*sin(x(3)-p(2)*p(3)),...
                    p(4)*sin(x(3)-p(2)*p(3)), cos(x(3)-p(2)*p(3)), 1];
            % Event map
            case 4
                X = x;
            % Event map Jacobian
            case 5
                switch l
                    % wrt x(t)
                    case 0
                        X = eye(3);
                end
            % Parameter Jacobian of g
            case 6
                X = zeros(3,5);
            % Mode change
            case 7
               switch l
                    % Incoming mode
                    case 1
                        X = 3;
                    % Outgoing mode
                    case 2
                        X = 2;
                end
        end

    % 7) Phase reset in forward sliding (1->1)
    case 7 
        switch type
            % Event condition
            case 1
                X = x(3)-2*pi;
            % Event condition Jacobian
            case 2
                switch l
                    % wrt x(t)
                    case 0
                        X = [0, 0, 1];
                end
            % Parameter Jacobian of h
            case 3
                X = zeros(1,5);
            % Event map
            case 4
                X = [x(1); x(2); 0];
            % Event map Jacobian
            case 5
                switch l
                    % wrt x(t)
                    case 0
                        X = [1, 0, 0;
                            0, 1, 0;
                            0, 0, 0];
                end
            % Parameter Jacobian of g
            case 6
                X = zeros(3,5);
            % Mode change
            case 7
               switch l
                    % Incoming mode
                    case 1
                        X = 1;
                    % Outgoing mode
                    case 2
                        X = 1;
                end
        end

    % 8) Phase reset in backward sliding (2->2)
    case 8
        switch type
            % Event condition
            case 1
                X = x(3)-2*pi;
            % Event condition Jacobian
            case 2
                switch l
                    % wrt x(t)
                    case 0
                        X = [0, 0, 1];
                end
            % Parameter Jacobian of h
            case 3
                X = zeros(1,5);
            % Event map
            case 4
                X = [x(1); x(2); 0];
            % Event map Jacobian
            case 5
                switch l
                    % wrt x(t)
                    case 0
                        X = [1, 0, 0;
                            0, 1, 0;
                            0, 0, 0];
                end
            % Parameter Jacobian of g
            case 6
                X = zeros(3,5);
            % Mode change
            case 7
               switch l
                    % Incoming mode
                    case 1
                        X = 2;
                    % Outgoing mode
                    case 2
                        X = 2;
                end
        end

    % 9) Phase reset in sticking (3->3)
    case 9 
        switch type
            % Event condition
            case 1
                X = x(3)-2*pi;
            % Event condition Jacobian
            case 2
                switch l
                    % wrt x(t)
                    case 0
                        X = [0, 0, 1];
                end
            % Parameter Jacobian of h
            case 3
                X = zeros(1,5);
            % Event map
            case 4
                X = [x(1); x(2); 0];
            % Event map Jacobian
            case 5
                switch l
                    % wrt x(t)
                    case 0
                        X = [1, 0, 0;
                            0, 1, 0;
                            0, 0, 0];
                end
            % Parameter Jacobian of g
            case 6
                X = zeros(3,5);
            % Mode change
            case 7
               switch l
                    % Incoming mode
                    case 1
                        X = 3;
                    % Outgoing mode
                    case 2
                        X = 3;
                end
        end
end

end