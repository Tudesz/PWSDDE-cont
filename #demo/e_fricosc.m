function X = e_fricosc(x,xd,p,id,type,l)
%E_FRICOSC Event definitions for an example dry friction oscillator problem
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
%   l: index of the requested delayed Jacobian J_h or J_g
% Output:
%   X: appropriate evaluation of the event functions of Jacobians

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

    % 1) Velocity reversal while sliding
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
                    % wrt x(t-tau_1)
                    case 1
                        X = zeros(1,3);
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
                    % wrt x(t-tau_1)
                    case 1
                        X = zeros(3);
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

    % 2) Velocity reversal while sliding
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
                    % wrt x(t-tau_1)
                    case 1
                        X = zeros(1,3);
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
                    % wrt x(t-tau_1)
                    case 1
                        X = zeros(3);
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

    % 3) Transition to sticking
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
                    % wrt x(t-tau_1)
                    case 1
                        X = zeros(1,3);
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
                    % wrt x(t-tau_1)
                    case 1
                        X = zeros(3);
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

    % 4) Separation from sticking
    case 4
        switch type
            % Event condition
            case 1
                X = -x(1) - 2*p(1)*x(2) + p(4)*cos(xd(3,1)) - p(5);
            % Event condition Jacobian
            case 2
                switch l
                    % wrt x(t)
                    case 0
                        X = [-1, -2*p(1), 0];
                    % wrt x(t-tau_1)
                    case 1
                        X = [0, 0, -p(4)*sin(xd(3,1))];
                end
            % Parameter Jacobian of h
            case 3
                X = [-2*x(2), 0, 0, cos(xd(3,1)), -1];
            % Event map
            case 4
                X = x;
            % Event map Jacobian
            case 5
                switch l
                    % wrt x(t)
                    case 0
                        X = eye(3);
                    % wrt x(t-tau_1)
                    case 1
                        X = zeros(3);
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

    % 5) Transition to sticking
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
                    % wrt x(t-tau_1)
                    case 1
                        X = zeros(1,3);
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
                    % wrt x(t-tau_1)
                    case 1
                        X = zeros(3);
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

    % 6) Separation from sticking
    case 6
        switch type
            % Event condition
            case 1
                X = -x(1) - 2*p(1)*x(2) + p(4)*cos(xd(3,1)) + p(5);
            % Event condition Jacobian
            case 2
                switch l
                    % wrt x(t)
                    case 0
                        X = [-1, -2*p(1), 0];
                    % wrt x(t-tau_1)
                    case 1
                        X = [0, 0, -p(4)*sin(xd(3,1))];
                end
            % Parameter Jacobian of h
            case 3
                X = [-2*x(2), 0, 0, cos(xd(3,1)), 1];
            % Event map
            case 4
                X = x;
            % Event map Jacobian
            case 5
                switch l
                    % wrt x(t)
                    case 0
                        X = eye(3);
                    % wrt x(t-tau_1)
                    case 1
                        X = zeros(3);
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

    % 7) Phase reset in forward sliding
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
                    % wrt x(t-tau_1)
                    case 1
                        X = zeros(1,3);
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
                    % wrt x(t-tau_1)
                    case 1
                        X = zeros(3);
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

    % 7) Phase reset in forward sliding
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
                    % wrt x(t-tau_1)
                    case 1
                        X = zeros(1,3);
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
                    % wrt x(t-tau_1)
                    case 1
                        X = zeros(3);
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

    % 9) Phase reset in sticking
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
                    % wrt x(t-tau_1)
                    case 1
                        X = zeros(1,3);
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
                    % wrt x(t-tau_1)
                    case 1
                        X = zeros(3);
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