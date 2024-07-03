function X = bldosc_e(x,xd,p,id,type,l)
%BLDOSC_E Event definitions for an example delayed bilinear oscillator 
%problem (event function, event map, corresponding jacobians, mode changes)
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
%   x(3): cos(om*t)
%   x(4): sin(om*t)

% Sytem parameters
%   p(1): tau
%   p(2): k
%   p(3): zeta
%   p(4): beta
%   p(5): f
%   p(6): om


% Select the event which occured
switch id

    % 1) Establishing contact
    case 1 
        switch type
            % Event condition
            case 1
                X = x(1)-1;
            % Event condition Jacobian
            case 2
                switch l
                    % wrt x(t)
                    case 0
                        X = [1 0 0 0];
                    % wrt x(t-tau_1)
                    case 1
                        X = zeros(1,4);
                end
            % Parameter Jacobian of h
            case 3
                X = zeros(1,6);
            % Event map
            case 4
                X = x;
            % Event map Jacobian
            case 5
                switch l
                    % wrt x(t)
                    case 0
                        X = eye(4);
                    % wrt x(t-tau_1)
                    case 1
                        X = zeros(4);
                end
            % Parameter Jacobian of g
            case 6
                X = zeros(4,6);
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

    % 2) Loosing contact
    case 2 
        switch type
            % Event condition
            case 1
                X = x(1)-1;
            % Event condition Jacobian
            case 2
                switch l
                    % wrt x(t)
                    case 0
                        X = [1 0 0 0];
                    % wrt x(t-tau_1)
                    case 1
                        X = zeros(1,4);
                end
            % Parameter Jacobian of h
            case 3
                X = zeros(1,6);
            % Event map
            case 4
                X = x;
            % Event map Jacobian
            case 5
                switch l
                    % wrt x(t)
                    case 0
                        X = eye(4);
                    % wrt x(t-tau_1)
                    case 1
                        X = zeros(4);
                end
            % Parameter Jacobian of g
            case 6
                X = zeros(4,6);
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

    % 3) Grazing event
    case 3 
        switch type
            % Event condition
            case 1
                X = x(1)-1;
            % Event condition Jacobian
            case 2
                switch l
                    % wrt x(t)
                    case 0
                        X = [1 0 0 0];
                    % wrt x(t-tau_1)
                    case 1
                        X = zeros(1,4);
                end
            % Parameter Jacobian of h
            case 3
                X = zeros(1,6);
            % Event map
            case 4
                X = x;
            % Event map Jacobian
            case 5
                switch l
                    % wrt x(t)
                    case 0
                        X = eye(4);
                    % wrt x(t-tau_1)
                    case 1
                        X = zeros(4);
                end
            % Parameter Jacobian of g
            case 6
                X = zeros(4,6);
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

end