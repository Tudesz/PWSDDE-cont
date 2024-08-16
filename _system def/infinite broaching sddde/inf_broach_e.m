function X = inf_broach_e(x,xd,p,id,type,l)
%INF_BROACH_E Event definitions for an example "infinite" broaching
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

% Select the event which occured
switch id

    % 1) Tooth entry event
    case 1 
        switch type
            % Event condition
            case 1
                X = x(6) + x(3) - 1;
            % Event condition Jacobian
            case 2
                switch l
                    % wrt x(t)
                    case 0
                        X = [0 0 1 0 0 1];
                    % wrt x(t-tau_1)
                    case 1
                        X = zeros(1,6);
                end
            % Parameter Jacobian of h
            case 3
                X = zeros(1,8);
            % Event map
            case 4
                X = x;
                X(6) = x(6) - 1;
                % X(5) = 1 + xd(5,1) - x(5);
            % Event map Jacobian
            case 5
                switch l
                    % wrt x(t)
                    case 0
                        X = eye(6);
                        % X(5,5) = -1;
                    % wrt x(t-tau_1)
                    case 1
                        X = zeros(6);
                        % X(5,5) = -1;
                end
            % Parameter Jacobian of g
            case 6
                X = zeros(6,8);
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

    % 2) Tooth exit event
    case 2 
        switch type
            % Event condition
            case 1
                X = x(6) + x(3) - p(8);
            % Event condition Jacobian
            case 2
                switch l
                    % wrt x(t)
                    case 0
                        X = [0 0 1 0 0 1];
                    % wrt x(t-tau_1)
                    case 1
                        X = zeros(1,6);
                end
            % Parameter Jacobian of h
            case 3
                X = [zeros(1,7) -1];
            % Event map
            case 4
                X = x;
            % Event map Jacobian
            case 5
                switch l
                    % wrt x(t)
                    case 0
                        X = eye(6);
                    % wrt x(t-tau_1)
                    case 1
                        X = zeros(6);
                end
            % Parameter Jacobian of g
            case 6
                X = zeros(6,8);
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
end

end

