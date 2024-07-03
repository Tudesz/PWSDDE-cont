function X = nlne_rob_tau(p,ind,type)
%NLNE_ROB_TAU Time delay definition for an example nonlinear robotic arm
% Input:
%   p: system parameter vector
%   ind: index of the requested time delay
%   type: flag for which property is requested:
%       1, type of the delayed term:
%           -> 1, fix point delay (normal)
%           -> 2, fix point delay (neutral)
%       2, value of the time delay tau(p)
%       3, parameter jacobian of the time delay Jp_tau(p)
% Output:
%   X: appropriate evaluation of the delay or its parameter Jacobian

% State variables
%   x(1): x1(t)
%   x(2): x1'(t)
%   x(3): x2(t)
%   x(4): x2'(t)

% Sytem parameters
%   p(1): khi
%   p(2): r
%   p(3): gamma
%   p(4): mu
%   p(5): k
%   p(6): knl
%   p(7): tau

% Select the requested time delay
switch ind

    % Single constant neutral point delay
    case 1
        switch type
            % Delay type
            case 1
                X = 2;
            % Time delays
            case 2
                X = p(7);
            % Parameter Jacobians
            case 3
                X = [zeros(1,6) 1];
        end

end