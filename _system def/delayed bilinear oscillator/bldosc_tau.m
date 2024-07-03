function X = bldosc_tau(p,ind,type)
%BLDOSC_TAU Time delay definition for an example delayed bilinear oscillator 
% problem
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

% Sytem parameters
%   p(1): tau
%   p(2): k
%   p(3): zeta
%   p(4): beta
%   p(5): f
%   p(6): om

% Select the requested time delay
switch ind

    % Single point delay present
    case 1
        switch type
            % Delay type
            case 1
                X = 1;
            % Time delays
            case 2
                X = p(1);
            % Parameter Jacobians
            case 3
                X = [1 zeros(1,5)];
        end

end