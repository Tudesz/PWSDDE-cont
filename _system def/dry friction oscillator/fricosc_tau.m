function X = fricosc_tau(p,ind,type)
%FRICOSC_TAU Time delay definition for an example dry friction oscillator 
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
%   p(1): zeta
%   p(2): omega
%   p(3): tau
%   p(4): f_0
%   p(5): eta


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
                X = p(3);
            % Parameter Jacobians
            case 3
                X = [0, 0, 1, 0, 0];
         end

end