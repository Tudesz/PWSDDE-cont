function X = tau_bldosc(p,type)
%TAU_IMPOSC Time delay definition for an example delayed bilinear oscillator 
% problem
% Input:
%   p: system parameter vector
%   type: flag for which property is requested:
%       1, time delay array tau(p)
%       2, parameter jacobian of time delays Jp_tau(p)
% Output:
%   X: appropriate evaluation of the delay or its parameter Jacobain

% Sytem parameters
%   p(1): tau
%   p(2): k
%   p(3): zeta
%   p(4): beta
%   p(5): f
%   p(6): om

switch type
    % Time delays
    case 1
        X = [p(1)];
    % Parameter Jacobians
    case 2
        X = [1 zeros(1,5)];
end

end