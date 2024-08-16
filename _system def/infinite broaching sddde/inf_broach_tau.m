function X = inf_broach_tau(x,p,ind,type)
%INF_BROACH_TAU Time delay definition for an example "infinite" broaching
% problem
% Input:
%   x: state of the system at t
%   p: system parameter vector
%   ind: index of requested time delay
%   type: flag for which property is requested:
%       1, type of the delayed term:
%           -> 1, state dependent point delay (normal)
%           -> 2, state dependent point delay (neutral)
%       2, value of the time delay tau(x,p)
%       3, parameter jacobian of the time delay Jp_tau(x,p)
%       4, derivatives of the delay with respect to x Jx_tau(x,p)
% Output:
%   X: appropriate evaluation of the delay or its parameter Jacobian

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

% Single state dependent point delay
switch type
    % Delay type
    case 1
        X = 1;
    % Time delays
    case 2
        X = x(5);
    % Parameter Jacobians
    case 3
        X = zeros(1,8);
    % State Jacobian
    case 4
        X = [0 0 0 0 1 0];
end

end
