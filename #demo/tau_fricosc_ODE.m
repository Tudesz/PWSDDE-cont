function X = tau_fricosc_ODE(p,type)
%TAU_FRICOSC_ODE Time delay definition for an example dry friction oscillator 
% problem
% --- REFERENCE ODE CASE ---
% Input:
%   p: system parameter vector
%   type: flag for which property is requested:
%       1, time delay array tau(p)
%       2, parameter jacobian of time delays Jp_tau(p)
% Output:
%   X: appropriate evaluation of the delay or its parameter Jacobain

% Sytem parameters
%   p(1): zeta
%   p(2): omega
%   p(3): tau
%   p(4): f_0
%   p(5): eta

switch type
    % Time delays
    case 1
        X = [];
    % Parameter Jacobians
    case 2
        X = [];
end

end