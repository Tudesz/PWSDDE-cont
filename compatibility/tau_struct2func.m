function X = tau_struct2func(sys,p,type)
%F_STRUCT2FUNC Convert time delay definitions from old struct data to new
%function definition form
% Input:
%   sys: old system definition data structure
%       -> tau: system time delays
%       -> tau_dp: derivative of time delays wrt system parameters
%   p: system parameter vector
%   type: flag for which property is requested:
%       1, time delay array tau(p)
%       2, parameter jacobian of time delays Jp_tau(p)
% Output:
%   X: appropriate evaluation of the delay or its parameter Jacobain

% Select vector field mode

switch type
    case 1 % Time delays
        X = sys.tau(p);
    case 2 % Time delay derivatives
        X = sys.tau_dp(p);
end

end

