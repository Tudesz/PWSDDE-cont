function sys = ddebiftool_sys(sys_rhs,sys_tau,sys_deri,uid)
%DDEBIFTOOL_SYS Convert DDE-Biftool compatible function definitions to work
%with PWSDDE-cont
% Input:
%   sys_rhs: function responsible for evaluating the DDE right hand side
%   sys_tau: function name for evaluating the time delays
%   sys_deri: derivatives of the RHS for continuation
%   uid: index of state variable for whose zero crossing the dummy event is
%       defined (default last)
% Output:
%   sys: functions compatible with the new continuation framework
%     -> f: vector field definitions
%     -> e: event functions and maps
%     -> tau: time delays
%     -> mode_no: number of vector field modes
%     -> event_no: number of possible events (one dummy event)

sys.f = @(x,xd,p,type,mode,l) biftool_f(sys_rhs,sys_deri,x,xd,p,type,l);
sys.mode_no = 1;
if nargin < 4
    sys.e = @(x,xd,p,id,type,l) sys_dummy_e(x,xd,p,id,type,l);
else
    sys.e = @(x,xd,p,id,type,l) sys_dummy_e(x,xd,p,id,type,l,uid);
end
sys.event_no = 1;
tau_ind = feval(sys_tau);
sys.tau = @(p,ind,type) biftool_tau(tau_ind,p,ind,type);
sys.tau_no = length(tau_ind);

end

% Auxiliary functions

function X = biftool_f(sys_rhs,sys_deri,x,xd,p,type,l)
%BIFTOOL_F Convert vector field definitions used by dde-biftool to a new
%function definition form(rhs of the pws-dde and its Jacobians)
% Input:
%   sys_rhs: function used by dde-biftool to evaluate the DDE RHS
%   sys_tau: function used by dde-biftool to find Jacobians of the DDE RHS
%   x: state vector at time t (n)
%   xd: delayed state vectors at t-tau_l (n x L)
%   p: system parameter vector
%   type: flag for which property is requested:
%       1, The function evaluation itself f(x,xd,p)
%       2, Jacobian wrt x or xd J_x(x,xd,p)
%       3, Parameter Jacobians J_p(x,xd,p)
%   l: index of the requested delayed Jacobian J_x
% Output:
%   X: appropriate evaluation of the vector field of one of its Jacobians

xx = [x xd]; % combined state variable vector

switch type
    case 1  % Vector field
        X = feval(sys_rhs,xx,p);
    case 2 % Jacobians wrt components of xx
        X = feval(sys_deri,xx,p,l,[],[]);
    case 3 % Jacobians wrt p
        X = zeros(length(x),length(p));
        for i = 1:length(p)
            X = feval(sys_deri,xx,p,[],i,[]).';
        end
end

end

function X = biftool_tau(tau_ind,p,ind,type)
%BIFTOOL_TAU Convert time delay definitions from dde-biftool to a new
%function definition form
% Input:
%   tau_ind: indices of the delays in p
%   p: system parameter vector
%   ind: index of requested time delay
%   type: flag for which property is requested:b
%       1, type of the delayed term:
%           -> 1, fix point delay (normal)
%           -> 2, fix point delay (neutral)
%       2, value of the time delay tau(p)
%       3, parameter jacobian of the time delay Jp_tau(p)
% Output:
%   X: appropriate evaluation of the delay or its parameter Jacobain

switch type
    case 1 % Delay type
        X = 1;
    case 2 % Time delays
        X = p(tau_ind(ind));
    case 3 % Time delay derivatives
        X = zeros(1,length(p));
        X(tau_ind(ind)) = 1;
end

end

