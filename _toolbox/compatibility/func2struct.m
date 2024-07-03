function struct = func2struct(sys,dims,type)
%STRUCT2FUNC Convert new function form to old system definition structure
% Input:
%   sys: functions compatible with the new continuation framework
%     -> f: vector field definitions
%     -> e: event functions and maps
%     -> tau: time delays
%     -> mode_no: number of vector field modes
%     -> event_no: number of possible events
%   dims: orbit dimensions [M n]
%   type: problem type: 1) smooth ODE 2) non-smooth ODE 
%       3) switched non-autonomous DDE 4) non-smooth DDE
%       (default 4)
% Output:
%   struct: old system definition data structure
%     -> f{m}: governing autonomous vector spaces R^n -> R^n (DDE rhs)
%     -> h{i}: switching or event surfaces R^n -> R
%     -> g{i}: event maps R^n->R^n
%     -> Ju{m}: vector space Jacobians R^n -> R^(n x n)(DDE rhs)
%     -> Jud{m,k}: vector space Jacobians wrt x(t-tau(k))
%     -> dh{i}: event surface Jacobians R^n -> R^(1 x n)
%     -> dg{i}: event map Jacobians R^n -> R^(n x n)
%     -> dhd{i,k}: event surface Jacobians wrt x(t-tau(k))
%     -> dgd{i,k}: event map Jacobians wrt x(t-tau(k))
%     -> Jp{m}: Jacobian of f R^n -> R^(n x l)
%     -> dhp{i}: Jacobian of h R^n -> R^(1 x l)
%     -> dgp{i}: Jacobian of g R^n -> R^(n x l)
%     -> tau: system time delays
%     -> tau_dp: derivative of time delays wrt system parameters
%     -> pm{i}: mode changes on event surfaces

if nargin<3
    type = 4;   % PWSDDE definition assumed by default
end

% auxiliary data
struct.dims = dims;
struct.type = type;

% vector field definitions
struct.f = cell(1,sys.mode_no);
struct.Ju = cell(1,sys.mode_no);
struct.Jud = cell(sys.mode_no,sys.tau_no);
struct.Jp = cell(1,sys.mode_no);
for i=1:sys.mode_no
    struct.f{i} = @(x,xd,p) feval(sys.f,x,xd,p,i,1,0);
    struct.Ju{i} = @(x,xd,p) feval(sys.f,x,xd,p,i,2,0);
    struct.Jp{i} = @(x,xd,p) feval(sys.f,x,xd,p,i,3,0);
    for j = sys.tau_no
        struct.Jud{i,j} = @(x,xd,p) feval(sys.f,x,xd,p,i,2,j);
    end
end

% event map definitions
struct.g = cell(1,sys.event_no);
struct.dg = cell(1,sys.event_no);
struct.dgd = cell(sys.event_no,sys.tau_no);
struct.dgp = cell(1,sys.event_no);
for i=1:sys.event_no
    struct.g{i} = @(x,xd,p) feval(sys.e,x,xd,p,i,4,0);
    struct.dg{i} = @(x,xd,p) feval(sys.e,x,xd,p,i,5,0);
    struct.dg{i} = @(x,xd,p) feval(sys.e,x,xd,p,i,6,0);
    for j = sys.tau_no
        struct.dgd{i,j} = @(x,xd,p) feval(sys.e,x,xd,p,i,5,j);
    end
end

% event function definitions
struct.h = cell(1,sys.event_no);
struct.dh = cell(1,sys.event_no);
struct.dhd = cell(sys.event_no,sys.tau_no);
struct.dhp = cell(1,sys.event_no);
for i=1:sys.event_no
    struct.h{i} = @(x,xd,p) feval(sys.e,x,xd,p,i,1,0);
    struct.dh{i} = @(x,xd,p) feval(sys.e,x,xd,p,i,2,0);
    struct.dhp{i} = @(x,xd,p) feval(sys.e,x,xd,p,i,3,0);
    for j = sys.tau_no
        struct.dgd{i,j} = @(x,xd,p) feval(sys.e,x,xd,p,i,2,j);
    end
end

% time delay definitions
struct.tau = @(p) f_tau(sys,p);
struct.tau_dp = @(p) f_tau_dp(sys,p);

% mode change identifiers
struct.pm = cell(1,sys.event_no);
for i = 1:sys.event_no
    struct.pm{i} = [feval(sys.e,[],[],[],i,7,1) feval(sys.e,[],[],[],i,7,2)]; 
end

% Auxiliary functions
function tau = f_tau(sys,p)
    tau = zeros(1,sys.tau_no);
    for k = 1:tau_no
        tau(k) = feval(sys.tau,p,k,2);
    end
end
function tau_dp = f_tau_dp(sys,p)
    tau_dp = zeros(sys.tau_no,length(p));
    for k = 1:tau_no
        tau_dp(k,:) = feval(sys.tau,p,k,3);
    end
end

end
