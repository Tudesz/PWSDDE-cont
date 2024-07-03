function [gr,int,ind] = po_graze_det(y0,y1,orb,sys,opts,bifs)
%PO_GRAZE_DET Grazing event detection in continuation of periodic orbits based
%on changes in the signs of event functions
% Input:
%   y0: previous state vector with its bifurcation parameter [u0; Ti0; pi0]
%   y1: current state vector with its bifurcation parameter [u0; Ti0; pi0]
%   orb: periodic orbit data structure (only the metadata part is used)
%    -> sig: sloution signature (event list)
%    -> n: number of degrees of freedom
%    -> M: Chebyshev mesh resolution
%   sys: names of the functions that define the system
%    -> f: vector field and its Jacobians
%    -> e: event function, map and corresponding Jacobians
%    -> tau: time delay and its parameter Jacobian
%   opts: solver options
%    -> pi: index of continuation parameter
%    -> gr_tol: tolerance for grazing event detection (default 1e-10)
%   bifs: extra data for two parameter continuation (optional input)
%    -> type: 1) grazing 2) sliding
%    -> ind: index of bifurcation event in the solution signature
% Output:
%   gr: bolean flag for grazing bifurcations
%   int: if true: interior grazing, else: exterior grazing
%   ind: index of h which underwent a grazing event or its position in sig
%   for external grazing events

if ~isfield(opts,'gr_tol')
    opts.gr_tol = 1e-10; % default detection tolerance limit
end

% Initialization
gr = false;             % no bifurcation by default
int = [];               % no bifurcation type by default
ind = [];               % no index by default
M = orb.M;              % mesh resolution
N = length(orb.sig);    % number of segments
m = sys.event_no;       % number of contact surfaces

% Function definitions
h_ej = @(j,x,xd,p) feval(sys.e,x,xd,p,j,1,0);  % event condition at ej

% Unpack solution vectors
lp = length(opts.pi);
p0 = pl_insert(orb.p,y0(end-lp+1:end),opts.pi);
Ti0 = y0(end-N-lp+1:end-lp);
p1 = pl_insert(orb.p,y1(end-lp+1:end),opts.pi);
Ti1 = y1(end-N-lp+1:end-lp);

% Evaluate current and delayed terms
[~,x0] = bvp2sig(y0(1:end-N-lp),Ti0,M);
[~,x1] = bvp2sig(y1(1:end-N-lp),Ti1,M);
[x0_tau,~,~,~] = po_delay_interp(y0(1:end-N-lp),Ti0,p0,M,sys);
[x1_tau,~,~,~] = po_delay_interp(y1(1:end-N-lp),Ti1,p1,M,sys);

% Differentiation matricies (for n=1)
D0 = zeros(M*N);
D1 = zeros(M*N);
DM = cheb_diff(M);
for i = 1:N
     ii = (i-1)*M+1:i*M; % indicies in D
     D0(ii,ii) = 1/Ti0(i)*DM;
     D1(ii,ii) = 1/Ti1(i)*DM;
end

% Evaluate h{i} for all points of x0 and x1
h0 = zeros(M*N,m); h1 = h0;
for i = 1:M*N
    for j = 1:m
        h0(i,j) = h_ej(j,x0(:,i),squeeze(x0_tau(:,:,i)),p0);
        h1(i,j) = h_ej(j,x1(:,i),squeeze(x1_tau(:,:,i)),p1);
    end
end

% Exterior grazing event conditions (based on dh(t))
dh1 = D1*h1(:,orb.sig);
dh0 = D0*h0(:,orb.sig);
ext_gr = diag(dh0(M:M:end,:)).*diag(dh1(M:M:end,:));
if any(ext_gr<opts.gr_tol,'all')
    [~,ind] = min(min(ext_gr));
    if lp == 1 || ind ~= bifs.ind
        % Omit false grazing events when following bifurcations
        gr = true;
        int = false;
    end
end

% Interior grazig (has higher priority than exterior grazing)
% TODO: option for using more points in the detection algorithm
int_gr = h0.*h1;
int_gr([1:M:N*M M:M:N*M],:) = [];
if any(int_gr<-opts.gr_tol,'all')
    gr = true;
    int = true;
    [~,ind] = min(min(int_gr));
end     

end

