function plot_orb_events(orb,sys,res)
%PLOT_ORB_EVENTS Plot event function evaluations over a periodic orbit
% Input:
%   orb: periodic orbit data structure
%    -> sig: sloution signature (event list)
%    -> U: state variable vector (M*N*n)
%    -> T: segment lengths (N)
%    -> p: parameter vector
%    -> n: number of degrees of freedom
%    -> M: Chebyshev mesh resolution
%   sys: names of the functions that define the system
%    -> e: event function, map and corresponding Jacobians
%    -> tau: time delay and its parameter Jacobians
%    -> event_no: number of disctinc events
%    -> tau_no: number of distinct time delays
%   res: interpolation resolution (default M)
%   fs: functions for plotting, default: fs1(t,y)=t wrt fs2(t,y)=y
%   res: output interpolation resolution (default 100)

% Function definitions
h_ej = @(j,x,xd) feval(sys.e,x,xd,orb.p,j,1,0);  % event condition at ej

% Initialization
if nargin<3
    res = 100;
end
[t0,x0] = bvp2sig(orb.U,orb.T,orb.M);
del0 = po_delay_interp(orb.U,orb.T,orb.p,orb.M,sys);
xd0 = del0.ud;
[t1,x1] = bvp2sig(orb.U,orb.T,orb.M,res);
del1 = po_delay_interp(orb.U,orb.T,orb.p,orb.M,sys,res);
xd1 = del1.ud;

% evaluate and plot event functions
yline(0,':k'); hold on
for i=1:sys.event_no
    ev1 = zeros(1,length(t1));
    for j = 1:length(t1)
        ev1(j) = h_ej(i,x1(:,j),squeeze(xd1(:,:,j)));
    end
    plot(t1,ev1)
end
for i=1:sys.event_no
    ev0 = zeros(1,length(t0));
    for j = 1:length(t0)
        ev0(j) = h_ej(i,x0(:,j),squeeze(xd0(:,:,j)));
    end
    plot(t0,ev0,'.k')
end
xlabel('$t$'); ylabel('$h_i(t)$')
hold off;
box on;

end

