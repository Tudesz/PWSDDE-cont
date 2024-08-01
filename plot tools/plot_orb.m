function plot_orb(orb,sys,fs,res)
%PLOT_ORB Plot the state or the output of user defined functions over a
% periodic orbit
% Input:
%   orb: periodic orbit data structure
%    -> sig: sloution signature (event list)
%    -> U: state variable vector (M*N*n)
%    -> T: segment lengths (N)
%    -> p: parameter vector
%    -> n: number of degrees of freedom
%    -> M: Chebyshev mesh resolution
%   sys: names of the functions that define the system
%    -> tau: time delay and its parameter Jacobians
%    -> tau_no: number of distinct time delays
%   fs: functions for plotting, default: fs1(t,y)=t wrt fs2(t,y)=y
%   res: output interpolation resolution (default 100)

if nargin<4
    res = 100; % default plot resolution
end

% Unpack orbit data
[t0,x0] = bvp2sig(orb.U,orb.T,orb.M);
del0 = po_delay_interp(orb.U,orb.T,orb.p,orb.M,sys);
xd0 = del0.ud;
[t1,x1] = bvp2sig(orb.U,orb.T,orb.M,res);
del1 = po_delay_interp(orb.U,orb.T,orb.p,orb.M,sys,res);
xd1 = del1.ud;

% Plot the state wrt t or the outputs of fs
if nargin<3 || isempty(fs)
    plot(t1,x1); hold on
    plot(t0,x0,'k.');
    xlabel('$t$'); ylabel('$u$');
    for i=1:length(orb.T)+1
        xline(sum(orb.T(1:i-1)),':k');
    end
    xlim('tight');
else
    [px0,py0,xlab,ylab] = feval(fs,t0,x0,xd0,orb.p);
    [px1,py1] = feval(fs,t1,x1,xd1,orb.p);
    plot(px1,py1,'k'); hold on
    plot(px0,py0,'k.');
    xlabel(xlab); ylabel(ylab);
    axis('padded');
end

hold off;
box on;

end

