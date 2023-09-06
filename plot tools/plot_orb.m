function plot_orb(orb,fs,res)
%PLOT_ORB Plot the state or defined functions on the bvp solution
% Input:
%   orb: periodic orbit data structure
%    -> sig: sloution signature (event list)
%    -> U: state variable vector (M*N*n)
%    -> T: segment lengths (N)
%    -> p: parameter vector
%    -> n: number of degrees of freedom
%    -> M: Chebyshev mesh resolution
%   fs: functions for plotting, default: fs1(t,y)=t wrt fs2(t,y)=y
%   res: output interpolation resolution (default 100)

if nargin<3
    res = 100;
end
[t0,x0] = bvp2sig(orb.U,orb.T,orb.M);
[t1,x1] = bvp2sig(orb.U,orb.T,orb.M,res);
if nargin<2 || isempty(fs)
    plot(t1,x1); hold on
    plot(t0,x0,'k.');
    xlabel('$t$'); ylabel('$u$');
    for i=1:length(orb.T)+1
        xline(sum(orb.T(1:i-1)),':k');
    end
else
    [px0,py0] = feval(fs,t0,x0,orb.p);
    [px1,py1] = feval(fs,t1,x1,orb.p);
    plot(px1,py1,'k'); hold on
    plot(px0,py0,'k.');
    xlabel('$f_1(t,x)$'); ylabel('$f_2(t,x)$');
end

hold off;
box on;

end

