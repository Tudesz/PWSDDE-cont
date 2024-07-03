function plot_orb_delay_interp(orb,sys,u_i,res)
%PLOT_ORB_DELAY_INTERP Plot the state and the interpolation of its delayed 
% terms shifted by tau on the same plot for validation
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
%   res: interpolation resolution (default M)
%   u_i: index of state dimension to plot (default all)
%   res: output interpolation resolution (default 100)

if nargin<3
    u_i = 1:orb.n;
end
if nargin<4
    res = 100;
end
[t0,x0] = bvp2sig(orb.U,orb.T,orb.M);
xd0 = po_delay_interp(orb.U,orb.T,orb.p,orb.M,sys);
[t1,x1] = bvp2sig(orb.U,orb.T,orb.M,res);
xd1 = po_delay_interp(orb.U,orb.T,orb.p,orb.M,sys,res);

tau = zeros(1,sys.tau_no);
for i = 1:sys.tau_no
    tau(i) = feval(sys.tau,orb.p,i,2); % time delays
end
T = sum(orb.T);

plot(t1,x1(u_i,:),'b'); hold on
plot(t0,x0(u_i,:),'b.');
for i=1:sys.tau_no
    plot(mod(t1-tau(i),T),squeeze(xd1(u_i,i,:)),'k--'); hold on
    plot(mod(t0-tau(i),T),squeeze(xd0(u_i,i,:)),'k.');
end
xlabel('$t$'); ylabel('$u$');
for i=1:length(orb.T)+1
    xline(sum(orb.T(1:i-1)),':k');
end

hold off;
box on;

end

