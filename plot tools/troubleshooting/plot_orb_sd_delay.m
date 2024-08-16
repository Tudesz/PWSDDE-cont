function plot_orb_sd_delay(orb,sys,res)
%PLOT_ORB_SD_DELAY Plot the evaluations of state dependent delays in a 
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
%   res: output interpolation resolution (default 100)

% Initialization
if nargin < 3
    res = 100;
end

% Evaluate time delays
[t0,x0] = bvp2sig(orb.U,orb.T,orb.M);
tau0 = zeros(length(t0),sys.tau_no);
for j = 1:sys.tau_no
    for i = 1:length(t0)
        tau0(i,j) = feval(sys.tau,x0(:,i),orb.p,j,2);
    end
end
[t1,x1] = bvp2sig(orb.U,orb.T,orb.M,res);
tau1 = zeros(length(t1),sys.tau_no);
for j = 1:sys.tau_no
    for i = 1:length(t1)
        tau1(i,j) = feval(sys.tau,x1(:,i),orb.p,j,2);
    end
end

plot(t1,tau1,'b'); hold on
plot(t0,tau0,'b.');

xlabel('$t$'); ylabel('$\tau_i$');
for i=1:length(orb.T)+1
    xline(sum(orb.T(1:i-1)),':k');
end
hold off;
box on;

end

