function plot_orb_discprop(orb,sys)
%PLOT_ORB_DISCPROP Display how discontinuities propagate in periodic orbits
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

% points of discontinuity
tc = cumsum(orb.T).';
T = sum(orb.T);

% plots points of discontinuity and their propagation with tau_l
xline([0 tc],'k--');hold on
tau_label = cell(1,3*sys.tau_no);
[t_mesh,~] = bvp2sig(orb.U,orb.T,orb.M);
em = ones(size(t_mesh));
e = ones(size(tc));
for i = 1:sys.tau_no
    % time delay
    tau = feval(sys.tau,orb.p,i,2);
    % current time mesh
    scatter(t_mesh,i*em,20,'bo');
    scatter([0 tc],i*[1 e],50,'bo','filled');
    tau_label{3*(i-1)+2} = sprintf('$\\tau_%i$ = %0.3f',i,tau);
    % t_mesh - tau
    scatter(mod(t_mesh-tau,T),(i-0.25)*em,20,'kx');
    scatter(mod(tc-tau,T),(i-0.25)*e,100,'kx');
    tau_label{3*(i-1)+1} = sprintf('$t-\\tau_%i$',i);
    % t_mesh + tau
    scatter(mod(t_mesh+tau,T),(i+0.25)*em,20,'kx');
    scatter(mod(tc+tau,T),(i+0.25)*e,100,'kx');
    tau_label{3*(i-1)+3} = sprintf('$t+\\tau_%i$',i);
end
hold off; box on;
xlabel('$t$'); ylabel('$\tau$')
xlim([-0.05*T, 1.05*T]);
ylim([0.5 sys.tau_no+0.5]);
yticks(sort([1:sys.tau_no (1:sys.tau_no)-0.25 (1:sys.tau_no)+0.25]));
yticklabels(tau_label)

end

