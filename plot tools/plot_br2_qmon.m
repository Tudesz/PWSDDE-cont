function plot_br2_qmon(branch,p_i,q_i,annotate,ls)
%PLOT_BR1_QMON Display a user defined monitor function value wrt 2
% continuation parameters
% Input:
%   branch: continuation run output structure
%    -> u: solution vectors (M*n)
%    -> T: segment lengths (N)
%    -> p:  bifurcation parameters
%    -> sig: solution signature (n)
%   p_i: bifurcation parameter index
%   q_i: index of the monitor function to be displayed
%   annotate: level of annotations on the bifurcation plot (default 1)
%     -> 0: no extra information (ls options are applied)
%     -> 1: display stability and bifurcation points
%     -> 2: display only bifurcation points
%     -> 3: display only stability
%   ls: linestyle data forwarded to the plot function (default empty)

if nargin<4
    annotate = 1; % display stability and bifurcation points
end
if nargin<5 
    ls = {}; % no listyle data
elseif ~iscell(ls)
    ls = {ls};
end

% Preprocess output data
[p,~,~,mu_c,ibif] = get_br_data(branch,p_i);
q = cell2mat({branch.q});
q_mon = q(q_i,:);

% Plot results
switch annotate
    case 0
        plot3(p(1,:),p(2,:),q_mon,ls{:});
    case 1
        pi = p; pi(abs(mu_c)<=1) = NaN;
        plot3(p(1,:),p(2,:),q_mon,'b',pi(1,:),pi(2,:),q_mon,'r'); hold on
        mark_br_bifs(ibif,[p; q_mon]);
    case 2
        plot3(p(1,:),p(2,:),q_mon,'k'); hold on;
        mark_br_bifs(ibif,[p; q_mon]);
    case 3
        pi = p; pi(abs(mu_c)<=1) = NaN;
        plot3(p(1,:),p(2,:),q_mon,'b',pi(1,:),pi(2,:),q_mon,'r');
end
xlabel(sprintf('$\\lambda_%i$',p_i(1))); 
ylabel(sprintf('$\\lambda_%i$',p_i(2))); 
zlabel(sprintf('$q_%i$',q_i));
box on;

end

