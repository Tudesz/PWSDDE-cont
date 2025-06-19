function plot_br1_qmon(branch,p_i,q_i,annotate)
%PLOT_BR1_QMON Display a user defined monitor function value wrt 1 
% continuation parameter
% Input:
%   branch: continuation run output structure
%    -> u: solution vectors (M*n)
%    -> T: segment lengths (N)
%    -> p:  bifurcation parameters
%    -> sig: solution signature (n)
%   p_i: bifurcation parameter index
%   q_i: index of the monitor function to be displayed
%   annotate: level of annotations on the bifurcation plot (default 1)
%     -> 0: no extra information
%     -> 1: display stability and bifurcation points
%     -> 2: display only bifurcation points
%     -> 3: display only stability

if nargin<4
    annotate = 1; % display stability and bifurcation points
end

% Preprocess output data
[p,~,~,mu_c,ibif] = get_br_data(branch,p_i(1));
q = cell2mat({branch.q});
q_mon = q(q_i,:);

% Plot results
switch annotate
    case 0
        plot(p,q_mon);
    case 1
        pi = p; pi(abs(mu_c)<=1) = NaN;
        plot(p,q_mon,'b',pi,q_mon,'r'); hold on
        mark_br_bifs(ibif,[p; q_mon]);
    case 2
        plot(p,q_mon); hold on;
        mark_br_bifs(ibif,[p; q_mon]);
    case 3
        pi = p; pi(abs(mu_c)<=1) = NaN;
        plot(p,q_mon,'b',pi,q_mon,'r');
end
xlabel(sprintf('$\\lambda_%i$',p_i)); 
ylabel(sprintf('$q_%i$',q_i));
box on;

end

