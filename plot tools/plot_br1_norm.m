function plot_br1_norm(branch,p_i,annotate)
%PLOT_BR1_NORM Display continuation result norms wrt 1 continuation
%parameter
% Input:
%   branch: continuation run output structure
%    -> u: solution vectors (M*n)
%    -> T: segment lengths (N)
%    -> p: bifurcation parameters
%    -> tg: solution tangents in norm
%   p_i: bifurcation parameter index
%   annotate: level of annotations on the bifurcation plot (default 1)
%     -> 0: no extra information
%     -> 1: display stability and bifurcation points
%     -> 2: display tangent vectors and number of steps (with default p_i!)
%     -> 3: display only bifurcation points
%     -> 4: display only stability

if nargin<3
    annotate = 1; % display stability and bifurcation points
end

% Preprocess output data
[p,un,tg,mu_c,ibif] = get_br_data(branch,p_i(1),0);

% Plot results
switch annotate
    case 0
        plot(p,un);
    case 1
        pi = p; pi(abs(mu_c)<=1) = NaN;
        plot(p,un,'b',pi,un,'r'); hold on
        mark_br_bifs(ibif,[p; un]);
    case 2
        plot(p,un,'b'); hold on
        text(p(1:10:end),un(1:10:end),string(0:10:length(un)-1),...
        'VerticalAlignment','bottom','HorizontalAlignment','left')
        quiver(p,un,tg(2,:),tg(1,:),0,'g')
    case 3
        plot(p,un); hold on
        mark_br_bifs(ibif,[p; un]);
    case 4
        pi = p; pi(abs(mu_c)<=0) = NaN;
        plot(p,un,'b',pi,un,'r');
end
xlabel('$\lambda$'); ylabel('$|u|$'); box on;

end

