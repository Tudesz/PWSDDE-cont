function plot_br2_par(branch,p_i,annotate)
%PLOT_BR2_PAR  Display continuation results in 2 continuation parameters
% on a 2D plot
% Input:
%   branch: continuation run output structure
%    -> p: bifurcation parameters (lp)
%    -> tg: solution tangents in norm
%   p_i: bifurcation parameter indicies
%   annotate: level of annotations on bifurcation plot (default 1)
%     -> 0: no extra information
%     -> 1: display stability and bifurcation points
%     -> 2: display tangent vectors and number of steps (with default p_i!)
%     -> 3: display only bifurcation points
%     -> 4: display only stability

if nargin<3
    annotate = 1; % display stability and bifurcation points
end

% Preprocess output data
[p,~,tg,mu_c,ibif] = get_br_data(branch,p_i);

% Plot results
switch annotate
    case 0
        plot(p(1,:),p(2,:));
    case 1
        pi = p; pi(:,abs(mu_c)<=1) = NaN;
        plot(p(1,:),p(2,:),'b',pi(1,:),pi(2,:),'r'); hold on
        mark_br_bifs(ibif,p);
    case 2
        plot(p(1,:),p(2,:),'b'); hold on
        text(p(1,1:10:end),p(2,1:10:end),string(0:10:length(p)-1),...
            'VerticalAlignment','bottom','HorizontalAlignment','left')
        quiver(p(1,:),p(2,:),tg(2,:),tg(3,:),0,'g')
    case 3
        plot(p(1,:),p(2,:)); hold on;
        mark_br_bifs(ibif,p);
    case 4
        pi = p; pi(:,abs(mu_c)<=1) = NaN;
        plot(p(1,:),p(2,:),'b',pi(1,:),pi(2,:),'r');
end
xlabel(sprintf('$\\lambda_%i$',p_i(1))); 
ylabel(sprintf('$\\lambda_%i$',p_i(2)));
box on;

end

