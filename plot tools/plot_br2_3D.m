function plot_br2_3D(branch,p_i,u_i,annotate)
%PLOT_BR2_3D  Display continuation result norms or amplitudes wrt 2 
% continuation parameters on a 3D plot
% Input:
%   branch: continuation run output structure
%    -> u: solution vectors (M*n)
%    -> T: segment lengths (N)
%    -> p:  bifurcation parameters
%   p_i: bifurcation parameter indeicies
%   u_i: index of state whose amplitudes are displayed (use norm by default)
%   annotate: level of annotations on bifurcation plot (default 1)
%     -> 0: norms with no extra information
%     -> 1: norms with stability and bifurcation points
%     -> 2: amplitudes with no extra information
%     -> 3: amplitudes with stability and bifurcation points
%     -> 4: norms with only bifurcation points
%     -> 5: norms with only stability information
%     -> 6: amplitudes with only bifurcation points
%     -> 7: amplitudes with only stability information

if nargin<4
    annotate = 1; % display norms, stability and bifurcation points
end

% Preprocess output data
switch annotate
    case {0,1,4,5}
        [p,Un,~,mu_c,ibif] = get_br_data(branch,p_i,0);
    otherwise
        [p,Un,~,mu_c,ibif] = get_br_data(branch,p_i,u_i);
end

% Plot results
switch annotate
    case {0,2}
        plot3(p(1,:),p(2,:),Un);
    case {1,3}
        pi = p; pi(:,abs(mu_c)<=1) = NaN;
        plot3(p(1,:),p(2,:),Un,'b',pi(1,:),pi(2,:),Un,'r'); hold on;
        mark_br_bifs(ibif,[p; Un]);
    case {4,6}
        plot3(p(1,:),p(2,:),Un); hold on;
        mark_br_bifs(ibif,[p; Un]);
    case {5,7}
        pi = p; pi(:,abs(mu_c)<=1) = NaN;
        plot3(p(1,:),p(2,:),Un,'b',pi(1,:),pi(2,:),Un,'r');
end
xlabel('$\lambda_1$'); ylabel('$\lambda_2$');
switch annotate
    case {0,1,4,5}
        zlabel('$|u|$');
    otherwise
        zlabel(sprintf('$(\\max(u_%i) - \\min(u_%i))/2$',u_i,u_i));
end
box on;

end

