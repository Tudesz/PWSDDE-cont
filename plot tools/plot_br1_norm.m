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
%   annotate: level of annotations on bifurcation plot (default 1)
%     -> 0: no extra information
%     -> 1: display stability and bifurcation points
%     -> 2: display tangent vectors and number of steps

if nargin<3
    annotate = 1; % display stability and bifurcation points
end

% Preprocess output data
un = zeros(1,length(branch));
p = zeros(1,length(branch));
tg = zeros(2,length(branch));
stab = zeros(1,length(branch));
i_gr = []; % index of grazing orbits
i_sl = []; % index of sliding orbits
i_va = []; % index of vanishing segments
i_sc = []; % index of stability change points
for i = 1:length(branch)
    un(i) = norm([branch(i).U; branch(i).T]);
    p(i) = branch(i).p(p_i);
    if i<length(branch)
        tg(:,i) = branch(i).tg;
    end
    if max(abs(branch(i).mu))>1+1e-7
        stab(i) = 1;
    end
    if ~isempty(strfind(branch(i).bif_type,'grazing'))
        i_gr(end+1) = i;
    end
    if ~isempty(strfind(branch(i).bif_type,'sliding'))
        i_sl(end+1) = i;
    end
    if ~isempty(strfind(branch(i).bif_type,"Vanishing"))
        i_va(end+1) = i;
    end
    if ~isempty(strfind(branch(i).bif_type,'Stability'))
        i_sc(end+1) = i;
    end
end

% Plot results
switch annotate
    case 0
        plot(p,un);
    case 1
        pi = p; pi(stab==0) = NaN;
        plot(p,un,'b',pi,un,'r'); hold on
        scatter(p(i_gr),un(i_gr),'ob');
        scatter(p(i_sl),un(i_sl),'oy');
        scatter(p(i_sc),un(i_sc),'or');
        scatter(p(i_va),un(i_va),'oc');
    case 2
        plot(p,un,'b'); hold on
        text(p(1:10:end),un(1:10:end),string(0:10:length(un)-1),...
        'VerticalAlignment','bottom','HorizontalAlignment','left')
        quiver(p,un,tg(2,:),tg(1,:),0,'g')
end
xlabel('$\lambda$'); ylabel('$|u|$'); box on;

end

