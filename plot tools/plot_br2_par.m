function plot_br2_par(branch,p_i,annotate)
%PLOT_BR2_PAR Display 2 parameter continuation results in bifurcation
%parameters
% Input:
%   branch: continuation run output structure
%    -> p: bifurcation parameters (lp)
%    -> tg: solution tangents in norm
%   p_i: bifurcation parameter indicies
%   annotate: level of annotations on bifurcation plot (default 1)
%     -> 0: no extra information
%     -> 1: display stability and bifurcation points
%     -> 2: display tangent vectors and number of steps

if nargin<3
    annotate = 1; % display stability and bifurcation points
end

% Preprocess output data
p = zeros(2,length(branch));
tg = zeros(2,length(branch));
stab = zeros(1,length(branch));
i_gr = []; % index of grazing orbits
i_sl = []; % index of sliding orbits
i_va = []; % index of vanishing segments
i_sc = []; % index of stability change points
for i = 1:length(branch)
    p(:,i) = branch(i).p(p_i);
    if i<length(branch)
        tg(:,i) = branch(i).tg(2:end);
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
        plot(p(1,:),p(2,:));
    case 1
        pi = p; pi(:,stab==0) = NaN;
        plot(p(1,:),p(2,:),'b',pi(1,:),pi(2,:),'r'); hold on
        scatter(p(1,i_gr),p(2,i_gr),'ob');
        scatter(p(1,i_sl),p(2,i_sl),'oy');
        scatter(p(1,i_sc),p(2,i_sc),'or');
        scatter(p(1,i_va),p(2,i_va),'oc');
    case 2
        plot(p(1,:),p(2,:),'b'); hold on
        text(p(1,1:10:end),p(2,1:10:end),string(0:10:length(p)-1),...
            'VerticalAlignment','bottom','HorizontalAlignment','left')
        quiver(p(1,:),p(2,:),tg(1,:),tg(2,:),0,'g')
end
xlabel('$\lambda_1$'); ylabel('$\lambda_2$'); box on;

end

