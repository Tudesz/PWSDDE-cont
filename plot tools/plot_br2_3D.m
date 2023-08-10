function plot_br2_3D(branch,p_i,u_i,annotate)
%PLOT_BR2_3D Display 2 parameter continuation results in bifurcation
%parameters and solution norms in 3D
% Input:
%   branch: continuation run output structure
%    -> u: solution vectors (M*n)
%    -> T: segment lengths (N)
%    -> p:  bifurcation parameters
%   p_i: bifurcation parameter indeicies
%   u_i: index of state whose maxima are displayed (use norm by default)
%   annotate: level of annotations on bifurcation plot (default 1)
%     -> 0: norms with no extra information
%     -> 1: norms with stability and bifurcation points
%     -> 2: amplitudes with no extra information
%     -> 3: amplitudes with stability and bifurcation points

if nargin<4
    annotate = 1; % display norms, stability and bifurcation points
end

% Preprocess output data
ampl = zeros(1,length(branch));
p = zeros(2,length(branch));
stab = zeros(1,length(branch));
i_gr = []; % index of grazing orbits
i_sl = []; % index of sliding orbits
i_va = []; % index of vanishing segments
i_sc = []; % index of stability change points
for i = 1:length(branch)
    p(:,i) = branch(i).p(p_i);
    if nargin<3 || annotate < 2
        ampl(i) = norm([branch(i).U; branch(i).T]);
    else
        [~,u_temp] = bvp2sig(branch(i).U,branch(i).T,branch(i).M);
        ampl(i) = max(u_temp(u_i,:))-min(u_temp(u_i,:));
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
    case {0,2}
        plot3(p(1,:),p(2,:),ampl);
    case {1,3}
        pi = p; pi(:,stab==0) = NaN;
        plot3(p(1,:),p(2,:),ampl,'b',pi(1,:),pi(2,:),ampl,'r'); hold on;
        scatter3(p(1,i_gr),p(2,i_gr),ampl(i_gr),'ob');
        scatter3(p(1,i_sl),p(2,i_sl),ampl(i_sl),'oy');
        scatter3(p(1,i_sc),p(2,i_sc),ampl(i_sc),'or');
        scatter3(p(1,i_va),p(2,i_va),ampl(i_va),'oc');
end
xlabel('$\lambda_1$'); ylabel('$\lambda_2$');
if nargin<3 || annotate < 2
    zlabel('$|u|$');
else
    zlabel(sprintf('$\\max(u_%i) - \\min(u_%i)$',u_i,u_i));
end
box on;

end

