function plot_br1_ampl(branch,p_i,u_i,p_bif)
%PLOT_BR1_AMPL Display continuation result norms wrt 1 continuation
%parameter
% Input:
%   branch: continuation run output structure
%    -> u: solution vectors (M*n)
%    -> T: segment lengths (N)
%    -> p:  bifurcation parameters
%    -> sig: solution signature (n)
%   p_i: bifurcation parameter index
%   u_i: index of state whose maxima are displayed
%   M: segment mesh resolution
%   p_bif: if true mark bifurcation points as well (default true)

if nargin<4
    p_bif = true; %mark bifurcation points
end

% Preprocess output data
ampl = zeros(1,length(branch));
p = zeros(1,length(branch));
stab = zeros(1,length(branch));
i_gr = []; % index of grazing orbits
i_sl = []; % index of sliding orbits
i_va = []; % index of vanishing segments
i_sc = []; % index of stability change points
for i = 1:length(branch)
    [~,u_temp] = bvp2sig(branch(i).U,branch(i).T,branch(i).M);
    ampl(i) = max(u_temp(u_i,:))-min(u_temp(u_i,:));
    p(i) = branch(i).p(p_i);
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
pi = p; pi(stab==0) = NaN;
plot(p,ampl,'b',pi,ampl,'r');
xlabel('$\lambda$'); 
ylabel(sprintf('$\\max(u_%i) - \\min(u_%i)$',u_i,u_i));
if p_bif
    hold on
    scatter(p(i_gr),ampl(i_gr),'ob');
    scatter(p(i_sl),ampl(i_sl),'oy');
    scatter(p(i_sc),ampl(i_sc),'or');
    scatter(p(i_va),ampl(i_va),'oc');
end
box on;

end

