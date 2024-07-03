function [p,Un,tg,mu_c,ibif] = get_br_data(branch,p_i,u_i)
%GET_BR_DATA Extract indices of bifurcation points from branch data
% Input:
%   branch: continuation run output structure
%    -> bif_type: text data on encountered bifurcations
%   p_i: bifurcation parameter indicies
%   u_i: index of state dimension to return amplitude of (default 0, which
%       returns the norm of the solution vector)
% Output:
%   p: requested bifurcation parameters
%   Un: solution norm or amplitude depending on u_i
%   tg: tangent vectors from the predictor step of the pseudo-archlengt
%       method
%   mu_c: critical floquet multiplier of the orbits
%   ibif: index data structure marking encountered events
%    -> i_bd: boundary points of continuation
%    -> i_gr: indices of grazing points
%    -> i_sl: indices of sliding orbits
%    -> i_va: indices of vanishing segments
%    -> i_sc: indices of stability change points
%    -> i_fp: indices of detected fold points
%    -> i_q:  indices of user defined sign changes

if nargin < 3
    u_i = 0; % return norm of the solution vector by default
end
u_res = 100; % resolution for amplitude approximation

% fill up output variables
p = zeros(length(p_i),length(branch));
Un = zeros(1,length(branch));
mu_c = zeros(1,length(branch));
tg = zeros(length(branch(1).bif_p)+1,length(branch));
ibif = struct('i_bd',[],'i_gr',[],'i_sl',[],'i_va',[],'i_sc',[],...
    'i_fp',[],'i_q',[]);
for i = 1:length(branch)
    p(:,i) = branch(i).p(p_i);
    if u_i == 0
        Un(i) = norm([branch(i).U; branch(i).T]);
    else
        [~,u_temp] = bvp2sig(branch(i).U,branch(i).T,branch(i).M,u_res);
        Un(i) = (max(u_temp(u_i,:))-min(u_temp(u_i,:)))/2;
    end
    if ~isempty(branch(i).tg)
        tg(:,i) = branch(i).tg;
    end
    if ~isempty(strfind(branch(i).bif_type,'Starting')) || ...
            ~isempty(strfind(branch(i).bif_type,'domain'))
        ibif.i_bd(end+1) = i;
    end
    if ~isempty(branch(i).mu_crit)
        mu_c(i) = branch(i).mu_crit;
    else
        mu_c(i) = NaN;
    end
    if ~isempty(strfind(branch(i).bif_type,'grazing'))
        ibif.i_gr(end+1) = i;
    end
    if ~isempty(strfind(branch(i).bif_type,'sliding'))
        ibif.i_sl(end+1) = i;
    end
    if ~isempty(strfind(branch(i).bif_type,'segment'))
        ibif.i_va(end+1) = i;
    end
    if ~isempty(strfind(branch(i).bif_type,'Stability'))
        ibif.i_sc(end+1) = i;
    end
    if ~isempty(strfind(branch(i).bif_type,'Fold'))
        ibif.i_fp(end+1) = i;
    end
    if ~isempty(strfind(branch(i).bif_type,'user'))
        ibif.i_q(end+1) = i;
    end
end

end

