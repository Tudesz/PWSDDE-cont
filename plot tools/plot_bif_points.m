function plot_bif_points(branch,p_i,u_i,p_type,b_type,ms)
%PLOT_BIF_POINTS Mark selected type of bifurcation points on a selected
%type of diagram
% Input:
%   branch: continuation run output structure
%    -> u: solution vectors (M*n)
%    -> T: segment lengths (N)
%    -> p:  bifurcation parameters
%   p_i: bifurcation parameter indeicies
%   u_i: index of state whose amplitudes are displayed (use norm by 
%       default and if applicable)
%   p_type: type of bifurcation plot to mark
%    -> 1: plot_br1_ampl
%    -> 2: plot_br1_norm
%    -> 3: plot_br1_qmon
%    -> 4: plot_br2_3D
%    -> 5: plot_br2_par
%    -> 6: plot_br2_qmon
%   b_type: bifurcation type to mark
%    -> 0: parameter boundary points
%    -> 1: grazing points
%    -> 2: sliding orbits
%    -> 3: vanishing segments
%    -> 4: stability changes
%    -> 5: user defined sign changes
% ms: marker data passed on to the scatter function (default empty)

if nargin<6 
    ms = {}; % no marker style data
elseif ~iscell(ms)
    ms = {ms};
end

switch p_type
    case 1 % plot_br1_ampl
        [p,un,~,~,ibif] = get_br_data(branch,p_i(1),u_i);
    case 2 % plot_br1_norm
        [p,un,~,~,ibif] = get_br_data(branch,p_i(1),0);
    case 3 % plot_br1_qmon
        [p,~,~,~,ibif] = get_br_data(branch,p_i(1));
        q = cell2mat({branch.q});
        un = q(u_i,:);
    case 4 % plot_br2_3D
        [p,un,~,~,ibif] = get_br_data(branch,p_i,u_i);
    case 5 % plot_br2_par
        [pu,~,~,~,ibif] = get_br_data(branch,p_i);
        p = pu(1,:); un = pu(2,:);
    case 6 % plot_br2_qmon
        [p,~,~,~,ibif] = get_br_data(branch,p_i);
        q = cell2mat({branch.q});
        un = q(u_i,:);
end

switch b_type
    case 0 % domain boundary
        ib = ibif.i_bd;
    case 1 % grazing
        ib = ibif.i_gr;
    case 2 % sliding
        ib = ibif.i_sl;
    case 3 % vanishing
        ib = ibif.i_va;
    case 4 % stability
        ib = ibif.i_sc;
    case 5 % monitor
        ib = ibif.i_q;
end

if size(p,1) == 1
    scatter(p(ib),un(ib),ms{:});
else
    scatter3(p(1,ib),p(2,ib),un(ib),ms{:});
end


end