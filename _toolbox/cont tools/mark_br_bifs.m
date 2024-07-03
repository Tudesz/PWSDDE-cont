function mark_br_bifs(ibif,pdata)
%MARK_BR_BIFS Mark bifurcation points in plots of continuation results
% Input:
%   ibif: index data structure marking encountered events
%    -> i_gr: indices of grazing points (blue)
%    -> i_sl: indices of sliding orbits (yellow)
%    -> i_va: % indices of vanishing segments (cyan)
%    -> i_sc: % indices of stability change points (red)
%    -> i_q:  % indices of user defined sign changes (magenta)
%   pdata: data plotted on the bifurcation diagram [x, y, (z)]


if size(pdata,1) == 2
    scatter(pdata(1,ibif.i_bd),pdata(2,ibif.i_bd),'ob');
    scatter(pdata(1,ibif.i_gr),pdata(2,ibif.i_gr),'ob','filled');
    scatter(pdata(1,ibif.i_sl),pdata(2,ibif.i_sl),'oy','filled');
    scatter(pdata(1,ibif.i_sc),pdata(2,ibif.i_sc),'or','filled');
    scatter(pdata(1,ibif.i_va),pdata(2,ibif.i_va),'oc','filled');
    scatter(pdata(1,ibif.i_fp),pdata(2,ibif.i_fp),'ok');
    scatter(pdata(1,ibif.i_q),pdata(2,ibif.i_q),'om','filled');
else
    scatter3(pdata(1,ibif.i_gr),pdata(2,ibif.i_gr),...
        pdata(3,ibif.i_gr),'ob');
    scatter3(pdata(1,ibif.i_gr),pdata(2,ibif.i_gr),...
        pdata(3,ibif.i_gr),'ob','filled');
    scatter3(pdata(1,ibif.i_sl),pdata(2,ibif.i_sl),...
        pdata(3,ibif.i_sl),'oy','filled');
    scatter3(pdata(1,ibif.i_sc),pdata(2,ibif.i_sc),...
        pdata(3,ibif.i_sc),'or','filled');
    scatter3(pdata(1,ibif.i_va),pdata(2,ibif.i_va),...
        pdata(3,ibif.i_va),'oc','filled');
    scatter3(pdata(1,ibif.i_fp),pdata(2,ibif.i_fp),...
        pdata(3,ibif.i_fp),'ok');
    scatter3(pdata(1,ibif.i_q),pdata(2,ibif.i_q),...
        pdata(3,ibif.i_q),'om','filled');
end

end


