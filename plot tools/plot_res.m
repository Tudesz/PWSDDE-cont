function plot_res(res,ind,p)
%PLOT_RES Visualize the output of sim_ns_dde.m
% Input:
%   res: DDE solver output structure
%   ind: index of state dimension to plot (default all), if a name of a
%   funciton is given plot according to that
%   p: system parameter vector neccessary for the plotting function (only
%   needed if the name of a function is passed in ind)

colors = [0 0.4470 0.7410;
    0.8500 0.3250 0.0980;
    0.9290 0.6940 0.1250;
    0.4940 0.1840 0.5560;
    0.4660 0.6740 0.1880;
    0.3010 0.7450 0.9330;
    0.6350 0.0780 0.1840];
 
for i=1:length(res)
    ci = mod(res(i).mode,size(colors,1));
    if ci == 0
        color = colors(end,:);
    else
        color = colors(ci,:);
    end
    if nargin<2
        plot(res(i).x,res(i).y.','Color',color);
    elseif ischar(ind) && nargin>2
        [px,py,xlab,ylab] = feval(ind,res(i).x,res(i).y,p);
        plot(px,py,'Color',color);
    else
        plot(res(i).x,res(i).y(ind,:),'Color',color);
    end
    if i==1
        hold on
    end
end
if nargin<2 || ~ischar(ind)
   xlabel('$t$');  ylabel('$x$'); 
else
    xlabel(xlab); ylabel(ylab);
end
axis('padded')
box on; hold off

end