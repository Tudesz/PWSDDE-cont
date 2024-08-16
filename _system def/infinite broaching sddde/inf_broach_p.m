function [px,py,xlab,ylab] = inf_broach_p(t,x,xd,p)
%INF_BROACH_P Periodic orbit plot functions for an example "infinite" 
% broaching problem
% Input:
%   t: time mesh of periodic orbit
%   x: state vector on t
%   xd: delayed state vector on t
%   p: system parameter vector
% Output:
%   px: x coordinates to be plotted
%   py: y coordinates to be plotted
%   xlab: label on the x axis
%   ylab: label on the y axis

px = x(1,:);
py = x(2,:);
xlab = '$y(t)$';
ylab = '$\dot{y}(t)$';

end