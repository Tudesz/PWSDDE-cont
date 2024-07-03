function [px,py,xlab,ylab] = nlne_rob_p(t,x,xd,p)
%NLNE_ROB_P Periodic orbit plot functions for nonlinear robotic arm 
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
xlab = '$x(t)$';
ylab = '$\dot{x}(t)$';

end