function [px,py,xlab,ylab] = fricosc_p(t,x,xd,p)
%FRICOSC_P Periodic orbit plot functions for dry friction oscillator
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

% Sytem parameters
%   p(1): zeta
%   p(2): omega
%   p(3): tau
%   p(4): f_0
%   p(5): eta

% velocity
px = x(2,:);
xlab = '$v$(t)';
% sum of active forces
py = -x(1,:)-2*p(1)*x(2,:)+p(4)*cos(p(2)*(t-p(3)));
ylab = '$F_{\mathrm{act}}$(t)';

end