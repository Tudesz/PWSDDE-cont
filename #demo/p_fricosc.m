function [px,py] = p_fricosc(t,x,p)
%P_FRICOSC Periodic orbit plot functions for dry friction oscillator
% Input:
%   t: time mesh of periodic orbit
%   x: state vector on t
%   p: system parameter vector
% Output:
%   px: x coordinates to be plotted
%   py: y coordinates to be plotted

% Sytem parameters
%   p(1): zeta
%   p(2): omega
%   p(3): tau
%   p(4): f_0
%   p(5): eta

% velocity
px = x(2,:);
% sum of active forces
py = -x(1,:)-2*p(1)*x(2,:)+p(4)*cos(p(2)*(t-p(3)));

end